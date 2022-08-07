// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fftw3-mpi.h>
#include <list>
#include <map>
#include <string>
#include <utility>
#include "base.hpp"
#include "mpi.hpp"
#include "throwf.hpp"
#include "Broadcaster.hpp"
#include "Domain.hpp"
#include "CoulombCalculator.hpp"
#include "FunctionCollection.hpp"
#include "Graph.hpp"
#include "Parser.hpp"
#include "Pattern.hpp"
#include "Planner.hpp"
#include "Product.hpp"
#include "Timer.hpp"

#define VERSION "2.0"

//----------------------------------------------------------------------

void coulombo(const ParseResults& pr)
{
	Timer timer;
	std::string value;

	// initialize FFTW
	#ifdef _OPENMP
	int threadsPerNode = 1;
	if (pr.hasValue("threads-per-node", value)) {
		threadsPerNode = atoi(value.c_str());
		if (threadsPerNode <= 0) {
			throwfr("invalid value for threads-per-node");
		}
	}
	omp_set_num_threads(threadsPerNode);
	if (threadsPerNode > 1) {
		fftw_init_threads();
	}
	#endif
	fftw_mpi_init();
	#ifdef _OPENMP
	if (threadsPerNode > 1) {
		fftw_plan_with_nthreads(threadsPerNode);
	}
	#endif

	timer.start("Reading atom positions");

	// initialize input type
	int orbitalCount = 20;
	if (pr.hasValue("orbitals", value)) {
		orbitalCount = atoi(value.c_str());
		if (orbitalCount <= 0) {
			throwfr("invalid value for orbitals");
		}
	}
	int headerLinesToSkip = 0;
	if (pr.hasValue("skip-lines", value)) {
		headerLinesToSkip = atoi(value.c_str());
		if (headerLinesToSkip <= 0) {
			throwfr("invalid value for skip-lines");
		}
	}
	std::shared_ptr<FunctionCollection> functions;
	functions.reset(new FunctionCollection(pr.getValue("atoms"), orbitalCount, headerLinesToSkip));

	// parse various settings
	double dielectric = 1.0;
	if (pr.hasValue("dielectric", value)) {
		dielectric = atof(value.c_str());
	}
	std::string integrals = "****";
	pr.hasValue("integrals", integrals);

	std::string outputDir;
	if (pr.hasValue("output-dir", outputDir) && !outputDir.empty()) {
		outputDir += '/';
	}

	timer.start("Reading wavefunctions");

	// read input files
	int inputCount = pr.getArgCount(), hoStateCount = 0, elStateCount = 0;
	if (inputCount >= std::numeric_limits<short>::max()) {
	    throwfr("too many input files");
	}
	for (int input = 0; input<inputCount; ++input) {
		const std::string arg = pr.getArg(input);

		const char* basename = arg.c_str();
		auto last_slash = arg.find_last_of('/');
		if (last_slash != std::string::npos) {
			basename += last_slash + 1;
		}
		if (basename[0] == 'h') {
			if (elStateCount) {
				throwfr("hole states must appear before electron states");
			}
			++hoStateCount;
		}
		else if (basename[0] == 'e') {
			++elStateCount;
		}
		else {
			throwfr("invalid state file name: %s", arg.c_str());
		}
		functions->appendFile(pr.getArg(input));
	}

	timer.start("Preparing plan");

	Pattern pattern(integrals, hoStateCount);

	// generate list of product generators
	auto products = functions->createProducts();

	// generate list of integrals to be computed
	std::shared_ptr<Planner> planner;
	if (mpi::root()) {
		planner = std::make_shared<MasterPlanner>(products.size());
	} else {
		planner = std::make_shared<Planner>();
	}
	int integralCount = 0;
	std::vector<std::array<short, 4>> integralSpecs;
	if (mpi::root())
	for (short i1 = 1; i1<=inputCount; ++i1) {
		for (short i2 = 1; i2<=inputCount; ++i2) {
			for (short i3 = 1; i3<=inputCount; ++i3) {
				for (short i4 = 1; i4<=inputCount; ++i4) {
					if (pattern.match(i1, i2, i3, i4)) {
						planner->addIntegral(i1, i2, i3, i4);
						integralSpecs.push_back({i1, i2, i3, i4});
						++integralCount;
					}
				}
			}
		}
	}

	// compute the correct calculation plan
	planner->computePlan();

	Dimension dimension = functions->getPaddedDimension();

	timer.start("Initializing calculator");

	// get onsite value
	double onsite = 0.0;
	if (pr.hasValue("onsite", value)) {
		onsite = atof(value.c_str());
	}

	// initialize dielectric screening model
	Vector3D<double> stepXYZ = functions->getStepValues();
	std::unique_ptr<Interaction> interaction;
	if (pr.hasValue("tf-lattice", value)) {
		double latticeConstant = atof(value.c_str());
		if (!std::isgreater(latticeConstant, 0.0)) {
			throwfr("--tf-lattice parameter must be positive");
		}
		if (!std::isgreater(dielectric, 1.0)) {
			throwfr("dielectric constant must be >1 to use Thomas-Fermi model");
		}
		interaction.reset( new InteractionThomasFermi(stepXYZ, onsite, dielectric, latticeConstant) );
	} else {
		interaction.reset( new InteractionSimple(stepXYZ, onsite, dielectric) );
	}

	// initialize buffers
	CoulombCalculator calculator(dimension);
	calculator.initialize(*interaction);

	timer.start("Computing requested integrals");

	// perform the actual computation
	complex valueLast;
	std::vector<complex> integralValues(integralCount);
	int lastLeftProduct = -1;
	int lastRightProduct = -1;
	unsigned char lastRightConjugate = 0;

	PlannerStep step;
	while (planner->getNextStep(step)) {
		if (step.left.product!=lastLeftProduct) {
			products[step.left.product]->map(calculator.input, false);
			calculator.prepare();
			lastLeftProduct = step.left.product;
		}
		bool rightConjugate = (step.left.conjugate!=step.right.conjugate);
		if (step.right.product!=lastRightProduct || rightConjugate!=lastRightConjugate) {
			products[step.right.product]->map(calculator.input, rightConjugate);
			lastRightProduct = step.right.product;
			lastRightConjugate = rightConjugate;
			valueLast = calculator.calculate();
		}
		if (mpi::root()) {
			integralValues[step.id] = step.left.conjugate ? std::conj(valueLast) : valueLast;
		}
	}

	timer.start("Exporting results");

	// export results to output files
	if (mpi::root()) {
		std::map<std::array<short, 4>, complex> values;
		for (int i = 0; i<integralCount; ++i) {
			const auto& specs = integralSpecs[i];
			const complex& value = integralValues[i];
			values[specs] = value;
		}

		char core[5]; core[4] = 0;
		for (int ti=0; ti<2; ++ti)
			for (int tj=0; tj<2; ++tj)
				for (int tk=0; tk<2; ++tk)
					for (int tl=0; tl<2; ++tl) {
						FILE* file = nullptr;
						core[0] = ti ? 'e' : 'h';
						core[1] = tj ? 'e' : 'h';
						core[2] = tk ? 'e' : 'h';
						core[3] = tl ? 'e' : 'h';

						unsigned Ni = ti ? elStateCount : hoStateCount;
						unsigned Nj = tj ? elStateCount : hoStateCount;
						unsigned Nk = tk ? elStateCount : hoStateCount;
						unsigned Nl = tl ? elStateCount : hoStateCount;

						for (unsigned ni=1; ni<=Ni; ++ni)
							for (unsigned nj=1; nj<=Nj; ++nj)
								for (unsigned nk=1; nk<=Nk; ++nk)
									for (unsigned nl=1; nl<=Nl; ++nl) {

										std::array<short, 4> specs;
										specs[0] = ti ? hoStateCount + ni : hoStateCount + 1 - ni;
										specs[1] = tj ? hoStateCount + nj : hoStateCount + 1 - nj;
										specs[2] = tk ? hoStateCount + nk : hoStateCount + 1 - nk;
										specs[3] = tl ? hoStateCount + nl : hoStateCount + 1 - nl;

										auto it = values.find(specs);
										if (it != values.end()) {
											if (!file) file = fopen((outputDir + core + ".txt").c_str(), "w");
											const complex value = it->second;
											fprintf(file, "%2d %2d %2d %2d   %17.14f %17.14f\n",
												ni, nj, nk, nl, std::real(value), std::imag(value)
											);
										}
									}
						if (file) {
							fclose(file);
						}
					}
	}
}

int main(int argc, char** argv)
{
	try {
		mpi::init(argc, argv);

		Parser cmd;
		cmd.allowValue("atoms");
		cmd.allowValue("dielectric");
		cmd.allowValue("integrals");
		cmd.allowValue("onsite");
		cmd.allowValue("orbitals");
		cmd.allowValue("output-dir");
		cmd.allowValue("skip-lines");
		#ifdef _OPENMP
		cmd.allowValue("threads-per-node");
		#endif
		cmd.allowValue("tf-lattice");

		ParseResults pr = cmd.process(argc, argv);
		if (!pr.getArgCount()) {
			fprintf(stderr,
					"Coulombo v" VERSION " (c) 2022 [Computer Physics Communications] Rozanski & Zielinski:\n"
					"Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures\n"
					"\nUSAGE:\n"
					" coulombo [FLAGS/OPTIONS] data files ...\n"
					"\nOPTIONS:\n"
					"  --atoms=PATH           path to *.3d file with atoms' positions\n"
					"  --dielectric=VALUE     dielectric constant, default: 1\n"
					"  --integrals=LIST       comma-separated list of integrals to be computed\n"
					"                           (eg. \"eeee,hhhh,ehhe,eheh\"),\n"
					"                           default: all integrals are computed\n"
					"  --onsite=ENERGY        energy for on-site contribution, default: 0 (eV)\n"
					"  --orbitals=N           number of (spin-)orbitals per atom, default: 20\n"
					"  --output-dir=DIR       directory for output files, default: current\n"
					"  --skip-lines=N         number of lines to be skipped on top of each LCAO file, default: 0\n"
					#ifdef _OPENMP
					"  --threads-per-node=N   number of OpenMP threads per node, default: 1\n"
					#endif
					"  --tf-lattice=VALUE     lattice constant (Å) for Thomas-Fermi-Resta model,\n"
					"                           default: model not applied\n"
					"\n");
			exit(EXIT_FAILURE);
		}
		coulombo(pr);

		mpi::finalize();
	} catch (std::exception& ex) {
		fprintf(stderr, "ERROR: %s\n", ex.what());
		exit(EXIT_FAILURE);
	}
}
