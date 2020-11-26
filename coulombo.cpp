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

#define VERSION "1.0"

//----------------------------------------------------------------------

void coulombo(const ParseResults& pr)
{
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

	// initialize input type
	std::shared_ptr<FunctionCollection> functions;
	if (pr.hasFlag("spin")) {
		functions.reset(new SpinFunctionCollection());
	} else {
		functions.reset(new WaveFunctionCollection());
	}

	// parse various settings
	double gridStep = atof(pr.getValue("step").c_str());
	if (!std::isfinite(gridStep) || gridStep<=0) {
		throwfr("invalid step value");
	}
	double dielectric = 1.0;
	if (pr.hasValue("dielectric", value)) {
		dielectric = atof(value.c_str());
	}
	std::string integrals = "****";
	pr.hasValue("integrals", integrals);
	Pattern pattern(integrals);

	// read input files
	int inputCount = pr.getArgCount();
	if (inputCount >= std::numeric_limits<short>::max()) {
	    throwfr("too many input files");
	}
	for (int input = 0; input<inputCount; ++input) {
		functions->appendFile(pr.getArg(input));
	}

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

	// initialize dielectric screening model
	Vector3D<double> stepXYZ = {gridStep, gridStep, gridStep};
	std::unique_ptr<Interaction> interaction;
	if (pr.hasValue("tf-lattice", value)) {
		double latticeConstant = atof(value.c_str());
		if (!std::isgreater(latticeConstant, 0.0)) {
			throwfr("--tf-lattice parameter must be positive");
		}
		if (!std::isgreater(dielectric, 1.0)) {
			throwfr("dielectric constant must be >1 to use Thomas-Fermi model");
		}
		interaction.reset( new InteractionThomasFermi(stepXYZ, dielectric, latticeConstant) );
	} else {
		interaction.reset( new InteractionSimple(stepXYZ, dielectric) );
	}

	// initialize buffers
	CoulombCalculator calculator(dimension);
	calculator.initialize(*interaction);

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

	// export results to standard output
	if (mpi::root()) {
		for (int i = 0; i<integralCount; ++i) {
			const auto& specs = integralSpecs[i];
			const complex& value = integralValues[i];
			printf("%2d %2d %2d %2d  %12.9lf %12.9lf\n", specs[0], specs[1], specs[2], specs[3], std::real(value), std::imag(value));
		}
	}
}

int main(int argc, char** argv)
{
	try {
		mpi::init(argc, argv);

		Parser cmd;
		cmd.allowFlag("spin");
		cmd.allowValue("dielectric");
		cmd.allowValue("integrals");
		cmd.allowValue("step");
		#ifdef _OPENMP
		cmd.allowValue("threads-per-node");
		#endif
		cmd.allowValue("tf-lattice");

		ParseResults pr = cmd.process(argc, argv);
		if (!pr.getArgCount()) {
			fprintf(stderr,
					"Coulombo v" VERSION " (c) 2018 [Computer Physics Communications] Rozanski & Zielinski:\n"
					"Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures\n"
					"\nUSAGE:\n"
					" coulombo [FLAGS/OPTIONS] data files ...\n"
					"\nFLAGS:\n"
					"  --spin   if each wavefunction consists of two data files: for spin-up and down\n"
					"\nOPTIONS:\n"
					"  --dielectric=VALUE     dielectric constant, default: 1\n"
					"  --integrals=LIST       comma-separated list of integrals to be computed\n"
					"                           (eg. \"ijji,1ij1\"),\n"
					"                           default: all integrals are computed\n"
					"  --step=VALUE           grid step length (Å)\n"
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
