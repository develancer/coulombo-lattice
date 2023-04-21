// Coulombo Ⓒ 2023
// [Computer Physics Communications] Różański & Zieliński:
// Exploiting underlying crystal lattice for efficient computation of Coulomb matrix elements in multi-million atoms nanostructures
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

	std::string outputDir;
	if (pr.hasValue("output-dir", outputDir) && !outputDir.empty()) {
		outputDir += '/';
	}

	timer.start("Reading wavefunctions");

	// read input files
	int inputCount = pr.getArgCount();
	if (inputCount >= std::numeric_limits<short>::max()) {
	    throwfr("too many input files");
	}
	std::vector<std::string> fileNames;
	for (int input = 0; input<inputCount; ++input) {
		const std::string arg = pr.getArg(input);

		const char* basename = arg.c_str();
		auto last_slash = arg.find_last_of('/');
		if (last_slash != std::string::npos) {
			basename += last_slash + 1;
		}
		functions->appendFile(pr.getArg(input));
		fileNames.push_back(basename);
	}

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

	timer.start("Computing all quasi-potentials");

	// perform the actual computation
	ProductCollection products = functions->createSelfProducts();
	auto fileNameIterator = fileNames.begin();
	for (const auto& product : products) {
		product->map(calculator.input, false);
		calculator.prepare();

		std::vector<complex> fullPotentialValues = functions->extractAtomCellValues(calculator.potential());
		if (mpi::root()) {
			std::string outputFilePath = outputDir + "potential-" + *fileNameIterator;
			FILE* file = fopen(outputFilePath.c_str(), "w");
			if (!file) {
				throwfr("cannot open file %s for writing", outputFilePath.c_str());
			}
			for (complex x : fullPotentialValues) {
				fprintf(file, "%.12le\n", x.real());
			}
		}
		++fileNameIterator;
	}
}

int main(int argc, char** argv)
{
	try {
		mpi::init(argc, argv);

		Parser cmd;
		cmd.allowValue("atoms");
		cmd.allowValue("dielectric");
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
					"Coulombo v" VERSION " (c) 2023 [Computer Physics Communications] Rozanski & Zielinski:\n"
					"Exploiting underlying crystal lattice for efficient computation of Coulomb matrix elements in multi-million atoms nanostructures\n"
					"\nUSAGE:\n"
					" potentials [FLAGS/OPTIONS] data files ...\n"
					"\nOPTIONS:\n"
					"  --atoms=PATH           path to *.3d file with atoms' positions\n"
					"  --dielectric=VALUE     dielectric constant, default: 1\n"
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
