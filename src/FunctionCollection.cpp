// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "mpi.hpp"
#include "throwf.hpp"
#include "FunctionCollection.hpp"

//----------------------------------------------------------------------

Dimension FunctionCollection::getPaddedDimension() const
{
	if (!broadcaster) {
		throwfl("broadcaster is not yet initialized");
	}
	return broadcaster->getPaddedDimension();
}

std::shared_ptr<Domain<complex>> FunctionCollection::loadDomainFromFile(const std::string& path)
{
	arma::Cube<complex> input;
	if (mpi::root()) {
		// only the root process actually reads the file
		if (!input.load(path)) {
			// if loading complex cube fails, we try to load is as real-valued
			arma::Cube<double> inputReal;
			if (!inputReal.load(path)) {
				// if this doesn't work either, throw an exception
				throwfr("cannot load data file %s", path.c_str());
			}
			input = arma::conv_to<arma::Cube<complex >>::from(inputReal);
		}
	}
	if (!broadcaster) {
		// Broadcaster is created based on dimensions of the input
		Dimension dimension(input.n_rows, input.n_cols, input.n_slices);
		broadcaster.reset(new Broadcaster(dimension));
	}
	return broadcaster->broadcastData(input);
}

//----------------------------------------------------------------------

void WaveFunctionCollection::appendFile(const std::string& path)
{
	functions.push_back(loadDomainFromFile(path));
}

ProductCollection WaveFunctionCollection::createProducts() const
{
	ProductCollection products;
	int functionCount = functions.size();
	for (int fL = 0; fL<functionCount; ++fL) {
		for (int fR = 0; fR<=fL; ++fR) {
			std::shared_ptr<Product> product(
					new ProductFromWavefunctions(
							functions[fL].get(),
							functions[fR].get()
					)
			);
			products.push_back(product);
		}
	}
	return products;
}

//----------------------------------------------------------------------

void SpinFunctionCollection::appendFile(const std::string& path)
{
	(functionsU.size() < functionsD.size() ? functionsU : functionsD).push_back(loadDomainFromFile(path));
}

ProductCollection SpinFunctionCollection::createProducts() const
{
	ProductCollection products;
	if (functionsD.size() != functionsU.size()) {
		throw std::runtime_error("--spin requires even number of data files");
	}
	int functionCount = functionsU.size();
	for (int fL = 0; fL<functionCount; ++fL) {
		for (int fR = 0; fR<=fL; ++fR) {
			std::shared_ptr<Product> product(
					new ProductFromSpinfunctions(
							functionsD[fL].get(), functionsU[fL].get(),
							functionsD[fR].get(), functionsU[fR].get()
					)
			);
			products.push_back(product);
		}
	}
	return products;
}

//----------------------------------------------------------------------
