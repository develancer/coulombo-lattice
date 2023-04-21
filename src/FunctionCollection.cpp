// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "mpi.hpp"
#include "throwf.hpp"
#include "FunctionCollection.hpp"
#include <fstream>

//----------------------------------------------------------------------

FunctionCollection::FunctionCollection(const std::string& atomPositionsPath, int orbitalCount, int headerLinesToSkip)
		: headerLinesToSkip(headerLinesToSkip), orbitalCount(orbitalCount), totalAtomCount(0)
{
	arma::Mat<double> atomPositions;
	if (mpi::root()) {
		atomPositions = loadAtomsPositions(atomPositionsPath);
	}
	broadcaster.reset(new Broadcaster(atomPositions, orbitalCount));
	cellIndices = broadcaster->getLocalCellIndices();
}

Dimension FunctionCollection::getPaddedDimension() const
{
	if (!broadcaster) {
		throwfl("broadcaster is not initialized correctly");
	}
	return broadcaster->getPaddedDimension();
}

Vector3D<double> FunctionCollection::getStepValues() const
{
	if (!broadcaster) {
		throwfl("broadcaster is not initialized correctly");
	}
	return broadcaster->getStepValues();
}

arma::Mat<double> FunctionCollection::loadAtomsPositions(const std::string& path)
{
	std::vector<double> coordinates;
	double x, y, z;
	std::ifstream file(path);
	if (!file.is_open()) {
		throwfr("could not open file: %s", path.c_str());
	}
	std::string line;
	while (std::getline(file, line) && sscanf(line.c_str(), "%lf%lf%lf", &x, &y, &z) == 3) {
		++totalAtomCount;
		coordinates.push_back(x);
		coordinates.push_back(y);
		coordinates.push_back(z);
	}
	return arma::Mat<double>(coordinates.data(), 3, totalAtomCount);
}

std::shared_ptr<arma::Mat<complex>> FunctionCollection::loadFunctionFromFile(const std::string& path)
{
	arma::Mat<complex> input;
	if (!broadcaster) {
		throwfl("broadcaster is not initialized correctly");
	}
	if (mpi::root()) {
		// only the root process actually reads the file
		input.set_size(orbitalCount, totalAtomCount);
		std::ifstream file(path);
		if (!file.is_open()) {
			throwfr("could not open file: %s", path.c_str());
		}
		double re, im;
		std::string line;
		for (int i=0; i<headerLinesToSkip; ++i) {
			std::getline(file, line);
		}
		for (int a=0; a<totalAtomCount; ++a) {
			for (int orbital=0; orbital<orbitalCount; ++orbital) {
				if (!std::getline(file, line) || sscanf(line.c_str(), "%lf%lf", &re, &im) != 2) {
					throwfr("file is truncated: %s", path.c_str());
				}
				input(orbital, a) = complex(re, im);
			}
		}
	}
	return broadcaster->broadcastData(input);
}

void FunctionCollection::appendFile(const std::string& path)
{
	functions.push_back(loadFunctionFromFile(path));
}

ProductCollection FunctionCollection::createProducts() const
{
	ProductCollection products;
	int functionCount = functions.size();
	for (int fL = 0; fL<functionCount; ++fL) {
		for (int fR = 0; fR<=fL; ++fR) {
			products.push_back(std::make_shared<ProductFromTightBinding>(
				functions[fL].get(),
				functions[fR].get(),
				&cellIndices
			));
		}
	}
	return products;
}


ProductCollection FunctionCollection::createSelfProducts() const
{
	ProductCollection products;
	int functionCount = functions.size();
	for (int fL = 0; fL<functionCount; ++fL) {
		products.push_back(std::make_shared<ProductFromTightBinding>(
			functions[fL].get(),
			functions[fL].get(),
			&cellIndices
		));
	}
	return products;
}

const std::vector<complex> FunctionCollection::extractAtomCellValues(const SingleDomain<complex>& data) const
{
	std::vector<complex> localValues, rootValues;
	localValues.reserve(cellIndices.size());
	for (int index : cellIndices) {
		localValues.push_back(data[index]);
	}
	broadcaster->gatherArrayAtomWise(localValues, rootValues, MPI_DOUBLE_COMPLEX);
	return rootValues;
}

//----------------------------------------------------------------------
