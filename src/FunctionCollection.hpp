// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_FUNCTIONCOLLECTION_HPP
#define COULOMBO_FUNCTIONCOLLECTION_HPP

#include <memory>
#include <string>
#include <vector>
#include "Broadcaster.hpp"
#include "Dimension.hpp"
#include "Domain.hpp"
#include "Product.hpp"

using ProductCollection = std::vector<std::shared_ptr<Product>>;

//----------------------------------------------------------------------

class FunctionCollection {
public:
	FunctionCollection(const std::string& atomPositionsPath, int orbitalCount, int headerLinesToSkip);

	void appendFile(const std::string& path);

	ProductCollection createProducts() const;

	Dimension getPaddedDimension() const;

	Vector3D<double> getStepValues() const;

protected:
	arma::Mat<double> loadAtomsPositions(const std::string& path);

	std::shared_ptr<arma::Mat<complex>> loadFunctionFromFile(const std::string& path);

private:
	int headerLinesToSkip;
	int orbitalCount;
	int totalAtomCount;
	std::vector<int> cellIndices;
	std::unique_ptr<Broadcaster> broadcaster;
	std::vector<std::shared_ptr<arma::Mat<complex>>> functions;
};

//----------------------------------------------------------------------
#endif
