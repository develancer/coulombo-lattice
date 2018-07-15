// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
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
	virtual void appendFile(const std::string& path) =0;
	virtual ProductCollection createProducts() const =0;

	Dimension getPaddedDimension() const;

protected:
	std::shared_ptr<Domain<complex>> loadDomainFromFile(const std::string& path);

private:
	std::unique_ptr<Broadcaster> broadcaster;
};

//----------------------------------------------------------------------

class WaveFunctionCollection: public FunctionCollection {

	void appendFile(const std::string& path);
	ProductCollection createProducts() const;

private:
	std::vector<std::shared_ptr<Domain<complex>>> functions;
};

//----------------------------------------------------------------------

class SpinFunctionCollection: public FunctionCollection {

	void appendFile(const std::string& path);
	ProductCollection createProducts() const;

private:
	std::vector<std::shared_ptr<Domain<complex>>> functionsU;
	std::vector<std::shared_ptr<Domain<complex>>> functionsD;
};

//----------------------------------------------------------------------
#endif
