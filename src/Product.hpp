// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_PRODUCT_HPP
#define COULOMBO_PRODUCT_HPP

#include "Domain.hpp"

//----------------------------------------------------------------------

class Product {
public:
	virtual void map(Domain<complex>& F, bool conjugate) const =0;
	virtual ~Product();
};

//----------------------------------------------------------------------

class ProductFromTightBinding: public Product {
public:
	ProductFromTightBinding(const arma::Mat<complex>* left, const arma::Mat<complex>* right,
			const std::vector<int>* atomIndices);

	void map(Domain<complex>& F, bool conjugate) const;

private:
	const arma::Mat<complex>* left;
	const arma::Mat<complex>* right;
	const std::vector<int>* atomIndices;
};

//----------------------------------------------------------------------

#endif
