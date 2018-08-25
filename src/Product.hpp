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

class ProductFromWavefunctions: public Product {
public:
	ProductFromWavefunctions(const Domain<complex>* left, const Domain<complex>* right);

	void map(Domain<complex>& F, bool conjugate) const;

private:
	const Domain<complex>* left;
	const Domain<complex>* right;
};

//----------------------------------------------------------------------

class ProductFromSpinfunctions: public Product {
public:
	ProductFromSpinfunctions(const Domain<complex>* leftD, const Domain<complex>* leftU, const Domain<complex>* rightD,
			const Domain<complex>* rightU);

	void map(Domain<complex>& F, bool conjugate) const;

private:
	const Domain<complex>* leftU;
	const Domain<complex>* leftD;
	const Domain<complex>* rightU;
	const Domain<complex>* rightD;
};

//----------------------------------------------------------------------

#endif
