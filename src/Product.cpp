// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "Product.hpp"

//----------------------------------------------------------------------

Product::~Product() { }

//----------------------------------------------------------------------

ProductFromWavefunctions::ProductFromWavefunctions(const Domain<complex>* left, const Domain<complex>* right)
		:left(left), right(right) { }

void ProductFromWavefunctions::map(Domain<complex>& F, bool conjugate) const
{
	if (conjugate) {
		F = (*left)%arma::conj(*right);
	}
	else {
		F = arma::conj(*left)%(*right);
	}
}

//----------------------------------------------------------------------

ProductFromSpinfunctions::ProductFromSpinfunctions(
		const Domain<complex>* leftD, const Domain<complex>* leftU,
		const Domain<complex>* rightD, const Domain<complex>* rightU)
		:leftU(leftU), leftD(leftD), rightU(rightU), rightD(rightD) { }

void ProductFromSpinfunctions::map(Domain<complex>& F, bool conjugate) const
{
	if (conjugate) {
		F = (*leftU)%arma::conj(*rightU)+(*leftD)%arma::conj(*rightD);
	}
	else {
		F = arma::conj(*leftU)%(*rightU)+arma::conj(*leftD)%(*rightD);
	}
}

//----------------------------------------------------------------------
