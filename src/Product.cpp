// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "Product.hpp"

//----------------------------------------------------------------------

Product::~Product() { }

//----------------------------------------------------------------------

ProductFromTightBinding::ProductFromTightBinding(
		const arma::Mat<complex>* left, const arma::Mat<complex>* right, const std::vector<int>* atomIndices)
		:left(left), right(right), atomIndices(atomIndices) { }

void ProductFromTightBinding::map(Domain<complex>& F, bool conjugate) const
{
	F.zeros();
	const int atomCount = atomIndices->size();
	if (conjugate) {
		for (int i = 0; i < atomCount; ++i) {
			F[(*atomIndices)[i]] = arma::accu( left->col(i) % arma::conj(right->col(i)) );
		}
	}
	else {
		for (int i = 0; i < atomCount; ++i) {
			F[(*atomIndices)[i]] = arma::accu( arma::conj(left->col(i)) % right->col(i) );
		}
	}
}

//----------------------------------------------------------------------
