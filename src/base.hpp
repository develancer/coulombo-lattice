// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_BASE_HPP
#define COULOMBO_BASE_HPP

#include <cassert>
#include <complex>

#ifdef NDEBUG
#define ARMA_NO_DEBUG
#endif
#define ARMA_64BIT_WORD
#include <armadillo>
using arma::uword;

using real = double;
using complex = std::complex<real>;

template<typename T>
struct Vector3D {
	T x, y, z;
};

/**
 * e²/4πε₀ [eV·Å]
 */
const double E2_4PE0 = 14.39963737103201;

#endif
