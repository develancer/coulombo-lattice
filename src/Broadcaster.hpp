// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_BROADCASTER_HPP
#define COULOMBO_BROADCASTER_HPP

#include <memory>
#include "base.hpp"
#include "mpi.hpp"
#include "Dimension.hpp"
#include "Domain.hpp"

//----------------------------------------------------------------------

class Broadcaster {
public:
	Broadcaster(const Dimension& dimension);

	std::shared_ptr<Domain<complex>> broadcastData(const arma::Cube<complex>& input);

	Dimension getPaddedDimension() const;

private:
	Dimension dimensionRaw, dimensionPad;
	DistributedDimension dimensionLocal;

	mpi::Type typeSlicePad;
	mpi::Type typeSliceRaw;
	std::vector<int> zOffsets, zLengths;

	static mpi::Type createTypeDimension(void);

	static Dimension round4FFT(const Dimension& dimension);
};

//----------------------------------------------------------------------

#endif
