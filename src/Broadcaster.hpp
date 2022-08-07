// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
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
	Broadcaster(const arma::Mat<double>& atomPositions, int orbitalCount);

	std::shared_ptr<arma::Mat<complex>> broadcastData(const arma::Mat<complex>& input);

	std::vector<int> getLocalCellIndices() const;

	Dimension getPaddedDimension() const;

	Vector3D<double> getStepValues() const;

private:
	int orbitalCount;
	Dimension dimensionRaw, dimensionPad;
	DistributedDimension dimensionLocal;

	mpi::Type typeAtomCoeffs;
	std::vector<int> atomOffsets, atomCounts;
	std::vector<int> localCellIndices;
	Vector3D<double> step;

	static mpi::Type createTypeDimension(void);

	static mpi::Type createTypeVector3D(void);

	static Dimension round4FFT(const Dimension& dimension);
};

//----------------------------------------------------------------------

#endif
