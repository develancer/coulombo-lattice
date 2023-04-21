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

	template<typename T>
	void scatterArrayAtomWise(const std::vector<T>& rootInput, std::vector<T>& localOutput, MPI_Datatype mpiType) const
	{
		const int localAtomCount = atomCounts[mpi::rank()];
		localOutput.resize(localAtomCount);
		MPI_Scatterv(rootInput.data(), atomCounts.data(), atomOffsets.data(), mpiType,
				localOutput.data(), localAtomCount, mpiType, 0, MPI_COMM_WORLD);
	}

	template<typename T>
	void gatherArrayAtomWise(const std::vector<T>& localInput, std::vector<T>& rootOutput, MPI_Datatype mpiType) const
	{
		const int localAtomCount = atomCounts[mpi::rank()];
		rootOutput.resize(totalAtomCount); // totalAtomCount will be 0 in each non-root process
		MPI_Gatherv(localInput.data(), localAtomCount, mpiType,
				rootOutput.data(), atomCounts.data(), atomOffsets.data(), mpiType, 0, MPI_COMM_WORLD);
	}

private:
	int orbitalCount, totalAtomCount;
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
