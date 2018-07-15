// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#include "throwf.hpp"
#include "Broadcaster.hpp"
#include "Round.hpp"

//----------------------------------------------------------------------

Broadcaster::Broadcaster(const Dimension& dimension)
		:zOffsets(mpi::size()), zLengths(mpi::size())
{
	if (mpi::root()) {
		// root process
		dimensionRaw = dimension;
		dimensionPad = round4FFT(dimensionRaw);
	}
	mpi::Type typeDimension = createTypeDimension();
	mpi::broadcast(&dimensionRaw, 1, typeDimension);
	mpi::broadcast(&dimensionPad, 1, typeDimension);
	dimensionLocal = DualDimension(dimensionPad).real;

	int sizePad[2] = {dimensionPad.x, dimensionPad.y};
	int sizeRaw[2] = {dimensionRaw.x, dimensionRaw.y};
	int offsets[2] = {0, 0};
	typeSlicePad = mpi::subarray(2, sizePad, sizeRaw, offsets, MPI_DOUBLE_COMPLEX);
	typeSliceRaw = mpi::subarray(2, sizeRaw, sizeRaw, offsets, MPI_DOUBLE_COMPLEX);

	mpi::allgather(&dimensionLocal.z, 1, MPI_INT, zLengths.data(), 1);
	mpi::allgather(&dimensionLocal.zOffset, 1, MPI_INT, zOffsets.data(), 1);
	const int mpi_size = mpi::size();
	for (int i = 0; i<mpi_size; ++i) {
		zLengths[i] = std::max(0, std::min(zLengths[i], dimensionRaw.z-zOffsets[i]));
	}
}

std::shared_ptr<Domain<complex>> Broadcaster::broadcastData(const arma::Cube<complex>& input)
{
	if (mpi::root()) {
		// root process
		Dimension dimensionInput(input.n_rows, input.n_cols, input.n_slices);
		if (dimensionInput!=dimensionRaw) {
			throwfr("functions' dimensions differ");
		}
	}
	std::shared_ptr<Domain<complex>> domain = std::make_shared<SingleDomain<complex >> (dimensionLocal);
	int sliceCount = zLengths[mpi::rank()];
	MPI_Scatterv(const_cast<complex*>(input.memptr()), zLengths.data(), zOffsets.data(), typeSliceRaw,
			domain->memptr(), sliceCount, typeSlicePad, 0, MPI_COMM_WORLD);
	return domain;
}

Dimension Broadcaster::getPaddedDimension() const
{
	return dimensionPad;
}

mpi::Type Broadcaster::createTypeDimension(void)
{
	int dimensionBlockLengths[3] = {1, 1, 1};
	MPI_Aint dimensionDisplacements[3] = {offsetof(Dimension, x), offsetof(Dimension, y), offsetof(Dimension, z)};
	MPI_Datatype dimensionTypes[3] = {MPI_INT, MPI_INT, MPI_INT};
	return mpi::structure(3, dimensionBlockLengths, dimensionDisplacements, dimensionTypes);
}

Dimension Broadcaster::round4FFT(const Dimension& dimension)
{
	return Dimension(
			Round::up<2, 3, 5>(dimension.x),
			Round::up<2, 3, 5>(dimension.y),
			Round::up<2, 3, 5>(dimension.z)
	);
}

//----------------------------------------------------------------------
