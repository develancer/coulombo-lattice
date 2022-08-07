// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "throwf.hpp"
#include "Broadcaster.hpp"
#include "Round.hpp"

static double gcd(double a, double b)
{
	if (a < b) std::swap(a, b);
	while (true) {
		double r = fmod(a, b);
		if (r < 1.0e-10) {
			return b;
		} else if (r < 1.0e-2) {
			return 0.0;
		} else {
			a = b;
			b = r;
		}
	}
}

static void computeGrid(const arma::Mat<double>& atomPositions,
		int dimensionIndex, int& size, double& origin, double& step)
{
	size_t atomCount = atomPositions.n_cols;
	std::vector<double> buffer(atomCount);
	for (size_t a=0; a<atomCount; ++a) {
		buffer[a] = atomPositions(dimensionIndex, a);
	}
	std::sort(buffer.begin(), buffer.end());
	step = NAN;
	for (size_t a=1; a<atomCount; ++a) {
		double next = buffer[a] - buffer[a-1];
		if (next > 1.0e-10) {
			step = std::isnan(step) ? next : gcd(step, next);
			if (step == 0) throwfr("could not find common grid step");
		}
	}
	size = lrint((buffer[atomCount-1] - buffer[0]) / step) + 1;
	origin = buffer[0];
}

//----------------------------------------------------------------------

Broadcaster::Broadcaster(const arma::Mat<double>& atomPositions, int orbitalCount)
		:orbitalCount(orbitalCount), atomOffsets(mpi::size()), atomCounts(mpi::size())
{
	Vector3D<double> origin; // used only by root process
	if (mpi::root()) {
		// root process
		computeGrid(atomPositions, 0, dimensionRaw.x, origin.x, step.x);
		computeGrid(atomPositions, 1, dimensionRaw.y, origin.y, step.y);
		computeGrid(atomPositions, 2, dimensionRaw.z, origin.z, step.z);
		dimensionPad = round4FFT(dimensionRaw);
	}
	mpi::Type typeDimension = createTypeDimension();
	mpi::Type typeVector3D = createTypeVector3D();
	mpi::broadcast(&dimensionRaw, 1, typeDimension);
	mpi::broadcast(&dimensionPad, 1, typeDimension);
	mpi::broadcast(&step, 1, typeVector3D);
	dimensionLocal = DualDimension(dimensionPad).real;

	typeAtomCoeffs = mpi::contiguous(orbitalCount, MPI_DOUBLE_COMPLEX);

	const int mpi_size = mpi::size();
	std::vector<int> zOffsets, atomCellIndices; // used only by root process
	if (mpi::root()) {
		zOffsets.resize(mpi_size);
	}
	mpi::gather(&dimensionLocal.zOffset, 1, MPI_INT, zOffsets.data(), 1);
	if (mpi::root()) {
		const int atomCount = atomPositions.n_cols;
		atomCellIndices.reserve(atomCount);
		for (int node=0, a=0; a<atomCount; ++a) {
			int ix = lrint((atomPositions(0, a) - origin.x) / step.x);
			int iy = lrint((atomPositions(1, a) - origin.y) / step.y);
			int iz = lrint((atomPositions(2, a) - origin.z) / step.z);

			while (node+1 < mpi_size && zOffsets[node+1] > 0 && iz >= zOffsets[node+1]) {
				++node;
				atomOffsets[node] = a;
			}
			++atomCounts[node];
			size_t cellIndex = ((iz - zOffsets[node]) * dimensionLocal.y + iy) * dimensionLocal.x + ix;
			atomCellIndices.push_back(cellIndex);
		}
	}
	mpi::broadcast(atomOffsets.data(), mpi_size, MPI_INT);
	mpi::broadcast(atomCounts.data(), mpi_size, MPI_INT);

	const int localAtomCount = atomCounts[mpi::rank()];
	localCellIndices.resize(localAtomCount);
	MPI_Scatterv(atomCellIndices.data(), atomCounts.data(), atomOffsets.data(), MPI_INT,
			localCellIndices.data(), localAtomCount, MPI_INT, 0, MPI_COMM_WORLD);
}

std::shared_ptr<arma::Mat<complex>> Broadcaster::broadcastData(const arma::Mat<complex>& input)
{
	const int localAtomCount = atomCounts[mpi::rank()];
	std::shared_ptr<arma::Mat<complex>> coeffs = std::make_shared<arma::Mat<complex>>(orbitalCount, localAtomCount);
	MPI_Scatterv(const_cast<complex*>(input.memptr()), atomCounts.data(), atomOffsets.data(), typeAtomCoeffs,
			coeffs->memptr(), localAtomCount, typeAtomCoeffs, 0, MPI_COMM_WORLD);
	return coeffs;
}

std::vector<int> Broadcaster::getLocalCellIndices() const
{
	return localCellIndices;
}

Dimension Broadcaster::getPaddedDimension() const
{
	return dimensionPad;
}

Vector3D<double> Broadcaster::getStepValues() const
{
	return step;
}

mpi::Type Broadcaster::createTypeDimension(void)
{
	int dimensionBlockLengths[3] = {1, 1, 1};
	MPI_Aint dimensionDisplacements[3] = {offsetof(Dimension, x), offsetof(Dimension, y), offsetof(Dimension, z)};
	MPI_Datatype dimensionTypes[3] = {MPI_INT, MPI_INT, MPI_INT};
	return mpi::structure(3, dimensionBlockLengths, dimensionDisplacements, dimensionTypes);
}

mpi::Type Broadcaster::createTypeVector3D(void)
{
	int dimensionBlockLengths[3] = {1, 1, 1};
	MPI_Aint dimensionDisplacements[3] = {offsetof(Vector3D<double>, x), offsetof(Vector3D<double>, y), offsetof(Vector3D<double>, z)};
	MPI_Datatype dimensionTypes[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
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
