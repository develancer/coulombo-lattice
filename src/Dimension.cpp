// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <fftw3-mpi.h>
#include "Dimension.hpp"

static void checkDimensions(int x, int y, int z)
{
	if (x<0 || y<0 || z<0) {
		throw std::logic_error("dimensions are negative");
	}
	if (x && y && z && std::numeric_limits<uword>::max()/x/y/z==0) {
		throw std::logic_error("dimensions are too large");
	}
}

//----------------------------------------------------------------------

Dimension::Dimension(void) { }

Dimension::Dimension(int x, int y, int z)
		:Vector3D<int>{x, y, z}
{
	checkDimensions(x, y, z);
}

uword Dimension::cells(void) const
{
	return static_cast<uword>(x)*static_cast<uword>(y)*static_cast<uword>(z);
}

Dimension Dimension::plusone(void) const
{
	return Dimension(x+1, y+1, z+1);
}

Dimension Dimension::twice(void) const
{
	return Dimension(2*x, 2*y, 2*z);
}

bool Dimension::operator!=(const Dimension& other) const
{
	return x!=other.x || y!=other.y || z!=other.z;
}

//----------------------------------------------------------------------

DistributedDimension::DistributedDimension(void)
		:zOffset(0), zFull(0) { }

DistributedDimension::DistributedDimension(const Dimension& dimension)
		:Dimension(dimension), zOffset(0), zFull(dimension.z) { }

uword DistributedDimension::cellsFull(void) const
{
	return static_cast<uword>(x)*static_cast<uword>(y)*static_cast<uword>(zFull);
}

//----------------------------------------------------------------------

DualDimension::DualDimension(const Dimension& dimension)
{
	ptrdiff_t zRealLength, zRealOffset, yFreqLength, yFreqOffset;
	cells = fftw_mpi_local_size_3d_transposed(
			dimension.z, dimension.y, dimension.x,
			MPI_COMM_WORLD,
			&zRealLength, &zRealOffset,
			&yFreqLength, &yFreqOffset
	);

	real = DistributedDimension({
			dimension.x,
			dimension.y,
			static_cast<int>(zRealLength)
	});
	real.zOffset = zRealOffset;
	real.zFull = dimension.z;

	freq = DistributedDimension({
			dimension.x,
			dimension.z,
			static_cast<int>(yFreqLength)
	});
	// freq domain has y and z swapped
	freq.zOffset = yFreqOffset;
	freq.zFull = dimension.y;
}

//----------------------------------------------------------------------
