// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_DIMENSION_HPP
#define COULOMBO_DIMENSION_HPP

#include "base.hpp"

//----------------------------------------------------------------------

/**
 * Simple structure for storing dimensions of the 3-D domains,
 * with values of "x", "y" and "z".
 */
struct Dimension: public Vector3D<int> {

	/**
	 * Set all dimensions to zero.
	 */
	Dimension(void);

	/**
	 * Set all dimensions to given values.
	 * @throws logic_error if dimensions are invalid
	 */
	Dimension(int x, int y, int z);

	/**
	 * @return number of cells in domain of these dimensions
	 */
	uword cells(void) const;

	/**
	 * @return increased by 1 cell in every dimension
	 */
	Dimension plusone(void) const;

	/**
	 * @return augmented by a factor of 2 in every dimension
	 */
	Dimension twice(void) const;

	bool operator!=(const Dimension& other) const;
};

//----------------------------------------------------------------------

/**
 * Simple structure for storing dimensions of the 3-D domains
 * distributed among several MPI nodes. The domains are always
 * divided in the "z" dimension. This structure has two additional
 * fields for the information about memory fragmentation.
 *
 * Example:
 * Domain 400×400×400 is divided equally to four MPI nodes.
 * All of the nodes will have x=400, y=400, z=100, and zFull=400.
 * First node will have zOffset=0, second: zOffset=100, and so forth.
 */
struct DistributedDimension: public Dimension {
	int zOffset;
	int zFull;

	/**
	 * Set all dimensions to zero.
	 */
	DistributedDimension(void);

	/**
	 * Copy all dimensions from given source object.
	 * Set zOffset=0, zFull=z (as if the domain was not fragmented at all)
	 */
	DistributedDimension(const Dimension& dimension);

	/**
	 * @return number of cells in the entire (combined) domain
	 */
	uword cellsFull(void) const;
};

//----------------------------------------------------------------------

/**
 * Structure for storing dimension information for DualDomain.
 */
struct DualDimension {
	/**
	 * @return number of cells in the entire (combined) domain
	 */
	uword cells;

	/**
	 * Dimensions for the real-space view of the domain.
	 * View will be fragmented in z direction.
	 */
	DistributedDimension real;

	/**
	 * Dimensions for the Fourier-space view of this domain.
	 * Data will be transposed, i.e. although view will be fragmented
	 * in z direction (as always), this will correspond to
	 * the y direction in the original coordinate system.
	 */
	DistributedDimension freq;

	/**
	 * Compute a dimensions for DualDomain based on a given dimensions
	 * for the entire (combined) real-space domain.
	 *
	 * Example:
	 * If dimension=[100,200,400] and the dimension is to be divided
	 * equally among 4 processes, then dimension for real-space view
	 * will be [100,200,100] (fragmented in z direction) and the
	 * dimension for the Fourier-space view will be [100,400,50]
	 * as the y and z dimensions are swapped and fragmented z direction
	 * corresponds to the original direction of y.
	 *
	 * @param dimension  dimensions for the real-space domain
	 */
	explicit DualDimension(const Dimension& dimension);
};

//----------------------------------------------------------------------

#endif
