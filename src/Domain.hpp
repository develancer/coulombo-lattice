// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_DOMAIN_HPP
#define COULOMBO_DOMAIN_HPP

#include <vector>
#include "Dimension.hpp"
#include "DomainAllocator.hpp"

//----------------------------------------------------------------------

template<typename VALUE>
using DomainData = std::vector<VALUE, DomainAllocator<VALUE>>;

//----------------------------------------------------------------------

/**
 * Base class for pre-allocated 3-D domains (grids) of given type.
 * All methods of arma::Cube (see armadillo library) may be used.
 * Domains represent a fragmented (due to MPI parallelization)
 * 3-D grid of numbers (usually real or complex).
 *
 * @tparam VALUE  type of elements to be stored, e.g. double for real-valued domain
 */
template<typename VALUE>
class Domain: public arma::Cube<VALUE> {
public:
	const DistributedDimension dimension;

	Domain(VALUE* data, const DistributedDimension& dimension)
			:arma::Cube<VALUE>(data, dimension.x, dimension.y, dimension.z, false, true),
			 dimension(dimension) { }

	void operator=(const Domain& source)
	{
		arma::Cube<VALUE>::operator=(source);
	}

	template<class SOURCE>
	void operator=(const SOURCE& source)
	{
		arma::Cube<VALUE>::operator=(source);
	}
};

//----------------------------------------------------------------------

/**
 * Simple pre-allocated 3-D grid of given type and dimension.
 *
 * @tparam VALUE  type of elements to be stored, e.g. double for real-valued domain
 */
template<typename VALUE>
class SingleDomain: protected DomainData<VALUE>, public Domain<VALUE> {
protected:
	SingleDomain(const DistributedDimension& dimension, uword cells)
			:DomainData<VALUE>(cells),
			 Domain<VALUE>(DomainData<VALUE>::data(), dimension) { }

public:
	SingleDomain(const DistributedDimension& dimension)
			:DomainData<VALUE>(dimension.cells()),
			 Domain<VALUE>(DomainData<VALUE>::data(), dimension) { }

	VALUE& operator[](uword index)
	{
		return DomainData<VALUE>::operator[](index);
	}

	VALUE operator[](uword index) const
	{
		return DomainData<VALUE>::operator[](index);
	}

	using Domain<VALUE>::operator=;
};

//----------------------------------------------------------------------

/**
 * 3-D grid of given type and dimension, suitable for executing
 * Fourier transforms in-place. In addition to normal (real-space) access,
 * this object has also a Fourier-space view (freq) with potentially
 * different dimensions.
 *
 * @tparam VALUE  type of elements to be stored, e.g. double for real-valued domain
 */
template<typename VALUE>
class DualDomain: public SingleDomain<VALUE> {
public:
	/**
	 * Fourier-space view for this domain.
	 */
	Domain<VALUE> freq;

	DualDomain(const DualDimension& dims)
			:SingleDomain<VALUE>(dims.real, dims.cells),
			 freq(DomainData<VALUE>::data(), dims.freq) { }

	using Domain<VALUE>::operator=;
};

//----------------------------------------------------------------------

#endif
