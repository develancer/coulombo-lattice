// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_DOMAINALLOCATOR_HPP
#define COULOMBO_DOMAINALLOCATOR_HPP

#include <cstdlib>
#include <cstring>
#include <limits>
#include <fftw3.h>
#include "base.hpp"

//----------------------------------------------------------------------

/**
 * Provides an FFTW-friendly allocator for standard containers (e.g. std::vector).
 * (Inspired by Stephan Lavavej's Mallocator.)
 * @tparam T  type of elements to be allocated, e.g. double for real-valued arrays
 */
template<class T>
class DomainAllocator {
public:
	typedef T value_type;

	DomainAllocator() noexcept { }

	template<class U>
	DomainAllocator(const DomainAllocator<U>&) noexcept { }

	template<class U>
	bool operator==(const DomainAllocator<U>&) const noexcept
	{
		return true;
	}

	template<class U>
	bool operator!=(const DomainAllocator<U>&) const noexcept
	{
		return false;
	}

	T* allocate(const size_t n) const
	{
		if (!n) return nullptr;
		if (n>std::numeric_limits<size_t>::max()/sizeof(T)) {
			throw std::bad_array_new_length();
		}
		size_t size = n*sizeof(T);
		void* const result = fftw_malloc(size);
		if (!result) throw std::bad_alloc();
		memset(result, 0, size);
		return static_cast<T*>(result);
	}

	void deallocate(T* const p, size_t) const noexcept
	{
		fftw_free(p);
	}
};

//----------------------------------------------------------------------

#endif
