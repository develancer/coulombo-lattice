// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_FOURIERPLAN_HPP
#define COULOMBO_FOURIERPLAN_HPP

#include <cassert>
#include <memory>
#include <fftw3.h>
#include "base.hpp"

//----------------------------------------------------------------------

/**
 * RAII-style wrapper for FFTW plan structure.
 */
class FourierPlan {
	std::shared_ptr<fftw_plan_s> plan;

public:
	/**
	 * Create a wrapper for an existing FFT plan.
	 */
	inline FourierPlan(fftw_plan plan)
			:plan(plan, fftw_destroy_plan)
	{
		if (!plan) {
			throw std::runtime_error("failed to create plan for FFT");
		}
	}

	/**
	 * Execute plan.
	 */
	inline void execute(void) const
	{
		fftw_execute(plan.get());
	}
};

//----------------------------------------------------------------------

#endif
