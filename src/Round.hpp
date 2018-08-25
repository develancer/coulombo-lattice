// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_ROUND_HPP
#define COULOMBO_ROUND_HPP

#include <algorithm>
#include <cassert>

//----------------------------------------------------------------------

class Round {
public:
	/**
	 * @tparam first > 1
	 * @param n
	 * @param c ≥ 1
	 * @return smallest integer ≥ n of form (c · first^i) for some i ∊ ℕ
	 */
	template<int first>
	static int up(int n, int c = 1)
	{
		static_assert(first > 1, "must be > 1");
		assert(c >= 1);
		while (c < n) c *= first;
		return c;
	}

	/**
	 * @tparam first > 1
	 * @tparam second > 1
	 * @tparam other all > 1
	 * @param n
	 * @param c ≥ 1
	 * @return smallest integer ≥ n of form (c · first^i₁ · second^i₂ · ...)
	 *         for some i₁, i₂, ... ∊ ℕ
	 */
	template<int first, int second, int... other>
	static int up(int n, int c = 1)
	{
		static_assert(first > 1, "must be > 1");
		assert(c >= 1);
		int best = up<second, other...>(n, c);
		while (c < n) {
			c *= first;
			best = std::min(best, up<second, other...>(n, c));
		}
		return best;
	}
};

//----------------------------------------------------------------------

#endif
