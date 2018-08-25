// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_COULOMBCALCULATOR_HPP
#define COULOMBO_COULOMBCALCULATOR_HPP

#include "base.hpp"
#include "Domain.hpp"
#include "FourierPlan.hpp"
#include "Interaction.hpp"

//----------------------------------------------------------------------

/**
 * Calculator for Coulomb matrix elements based on quasi-densities on 3-D grid.
 *
 * It should be used as follows:
 * 1. Instantiate the object with the dimension (X×Y×Z) of the wavefunctions.
 * 2. Set values of the interaction function into domain G.
 * 3. Call initialize().
 *   4. Set values of the first quasi-density into domain F.
 *   5. Call prepare().
 *     6. Set values of the second quasi-density into domain F.
 *     7. Call calculate() and collect the value returned in MPI root process.
 *   (points 6-7 can be repeated to only change the second quasi-density)
 * (points 4-7 can be repeated to change both quasi-densities)
 */
class CoulombCalculator {
	DualDomain<complex> Ftemp;
	SingleDomain<real> Gfreq;
	SingleDomain<complex> F_;
	DualDomain<real> G_;
	SingleDomain<complex> V_;

	FourierPlan planForwardF;
	FourierPlan planInverseF;
	FourierPlan planForwardG;

	static DistributedDimension computeGfreqDimension(const DistributedDimension& freqDimension);

	std::vector<complex> phaseFactorX, phaseFactorY, phaseFactorZ;

public:
	/**
	 * 3-D domain for the input of quasi-densities.
	 * It will have the same dimensions (X×Y×Z) as passed to the constructor.
	 */
	Domain<complex>& input;

	/**
	 * Create a new calculator for the specific dimension.
	 *
	 * @param dimension dimensions of the X×Y×Z cube (i.e. F)
	 * @param measureFFT if true, constructor will take longer, but the generated FFT plans will be more efficient
	 * (should be set to true iff the number of FFT transforms to be computed is large)
	 */
	CoulombCalculator(const Dimension& dimension, bool measureFFT = true);

	/**
	 * Initialize the calculator with a given interaction function and step lengths.
	 *
	 * @param interaction interaction (G) function to be used
	 * @param step consists of step lengths in all directions
	 */
	void initialize(const Interaction& interaction);

	/**
	 * Prepare the calculator with one of the quasi-densities
	 * which will be later used for calculation of Coulomb integrals.
	 * Values of the first quasi-density (A) should be
	 * set into domain F prior to calling this method.
	 * Also, method initialize() should have been already called.
	 */
	void prepare(void);

	/**
	 * Calculate the value of the Coulomb integral
	 *   I = ∫∫ A(rₐ) G(|rₐ-rₑ|) E(rₑ) drₐ³ drₑ³
	 * using convolution and fast Fourier transform.
	 * Values of the second quasi-density (E) should be set into domain F
	 * prior to calling this method. Also, methods initialize() and prepare()
	 * should have been already called.
	 */
	complex calculate(void);
};

//----------------------------------------------------------------------

#endif
