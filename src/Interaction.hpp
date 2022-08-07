// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_INTERACTION_HPP
#define COULOMBO_INTERACTION_HPP

#include "base.hpp"
#include "Domain.hpp"

//----------------------------------------------------------------------

class Interaction {
public:
	virtual void map(Domain<real>& G) const =0;
};

//----------------------------------------------------------------------

template<class IMPL>
class InteractionBase : public Interaction {
public:
	InteractionBase(const Vector3D<double>& step, double onsite)
			:step(step), onsite(onsite) { }

	/**
	 * We assume the interaction to be symmetric, i.e.
	 *   G(±x,±y,±z) = G(x,y,z)
	 * so only one-eighth with positive coordinates has to be defined.
	 *
	 * @param G domain to map interaction values into
	 */
	void map(Domain<real>& G) const
	{
		const double constant = E2_4PE0;
		const int zStart = G.dimension.zOffset, zEnd = zStart+G.dimension.z;
		const IMPL* self = static_cast<const IMPL*>(this);
		for (int iz = zStart; iz<zEnd; ++iz) {
			double z = iz*step.z, z2 = z*z;
			for (int iy = 0; iy<G.dimension.y; ++iy) {
				double y = iy*step.y, y2z2 = z2+y*y;
				for (int ix = 0; ix<G.dimension.x; ++ix) {
					double x = ix*step.x, r2 = y2z2+x*x, r = std::sqrt(r2);
					G(ix, iy, iz-zStart) = constant / (self->dielectric(r) * r);
				}
			}
		}
		if (!zStart && zEnd) {
			// only one process has point (0,0,0)
			G(0, 0, 0) = onsite;
		}
	}

private:
	const Vector3D<double> step;
	const double onsite;
};

//----------------------------------------------------------------------

class InteractionSimple: public InteractionBase<InteractionSimple> {
public:
	InteractionSimple(const Vector3D<double>& step, double onsite, double dielectric);

	double dielectric(double r) const;

private:
	const double dielectric_;
};

//----------------------------------------------------------------------

class InteractionThomasFermi : public InteractionBase<InteractionThomasFermi> {
public:
	InteractionThomasFermi(const Vector3D<double>& step, double onsite, double dielectric, double latticeConstant);

	double dielectric(double r) const;

private:
	const double dielectric_;
	const double latticeConstant;
	double rTF_, qTF_;
};

//----------------------------------------------------------------------

#endif
