// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "Interaction.hpp"

//----------------------------------------------------------------------

InteractionSimple::InteractionSimple(const Vector3D<double>& step, double onsite, double dielectric)
		:InteractionBase<InteractionSimple>(step, onsite), dielectric_(dielectric) { }

double InteractionSimple::dielectric(double) const
{
	return dielectric_;
}

//----------------------------------------------------------------------

InteractionThomasFermi::InteractionThomasFermi(const Vector3D<double>& step, double onsite, double dielectric, double latticeConstant)
		:InteractionBase<InteractionThomasFermi>(step, onsite), dielectric_(dielectric), latticeConstant(latticeConstant) {
	const double eps = 1.0e-12;
	qTF_ = 2.0/sqrt(M_PI) * cbrt(96.0*M_PI*M_PI) / latticeConstant;
	double x = sqrt(6.0*(dielectric_ - 1.0)), dx;
	do
	{
		double F = sinh(x) - dielectric_*x;
		double dF = cosh(x) - dielectric_;
		dx = (fabs(F) < eps) ? 0 : F / dF;
		x -= dx;
	}
	while (dx/qTF_ > eps);
	rTF_ = x / qTF_;
}

double InteractionThomasFermi::dielectric(double r) const
{
	return (r<rTF_)
		   ? dielectric_ * ( qTF_ * rTF_ / ( sinh(qTF_*(rTF_-r)) + qTF_*r ) )
		   : dielectric_;
}

//----------------------------------------------------------------------
