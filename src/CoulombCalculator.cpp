// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#include <cassert>
#include <fftw3-mpi.h>
#include "mpi.hpp"
#include "CoulombCalculator.hpp"

//----------------------------------------------------------------------

CoulombCalculator::CoulombCalculator(const Dimension& dimension, bool measureFFT)
		: Ftemp(DualDimension(dimension)),
		Gfreq(computeGfreqDimension(Ftemp.freq.dimension)),
		F_(Ftemp.dimension),
		G_(DualDimension(dimension.plusone())),
		V_(Ftemp.dimension),
		planForwardF(fftw_mpi_plan_dft_3d(
				dimension.z, dimension.y, dimension.x,
				reinterpret_cast<fftw_complex*>(Ftemp.memptr()),
				reinterpret_cast<fftw_complex*>(Ftemp.memptr()),
				MPI_COMM_WORLD,
				FFTW_FORWARD,
				FFTW_MPI_TRANSPOSED_OUT | (measureFFT ? FFTW_MEASURE : FFTW_ESTIMATE)
		)),
		planInverseF(fftw_mpi_plan_dft_3d(
				dimension.z, dimension.y, dimension.x,
				reinterpret_cast<fftw_complex*>(Ftemp.memptr()),
				reinterpret_cast<fftw_complex*>(Ftemp.memptr()),
				MPI_COMM_WORLD,
				FFTW_BACKWARD,
				FFTW_MPI_TRANSPOSED_IN | (measureFFT ? FFTW_MEASURE : FFTW_ESTIMATE)
		)),
		planForwardG(fftw_mpi_plan_r2r_3d(
				G_.dimension.zFull, G_.dimension.y, G_.dimension.x,
				G_.memptr(), G_.memptr(),
				MPI_COMM_WORLD,
				FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,
				FFTW_MPI_TRANSPOSED_OUT | FFTW_ESTIMATE
		)),
		phaseFactorX(Ftemp.dimension.x),
		phaseFactorY(Ftemp.dimension.y),
		phaseFactorZ(Ftemp.dimension.z),
		input(F_)
{
	const double wx = M_PI/Ftemp.dimension.x, wy = M_PI/Ftemp.dimension.y, wz = M_PI/Ftemp.dimension.zFull;
	for (int ix = 0; ix<Ftemp.dimension.x; ++ix) {
		phaseFactorX[ix] = std::polar(1.0, wx*ix);
	}
	for (int iy = 0; iy<Ftemp.dimension.y; ++iy) {
		phaseFactorY[iy] = std::polar(1.0, wy*iy);
	}
	for (int iz = 0; iz<Ftemp.dimension.z; ++iz) {
		phaseFactorZ[iz] = std::polar(1.0, wz*(iz+Ftemp.dimension.zOffset));
	}
}

void CoulombCalculator::initialize(const Interaction& interaction)
{
	interaction.map(G_);

	planForwardG.execute();

	const int size = mpi::size();
	// in current node we have part of data determined by G_.dims.freq
	int weHaveMin = G_.freq.dimension.zOffset;
	int weHaveMax = weHaveMin+G_.freq.dimension.z-1;
	// but the part needed for convolution will be different
	int weNeedMin = Gfreq.dimension.zOffset;
	int weNeedMax = weNeedMin+Gfreq.dimension.z-1;
	// let's exchange these dimensions with other nodes
	int haveMin[size], haveMax[size], needMin[size], needMax[size];
	mpi::allgather(&weHaveMin, 1, MPI_INT, haveMin, 1);
	mpi::allgather(&weHaveMax, 1, MPI_INT, haveMax, 1);
	mpi::allgather(&weNeedMin, 1, MPI_INT, needMin, 1);
	mpi::allgather(&weNeedMax, 1, MPI_INT, needMax, 1);
	// we have to calculate which parts need to be sent and received
	int sendcounts[size], sendoffsets[size], recvcounts[size], recvoffsets[size];
	const int sliceSize = G_.freq.dimension.x*G_.freq.dimension.y;

	mpi::Type typeSlice = mpi::contiguous(sliceSize, MPI_DOUBLE);
	for (int node = 0; node<size; ++node) {
		// counting in number of slices (sliceSize each)
		sendcounts[node] = recvcounts[node] = 0;
		sendoffsets[node] = recvoffsets[node] = 0;
		if (weHaveMin<=needMax[node] && needMin[node]<=weHaveMax) {
			int first = std::max(weHaveMin, needMin[node]);
			int last = std::min(weHaveMax, needMax[node]);
			sendoffsets[node] = first-weHaveMin;
			sendcounts[node] = last-first+1;
		}
		if (weNeedMin<=haveMax[node] && haveMin[node]<=weNeedMax) {
			int first = std::max(weNeedMin, haveMin[node]);
			int last = std::min(weNeedMax, haveMax[node]);
			recvoffsets[node] = first-weNeedMin;
			recvcounts[node] = last-first+1;
		}
	}
	// now we gather needed data
	mpi::alltoallv(G_.memptr(), sendcounts, sendoffsets, typeSlice, Gfreq.memptr(), recvcounts, recvoffsets);

	// FFT-related normalization
	Gfreq /= (F_.dimension.cellsFull()*8u);
}

void CoulombCalculator::prepare(void)
{
	V_.zeros();
	for (unsigned round = 0; round<8; ++round) {
		bool kx = round & 1;
		bool ky = round & 2;
		bool kz = round & 4;

		Ftemp = F_;

		if (Ftemp.dimension.z) {
			#pragma omp parallel for schedule(static)
			for (int iz = 0; iz<Ftemp.dimension.z; ++iz) {
				const complex pfZ = kz ? std::conj(phaseFactorZ[iz]) : 1;
				for (int iy = 0; iy<Ftemp.dimension.y; ++iy) {
					const complex pfYZ = ky ? std::conj(phaseFactorY[iy])*pfZ : pfZ;
					for (int ix = 0; ix<Ftemp.dimension.x; ++ix) {
						Ftemp(ix, iy, iz) *= kx ? std::conj(phaseFactorX[ix])*pfYZ : pfYZ;
					}
				}
			}
		}

		planForwardF.execute();

		#pragma omp parallel for schedule(static)
		for (int iz = 0; iz<Ftemp.freq.dimension.z; ++iz) {
			int izG = 2*(iz+Ftemp.freq.dimension.zOffset)+ky;
			for (int iy = 0; iy<Ftemp.freq.dimension.y; ++iy) {
				int iyG = 2*iy+kz; // ky and kz are swapped due to data transposition
				for (int ix = 0; ix<Ftemp.freq.dimension.x; ++ix) {
					int ixG = 2*ix+kx;
					Ftemp.freq(ix, iy, iz) *= Gfreq(
							std::min(ixG, 2*Ftemp.freq.dimension.x-ixG),
							std::min(iyG, 2*Ftemp.freq.dimension.y-iyG),
							std::min(izG, 2*Ftemp.freq.dimension.zFull-izG)-Gfreq.dimension.zOffset
					);
				}
			}
		}

		planInverseF.execute();

		if (Ftemp.dimension.z) {
			#pragma omp parallel for schedule(static)
			for (int iz = 0; iz<Ftemp.dimension.z; ++iz) {
				const complex pfZ = kz ? phaseFactorZ[iz] : 1;
				for (int iy = 0; iy<Ftemp.dimension.y; ++iy) {
					const complex pfYZ = ky ? phaseFactorY[iy]*pfZ : pfZ;
					for (int ix = 0; ix<Ftemp.dimension.x; ++ix) {
						Ftemp(ix, iy, iz) *= kx ? phaseFactorX[ix]*pfYZ : pfYZ;
					}
				}
			}
		}

		V_ += Ftemp;
	}
}

complex CoulombCalculator::calculate(void)
{

	complex part = arma::accu(V_%F_), result;

	mpi::reduce(&part, &result, 1, MPI_DOUBLE_COMPLEX, MPI_SUM);
	return result;
}

DistributedDimension CoulombCalculator::computeGfreqDimension(const DistributedDimension& freqDimension)
{
	int weNeedMin = 2*freqDimension.zOffset;
	int weNeedMax = weNeedMin+2*freqDimension.z-1;
	int maxIndex = freqDimension.zFull;
	if (weNeedMax>maxIndex) {
		if (weNeedMin>maxIndex) {
			std::swap(weNeedMin, weNeedMax);
			weNeedMin = 2*maxIndex-weNeedMin;
			weNeedMax = 2*maxIndex-weNeedMax;
		}
		else {
			weNeedMin = std::min(weNeedMin, 2*maxIndex-weNeedMax);
			weNeedMax = maxIndex;
		}
	}
	DistributedDimension result(freqDimension.plusone());
	result.z = weNeedMax-weNeedMin+1;
	result.zOffset = weNeedMin;
	result.zFull = freqDimension.z+1;
	return result;
}

//----------------------------------------------------------------------
