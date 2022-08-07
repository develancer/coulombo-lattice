// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>

#include "mpi.hpp"
#include "Timer.hpp"

ElapsingTimer::ElapsingTimer(void)
		:hasMessage_(false)
{
	gettimeofday(&start_, 0);
}

ElapsingTimer::ElapsingTimer(const char* msg, va_list ap)
		:hasMessage_(msg)
{
	if (hasMessage_) {
		vfprintf(stderr, msg, ap);
		fputs("... ", stderr);
		fflush(stderr);
	}

	gettimeofday(&start_, 0);
}

ElapsingTimer::~ElapsingTimer(void)
{
	if (hasMessage_) {
		fprintf(stderr, "(%.3f s)\n", time());
		fflush(stderr);
	}
}

float ElapsingTimer::time(void) const
{
	struct timeval now;
	gettimeofday(&now, 0);
	return (now.tv_usec-start_.tv_usec)*1.0e-6+(now.tv_sec-start_.tv_sec);
}

//----------------------------------------------------------------------

Timer::Timer(void) { }

void Timer::start(void)
{
	mpi::barrier();
	if (mpi::rank()) return;
	timer_.reset();
	timer_.reset(new ElapsingTimer());
}

void Timer::start(const char* msg, ...)
{
	mpi::barrier();
	if (mpi::rank()) return;
	va_list ap;
	va_start(ap, msg);
	timer_.reset();
	timer_.reset(new ElapsingTimer(msg, ap));
	va_end(ap);
}

float Timer::time(void) const
{
	return timer_.get() ? timer_->time() : 0.0f;
}

void Timer::stop(void)
{
	mpi::barrier();
	if (mpi::rank()) return;
	timer_.reset();
}
