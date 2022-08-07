// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_TIMER_HPP
#define COULOMBO_TIMER_HPP

#include <memory>

//----------------------------------------------------------------------

class ElapsingTimer {
	const bool hasMessage_;

	struct timeval start_;

	ElapsingTimer(const ElapsingTimer&);
	ElapsingTimer& operator=(const ElapsingTimer&);

public:
	ElapsingTimer(void);

	ElapsingTimer(const char* msg, va_list ap);

	~ElapsingTimer(void);

	float time(void) const;
};

//----------------------------------------------------------------------

class Timer {
	std::unique_ptr<ElapsingTimer> timer_;

public:
	Timer(void);

	void start(void);

	void start(const char* msg, ...) __attribute__ (( format (printf, 2, 3)));

	float time(void) const;

	void stop(void);
};

//----------------------------------------------------------------------

#endif
