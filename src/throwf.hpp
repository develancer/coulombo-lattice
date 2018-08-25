// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_THROWF_HPP
#define COULOMBO_THROWF_HPP

#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#define throwfr throwf<std::runtime_error>
#define throwfl throwf<std::logic_error>

template<class E>
void throwf(const char* format, ...) __attribute__ (( format (printf, 1, 2), noreturn ));

template<class E>
void throwf(const char* format, ...)
{
	int result;
	va_list ap;
	char bufor[1024];

	va_start(ap, format);
	result = vsnprintf(bufor, sizeof bufor, format, ap);
	va_end(ap);

	if (result>0) {
		bufor[sizeof bufor-1] = 0;
		throw E(bufor);
	}
	else {
		throw E("error");
	}
}

#endif
