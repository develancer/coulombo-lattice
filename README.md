Coulombo
====

[Computer Physics Communications] Różański & Zieliński:  
Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures

## What is it?

Coulombo is a ready-to-use implementation for calculating Coulomb matrix
elements for a given set of input wavefunctions. This implementation is based on
the approach introduced in [1], using fast Fourier transform to compute
convolution. In this work we further significantly improved the method by
eliminating the need to extend the computational domain with padding, thus
reducing the memory consumption by a factor of 8. The optimal computational plan
for each run is prepared by calculating a minimal vertex cover on a graph
representing a subset of requested Coulomb matrix elements.

The implementation is fully parallelized in distributed-memory model, using MPI
and parallel routines from FFTW [2]. Minimal vertex cover is computed by a
greedy approximation algorithm, which we found to perform significantly better
than the standard text-book heuristic.

[1] P. T. Różański & M. Zieliński,  
“Linear scaling approach for atomistic calculation of excitonic properties of
10-million-atom nanostructures”,
Phys. Rev. B 94 (2016) 045440

[2] M. Frigo & Steven G. Johnson,  
“The Design and Implementation of FFTW3”,
Proceedings of the IEEE 93 (2), 216–231 (2005), Invited paper, Special Issue on
Program Generation, Optimization, and Platform Adaptation

## How to use it?

### Compilation

Run “make”. This will produce a single executable file called “coulombo”.
To make it available system-wide, run “sudo make install”, which will copy it to
/usr/local/bin/.

#### Tests

To run the unit tests, execute “make test”. This will compile the test
executable and run it, resulting (hopefully) in displaying the “OK” message.

#### Requirements

To compile the software, one needs an C++ compiler with support for 2011
standard and a MPI wrapper (mpic++). OpenMP support is also enabled by default.

Coulombo has two external library dependencies:

* Fast Fourier Transform implementation _FFTW_ (version 3),
* _armadillo_ C++ linear algebra library.

Both of the above can be installed in a standard way (e.g. apt) in the most of
the major Linux distributions.

### Running Coulombo

Coulombo should be run with a set of command line parameters. Results of the
computation (values of the Coulomb matrix elements) will be written to the
standard output, whereas all diagnostic and error message will be written to the
standard error stream.

Most important part of command line parameters are the paths to the wavefunction
data in form of binary files in _armadillo_ machine-dependent format (one file
per wavefunction, unless --spin flag is given, see below). At least one input
file path should be specified.

List of all possible parameters can be displayed by running coulombo with no
additional arguments. They are, in alphabetical order:

* --**dielectric**=VALUE specifies the relative dielectric constant for the
calculation of Coulomb potential. If this option is not given, the dielectric
constant of 1 (as for vacuum) is assumed.

> Example: `coulombo ... --dielectric=10.89`  
> for high-frequency dielectric constant of GaAs

* --**integrals**=LIST specifies the comma-separated list of integrals to be
computed. If this option is not given, all possible integrals (M⁴ where M is
the number of wavefunctions) will be computed.

> Example: `coulombo ... --integrals=1212 f1.bin f2.bin`  
> computes only the single integral 1212  
> (defined as ∫ f₁\*(r⃗) f₂\*(r⃗') G(r⃗-r⃗') f₁(r⃗') f₂(r⃗) dr³ dr'³)

> Example: `coulombo ... --integrals=ijji f1.bin f2.bin`  
> computes integrals 1111, 1221, 2112, 2222  

> Example: `coulombo ... --integrals=1*1* f1.bin f2.bin`  
> computes integrals 1111, 1112, 1211, 1212  

* --**spin** specifies that the calculations should include spin. Each
wavefunction should be specified with two consecutive file paths, representing
spin-up and spin-down parts of the wavefunction. Order of up/down parts is not
relevant, as long as it is consistent across all wavefunctions.

> Example: `coulombo ... --spin f1-U.bin f1-D.bin f2-U.bin f2-D.bin`  
> assumes that the calculation should be performed for two wavefunctions: first
> given by files f1-U.bin and f1-D.bin and the second given by files f2-U.bin
and f2-D.bin.

* --**step**=VALUE specifies the size of the grid step (in ångströms) of the
provided wavefunction files. This option is mandatory.

> Example: `coulombo --step=0.8 ...`  
> assumes that the wavefunction files have the grid step of 0.8 Å (0.08 nm)

* --**threads-per-node**=N specifies the requested number of OpenMP threads per each
MPI node. If this option is not given, one thread per node is assumed and no
OpenMP parallelization is used.

> Example: `coulombo --threads-per-node=4 ...`  
> runs four OpenMP threads per each MPI node

* --**tf-lattice**=VALUE specifies the lattice constant for the Thomas-Fermi-Resta
dielectric screening model. If the option is not given, only the simple
screening model (defined by _dielectric_ parameter) is used.

> Example: `coulombo --tf-lattice=5.65325 ...`  
> uses the Thomas-Fermi-Resta model for the lattice constant of GaAs

#### Example

Subdirectory “example” contains a small example consisting of two wavefunctions.
The wavefunctions represent the first electron state (LUMO) and the second hole
state (a Kramers-degenerate counterpart to the HOMO state) of a InAs/InP quantum
dot. The example is run as

```
coulombo --dielectric=12.4 --spin --step=10 \
         --integrals=1111,1122,1221,2222 \
         e1-D.arma e1-U.arma h2-D.arma h2-U.arma
```

by requesting four integrals:

* `1111` equivalent to J<sub>ee</sub> electron-electron Coulomb integral,
* `1122` equivalent to exchange integral ⟨e₁h₁e₂h₂⟩ due to Kramers degeneracy,
* `1221` equivalent to J<sub>eh</sub> electron-hole Coulomb integral,
* `2222` equivalent to J<sub>hh</sub> hole-hole Coulomb integral.

The results are as follows:

```
 1  1  1  1   0.017643339  0.000000000
 1  1  2  2  -0.000026347  0.000040578
 1  2  2  1   0.017413492  0.000000000
 2  2  2  2   0.017726404  0.000000000
```

which, by taking absolute values of the integrals, translate to

* J<sub>ee</sub> = 17.64 meV
* ⟨e₁h₁e₂h₂⟩ = 48.4 µeV
* J<sub>eh</sub> = 17.41 meV
* J<sub>hh</sub> = 17.72 meV

### Licence

This software is provided under a Creative Commons Attribution 4.0 International
Public License, which should be provided as the LICENCE file.
