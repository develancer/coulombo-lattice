Coulombo-Lattice
====

P. T. Różański & M. Zieliński:
Exploiting underlying crystal lattice for efficient computation of Coulomb matrix elements in multi-million atoms nanostructures,
Computer Physics Communications 287 (2023), pp. 108693,
DOI: [10.1016/j.cpc.2023.108693](https://doi.org/10.1016/j.cpc.2023.108693)

## What is it?

Coulombo-Lattice is a ready-to-use implementation for calculating Coulomb matrix
elements for a given set of input wave-functions given in the form of LCAO
(linear combination of atomic orbitals), typically resulting from the
empirical tight binding approach. This implementation is build on top of the method
introduced in [3], using fast Fourier transform without zero-padding the
computational domain, allowing to achieve the quasi-linear scaling with
minimal memory footprint.

In this work we adapt the method to use the LCAO functions directly, without a need
to introduce any auxiliary set of atomic orbitals (e.g. Slater orbitals). Instead,
the method identifies and explores the underlying grid associated with the regular crystal
lattice in order to perform efficient numerical calculations on a regular, three-dimensional
grid.

Similarly to the wavefunction-based version [3], this implementation is fully parallelized
in distributed-memory model, using MPI and parallel routines from FFTW [2]. 

[1] P. T. Różański & M. Zieliński,  
“Linear scaling approach for atomistic calculation of excitonic properties of
10-million-atom nanostructures”,
Phys. Rev. B 94 (2016) 045440

[2] M. Frigo & Steven G. Johnson,  
“The Design and Implementation of FFTW3”,
Proceedings of the IEEE 93 (2), 216–231 (2005), Invited paper, Special Issue on
Program Generation, Optimization, and Platform Adaptation

[3] P. T. Różański & M. Zieliński,
“Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures”,
Computer Physics Communications 238 (2019), pp. 254–261

### Package directory structure

* example — includes an example input/output specification with a shell script
* src — includes all auxiliary C++ source files (the header files as well)
* coulombo.cpp — main C++ source file (the one with `main()` function)
* LICENCE — licence for the program distribution
* Makefile — makefile script for program compilation and testing
* README.md — description of the package (this file)
* run-tests.cpp — additional source file consisting of all unit tests

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
computation (values of the Coulomb matrix elements) will be written to a set of text
files, named eeee.txt, hhhh.txt, ehhe.txt and so on,
whereas all diagnostic and error message will be written to the
standard error stream.

Most important part of command line parameters are the paths to input files:
* exactly one path to a text file representing atom positions, given as `--atoms` parameter,
where each line corresponds to the coordinates of a single atom in angstroms (Å) as
```
1st_atom_x 1st_atom_y 1st_atom_z
2nd_atom_x 2nd_atom_y 2nd_atom_z
...
```
* one or more paths to text files, each starting with `e` or `h`
representing wavefunction data in LCAO form. Each file consists of the vector of all
LCAO coefficients, one per line, in the following order:
```
1st_atom_1st_orbital_real_part 1st_atom_1st_orbital_imaginary_part
1st_atom_2nd_orbital_real_part 1st_atom_2nd_orbital_imaginary_part
...
2nd_atom_1st_orbital_real_part 2nd_atom_1st_orbital_imaginary_part
2nd_atom_2nd_orbital_real_part 2nd_atom_2nd_orbital_imaginary_part
...
```

Please see the `example` subdirectory for an actual example.

List of all possible parameters can be displayed by running coulombo with no
additional arguments. They are, in alphabetical order:

* --**dielectric**=VALUE specifies the relative dielectric constant for the
calculation of Coulomb potential. If this option is not given, the dielectric
constant of 1 (as for vacuum) is assumed.

> Example: `coulombo ... --dielectric=10.89`  
> for high-frequency dielectric constant of GaAs

* --**integrals**=LIST specifies the comma-separated list of integrals to be
computed. If this option is not given, the integrals will be restricted to `eeee,hhhh,ehhe,eheh`.
Each integral specification may consist of
  * a digit `1`-`9`, representing a single, specified state (`1` represent the first state etc.)
  * a letter `a`-`z` except `e` or `h`, representing any state, but all occurrences of the same letter
  correspond to each other
  (e.g. `ijki` represent all integrals in which first and fourth state is the same)
  * a letter `e` or `h`, representing any state depending on whether its name starts with e or h,
  but in contrary to the previous bullet, here all occurrences of `e` and `h` are independent
  * an asterisk `*`, representing any state, and all occurrences of `*` are independent

> Example: `coulombo ... --integrals=1212 h1.bin e1.bin`  
> computes only the single integral 1212  
> (defined as ∫ f₁*(r⃗) f₂*(r⃗') G(r⃗-r⃗') f₁(r⃗') f₂(r⃗) dr³ dr'³)

> Example: `coulombo ... --integrals=ijji h1.bin e1.bin`  
> computes integrals 1111, 1221, 2112, 2222  

> Example: `coulombo ... --integrals=1*1* h1.bin e1.bin`  
> computes integrals 1111, 1112, 1211, 1212  

> Example: `coulombo ... --integrals=ehhe h1.bin h2.bin e1.bin`  
> computes integrals 3113, 3123, 3213, 3223

* --**skip-lines**=N allows to skip the first N lines from each LCAO data file
(it affects only LCAO files and not the file with atoms’ positions). Regardless of this value,
extra columns in each input file will be ignored as well.

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

Subdirectory “example” contains a small example consisting of a single wavefunction
The wavefunctions represent the ground state of a single phosphorous dopant
in a small box of silicon atoms with a size of around (5 nm)³.
The example is run as

```
coulombo --atoms=X.3d --dielectric=11.4 --orbitals=10 --skip-lines=1 el_1.dat
```

which will compute a single integral `1111` equivalent to J<sub>ee</sub> electron-electron Coulomb integral.

The results are as follows:

```
 1  1  1  1    0.090340071977
```

which, by taking absolute values of the integrals, translate to

* J<sub>ee</sub> = 90.34 meV

### Licence

This software is provided under a Creative Commons Attribution 4.0 International
Public License, which should be provided as the LICENCE file.
