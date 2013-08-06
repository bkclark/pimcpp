# PIMC++

PIMC++ is a code designed to perform fully-correlated simulations of quantum systems in continuous space at finite temperature using Path Integral Monte Carlo.
Written in object-oriented C++, it is designed in a modular way to facilitate easy addition of new methods, observables, moves and techniques.

Computational quantum mechanics involves the use of complicated algorithms and sophisticated codes to calculate properties of quantum systems.  In many situations, 
people end up reinventing the wheel by reimplementing algorithms or techniques that already are in existing codes. 
PIMC++ is designed to alleviate this problem. It is designed to be sufficiently general to be of use for a wide variety 
of Path Integral calculations, optimized for high performance calculations (both in terms of speed and parallelization) 
to allow for quick results, and written in an object oriented way to facillitate fast prototyping and addition of 
new algorithms and techniques.  The goal is for this to become a community-wide code that speeds up the rate of research 
and that PIMC++ will grow from through the contributions of the community.

## Requirements

* [FFTW3](https://github.com/FFTW/fftw3)
* [HDF5](http://www.hdfgroup.org/HDF5/)
* [GSL](http://www.gnu.org/software/gsl/)
* [SPRNG](http://www.sprng.org/)
* [Blitz++](http://sourceforge.net/projects/blitz/)
* [CMake](http://www.cmake.org/)

## Installation

* Set the environmental variables to point to the respective installation directories:
  `FFTW_HOME`,`HDF5_HOME`, `GSL_HOME`, `SPRNG_HOME`, `BLITZ_HOME`
* Adjust the CMakeLists.txt file to fit the machine and compilers you are using.
* Run `mkdir build && cd build && cmake ..`.
* Run `make` (use `-j` for distributed compiling if possible). Built libraries are put in `$PROJECTDIR/bin/`.

## Running

Basic execution is done very simply by `pimc++ InputFile.in`.

## Further information

For help on any subject (including creating input files), please see the [documentation](http://etano.github.io/pimcpp).


