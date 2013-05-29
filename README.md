# PIMC++

Fully-correlated simulations of quantum systems in continuous space at finite temperature using Path Integral Monte Carlo.

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
* Run `mkdir build && cd build && cmake ..`
* Run `make` (use `-j` for distributed compiling if possible). Built libraries are put in `$PROJECTDIR/bin/`.


