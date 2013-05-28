# PIMC++

Fully-correlated simulations of quantum systems in continuous space at finite temperature using Path Integral Monte Carlo.

## Requirements

* [Blitz++](http://sourceforge.net/projects/blitz/)
* [SPRNG](http://www.sprng.org/)
* [CMake](http://www.cmake.org/)

## Installation

* Set the environmental variables `BLITZ_HOME` and `SPRNG_HOME` to point to the respective installation directories.
* Adjust the CMakeLists.txt file to fit the machine and compilers you are using.
* Run `mkdir build && cd build && cmake ..` to make sure things are configured correctly.
* Run `make` (use `-j` for distributed compiling if possible). Built libraries are put in `$PROJECTDIR/bin/`.


