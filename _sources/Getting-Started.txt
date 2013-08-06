.. Getting Started

Getting Started
===============

Requirements
------------

PIMC++ requires the following libraries to be installed:

* `FFTW3 <https://github.com/FFTW/fftw3>`_
* `HDF5 <http://www.hdfgroup.org/HDF5/>`_
* `GSL <http://www.gnu.org/software/gsl/>`_
* `SPRNG <http://www.sprng.org/>`_
* `Blitz++ <http://sourceforge.net/projects/blitz/>`_
* `CMake <http://www.cmake.org/>`_

Downloading
-----------

The easiest way to get the newest version of the code is through the github repository:

::

  git clone git://github.com/etano/pimcpp.git

One can also obtain a compressed file of the most recent version of the code `here <http://github.com/etano/pimcpp/tarball/master>`_.

Installation
------------

1. Set the environmental variables to point to the respective installation directories: ``FFTW_HOME``, ``HDF5_HOME``, ``GSL_HOME``, ``SPRNG_HOME``, ``BLITZ_HOME``
2. Adjust the CMakeLists.txt file to fit the machine and compilers you are using.
3. Run ``mkdir build && cd build && cmake ..``.
4. Run ``make`` (use ``-j`` for distributed compiling if possible). Built libraries are put in ``$PROJECTDIR/bin/``.

Usage
-----

To test for a successful build, try out one of the input files from the :doc:`Tutorials`.

For help creating input files, please refer to the :doc:`Input` documentation.

Serial
^^^^^^

If running PIMC++ in serial, simply run

::

  pimc++ /PATH/TO/INPUTFILE.in

Parallel
^^^^^^^^

Parallel execution of pimc++ typically requires ``mpiexec`` or ``mpirun``, where one must specify the precise number of MPI and OpenMP processes.

For example

::

  mpiexec -np 16 -npernode 8 pimc++ /PATH/TO/INPUTFILE.in

This will vary with the specific machine on which the code is being run.

