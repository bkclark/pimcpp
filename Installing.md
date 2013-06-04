---
layout: base
title: Installing
---

### Generic Installation Instructions

\
 PIMC++ has dependencies on a number of open source libraries and relies
on pkg-config to link against these libraries. We use open source make
utilities to construct Makefiles. As such, if you are unfamiliar with
these utilities, installing PIMC++ can seem complicated and tedious.
Nonetheless, it should be possible to install and run on any modern
system running a Unix-based operating system, including Mac OS X and any
Linux flavor.\
 \
 Where we have succeeded in installing PIMC++ on specific architectures,
platforms, and operating systems, specialized instructions are included.
See below if your system is included. We start with a generic set of
instructions meant to work on **any** Unix-based system.\
 \
 **Selecting Compilers**\
 The most important choice to make is whether you intend to run PIMC++
in serial or in parallel. In the latter case, you will need to have an
MPI implementation installed. We recommend Open-MPI
(http://www.open-mpi.org/) or MPICH2
(http://www-unix.mcs.anl.gov/mpi/mpich/).\
 If building in serial, you can use the generic GNU compilers (g++, gcc,
g77) or any specialized compilers for your system (e.g. the Intel
compilers icpc, icc, ifort). If building in parallel you will use the
MPI compilers.\
 \
 Whatever you choose, it is important to use them consistently
throughout the installation process. You can ensure this by setting the
compiler environment variables in your shell of choice. Check this by
issuing\
 echo \$CXX\
 echo \$CC\
 echo \$F77\
 If you need to set any of these variables, follow this example (bash
shell and GNU compilers):\
 export CXX=g++\
 export CC=gcc\
 export F77=g77\
 \
 **Get make utilities**\
 We rely on the following utilities. If you do not already have these
installed (and please ensure you have current versions), links are
provided for download of the source tarballs. As needed, please download
and install.\
 autoconf (http://www.gnu.org/software/autoconf/\#downloading)\
 automake (http://www.gnu.org/software/automake/\#downloading)\
 libtool (http://www.gnu.org/software/libtool/)\
 pkg-config (http://pkgconfig.freedesktop.org/releases/)\
 You will need to define the environment variable PKG\\_CONFIG\\_PATH; see
the next section.\
 \
 **Install Libraries**\
 Download source for the following and install, compiling with your
compilers of choice.\
 **Note**: It is strongly recommended to install all libraries (and
PIMC++) to a single path. For example, if you want everything installed
under /usr (so the executables will be in /usr/bin, libraries in
/usr/lib, etc.), always invoke the configure script with\
 ./configure --prefix=/usr\
 If you install libraries under /usr, you will want to set the
**PKG\\_CONFIG\\_PATH** variable like so:\
 export PKG\\_CONFIG\\_PATH=/usr/lib/pkgconfig\
 **If this variable is not set correctly, pkg-config will not be able to
find the libraries you installed and you will get problems.**\
 Make sure the path where your libraries are installed (e.g. /usr/lib)
is included in **LD\\_LIBRARY\\_PATH**.\
 \
 blitz++:
[http://www.oonumerics.org/blitz/download/](http://www.oonumerics.org/blitz/download/)\
 sprng-2.0:
[http://esler.physics.uiuc.edu/downloads/sprng-2.0.tar.gz](http://esler.physics.uiuc.edu/downloads/sprng-2.0.tar.gz)\
 hdf5:
[http://hdf.ncsa.uiuc.edu/HDF5/release/obtain5.html](http://hdf.ncsa.uiuc.edu/HDF5/release/obtain5.html)\
 fftw3:
[http://www.fftw.org/download.html](http://www.fftw.org/download.html)\
 gsl:
[http://directory.fsf.org/GNUsl.html](http://directory.fsf.org/GNUsl.html)\
 \
 Generally, installation of each of these packages goes like this:\
 Unpack the tarball (tar -zxvf myPackage.tar.gz) and cd into the newly
unpacked directory\
 ./configure --prefix=/usr\
 make\
 sudo make install (requires root access if you want to install under
/usr)\
 \
 \
 **Get and Install Common**\
 PIMC++ links against the Common library. Download a tarball here:\

[http://pathintegrals.info/oldIndex.html](http://pathintegrals.info/oldIndex.html)\
 Unpack it, as above.\
 cd to the top-level directory and generate the configure script:\
 autoreconf --install\
 If your installation of libraries went ok, the pkg-config utility
should be able to find everything it needs. However, this is where
idiosyncrasies of individual machines and setups can cause problems. To
do some preemptive troubleshooting, you should check that pkg-config can
find what it's going to look for. Issue the following commands and
compare the output against the typical output in parenthesis:\
 pkg-config blitz --libs (-L/usr/lib -lblitz)\
 pkg-config sprng2 --libs (-L/usr/lib -lsprng)\
 pkg-config gsl --libs (-L/usr/lib -lgsl -lgslcblas -lm)\
 pkg-config fftw3 --libs (-L/usr/lib -lfftw3 -lm)\
 If you installed somewhere other than /usr, the path shown will be
different. But if your output for any of these is dramatically different
(most importantly, if the path is **not** where your library is
installed **or** if the linker flags are not as shown), something has
gone wrong with that library installation, or pkg-config is not set up
correctly.\
 \
 Assuming everything looks ok, go ahead and configure:\
 **serial version:** ./configure --prefix=/usr\
 **parallel version:** ./configure --prefix=/usr --enable-mpi\
 Again, you can install under another path; just substitute for /usr.\
 If you installed the libraries in an unusual place, or if pkg-config is
not set up correctly, the configure script will complain that it can't
find certain libraries (often, it will not find where you installed the
hdf5 libraries). When this happens, it will ask you to explicitly give
the path (and linker flags) by specifying another configure flag. See
./configure --help for a list of the flags you can specify.\
 \
 If configure runs without error, it will generate the necessary
Makefiles. Now, issue\
 make\
 make install (may require root access)\
 \
 To test your build of common, try\
 pkg-config common --libs (serial)\
 pkg-config mpicommon --libs (parallel)\
 What you should see is an aggregation of all the paths and flags that
Common needed (included the pkg-config output from each library you
already installed), but the specific output will depend on your system
and installation options. If something looks amiss, try to identify the
source of problem now to save yourself some time and frustration.\
 \
 **Install PIMC+++**\
 Good! You've made it though the hard part.\
 Download the PIMC++ tarball here:\

[http://pathintegrals.info/oldIndex.html](http://pathintegrals.info/oldIndex.html)\
 Unpack it, as above.\
 cd to the top-level directory and generate the configure script:\
 autoreconf --install\
 Again, be sure that pkg-config can find the Common installation you
just completed:\
 pkg-config common --libs (serial)\
 pkg-config mpicommon --libs (parallel)\
 As with Common, you need to build either in serial or in parallel. The
minimal configure options go:\
 **serial version:** ./configure --prefix=/usr\
 **parallel version:** ./configure --prefix=/usr --enable-mpi\
 Add any other options you may want (do ./configure --help to see the
possibilities).\
 If this worked, then do\
 make\
 make install\
 You should now have the executable, pimc++, in the src directory and
installed wherever you specified. Time to run some simulations!\
 \
 **It Didn't Work!**\
 If you ran into trouble, consult the list below to see if there are
specialized installation instructions for your setup. If not, please
email the developers with a description of your problem. To facilitate
communication, tell us about your system (minimally chip architecture,
operating system), whether you are attempting to build the serial or
parallel version, and error output from configure or make and anything
else you think is relevant.

* * * * *

### Specialized Instructions and Known Issues

Listed by system specs\
 [Ubuntu on 64-bit Intel Core 2 Quad](quad)\
 [Fedora Core 7 on 32-bit AMD with GNU compilers](FC7)\

