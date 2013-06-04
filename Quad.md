---
layout: base
title: Quad
---

### Installing on 64-bit Ubuntu (quad core) with gnu compilers

add adept-manager software:

using adept-manager get: autoconf\
 automake libtool libhdf5-mpich-dev (mpich-compatible hdf5 library)
libmpich-mpd1.0c2 libmpich1.0-dev libmpich1.0c2 mpich-bin mpich-mpd-bin

blitz++ (must also move blitz.pc to /usr/lib/pkgconfig) fftw3-dev

g77 libgsl0-dev svn

**Installing spring tarball:** ./configure --prefix=/usr make sudo make
install

Check out Common: svn co
svn+ssh://esler.physics.uiuc.edu/home/svn/repos/Common/trunk/ Common/

autoreconf --install

(For distcc compilation modify quad.conf with distcc in front of CXX)
./quad.conf

make or make -j8 if using distcc

sudo make install

Check out and install PIMC++: svn co
svn+ssh://esler.physics.uiuc.edu/home/svn/repos/PIMC++/trunk/ PIMC++/
autoreconf --install ./configure --prefix=/usr --enable-mpi
CXXFLAGS=-DMPICH\_IGNORE\_CXX\_SEEK CXX="distcc mpicxx" CC=mpicc
F77=mpif77 make -j8 export
LD\_LIBRARY\_PATH="/usr/lib:/opt/intel/mkl/9.1/lib/em64t/"

Not for PIMC++ but useful cmake distcc

* * * * *

### Installing on 64-bit Ubuntu (quad core) with Intel Fortran compilers

\
 \
 The key here is **not** to rely on package manager software to deliver
libraries because they will likely be built with GNU compilers.\
 \
 The only guaranteed way to avoid this is to do the following steps.\
 \
 **Install Intel C++ and Fortran compilers**\
 These can be downloaded free for non-commercial use from Intel:\

[http://www.intel.com/cd/software/products/asmo-na/eng/340679.htm](http://www.intel.com/cd/software/products/asmo-na/eng/340679.htm)\
 On ubuntu, the best way to do this is to run the script (ignore the
warnings that it can't find the 32 bit libraries)\
 The Intel Math Kernel Libraries (mkl) can also be installed this way.\
 Once installed, it is necessary to add the bin directories to the
\$PATH environment variable and the lib directories to
\$LD\_LIBRARY\_PATH.\
 Also, ensure that you are using the intel compilers at all times (\$CXX
= icpc, \$CC=icc, \$F77=ifort).\
 \
 **Get make utilities**\
 Using adept-manager get:\
 autoconf\
 automake\
 libtool\
 \
 if you have not installed these already (it's fine if they're built
with GNU).\
 \
 **Install an MPI Implementation**\
 We have had problems with MPICH2, although it is probably fine if you
can get it working. Open-MPI is known to work on this platform and with
this software.\
 Install according to standard procedure, being sure to build with
**Intel compilers**!\
 \
 **Install Libraries**\
 Download source for the following and install, compiling with your
newly build MPI compilers.\
 **Note**: To be sure that your MPI compilers (mpic++, mpicc, mpif77)
wrap the Intel compilers, do\
 mpic++ --version\
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
unpack directory\
 ./configure --prefix=/usr\
 make\
 sudo make install\
 \
 \

After installing blitz there will be a file blitz-uninstall.pc (usually
in /usr/lib/pkgconfig/). This has to be removed (it confuses pkg-config
because it thinks those flags are relevant when being asked for blitz)\
 \
 Get svn if you don't have that already (package manager should be able
to do that)\
 Check out Common:\
 svn co svn+ssh://esler.physics.uiuc.edu/home/svn/repos/Common/trunk/
Common/\
 autoreconf --install\
 (For distcc compilation modify quad.conf with distcc in front of CXX)\
 ./quad.conf\
 make or make -j8 if using distcc\
 sudo make install\
 \
 \
 Check out and install PIMC++:\
 svn co svn+ssh://esler.physics.uiuc.edu/home/svn/repos/PIMC++/trunk/
PIMC++/\
 autoreconf --install\
 ./configure --prefix=/usr --enable-mpi
CXXFLAGS=-DMPICH\_IGNORE\_CXX\_SEEK CXX="distcc mpic++" CC=mpicc
F77=mpif77\
 make -j8\
 export LD\_LIBRARY\_PATH="/usr/lib:/opt/intel/mkl/9.1/lib/em64t/"\

* * * * *
