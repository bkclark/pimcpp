| **Install Libraries using Package Manager**
| Using the package manager, or yum, install the following packages,
| hdf5
| hdf5-devel
| blas-devel
| lapack
| lapack-devel
| gsl
| gsl-devel
| fftw3
| fftw3-devel
| compat-gcc-34-g77
| For the path visualization libraries you will also need,
| gtkmm
| gtkglextmm
| freeglut-devel
| libgnomecanvasmm26-devel
| gtkglext-devel
| python-HTMLgen
| **Install Libraries**
| Download source for the following and install, compiling with your
newly build MPI compilers.
| blitz++: http://www.oonumerics.org/blitz/download/
| sprng-2.0: http://esler.physics.uiuc.edu/downloads/sprng-2.0.tar.gz
| For the path visualization libraries you will also need the following
packages.
| GLE-tubing and extrusion
library\ `1 <http://sourceforge.net/project/showfiles.php?group_id=5539&package_id=5580&release_id=147184>`__
| REVEL-video encoding
library\ `2 <http://sourceforge.net/project/showfiles.php?group_id=118631&package_id=129179>`__
| XVid-xvidcore
libraries\ `3 <http://www.xvid.org/Downloads.43.0.html>`__
| Generally, installation of each of these packages goes like this:
| Unpack the tarball (tar -zxvf myPackage.tar.gz) and cd into the newly
unpack directory
| ./configure --prefix=/usr
| make
| sudo make install
| Make sure your PKG\_CONFIG\_PATH points to all the locations the
"lib".pc files live.
| **Installing PIMC-Common**
| As in the general case,
| autoreconf -i
| ./configure --prefix=/path/to/somewhere --withgsl-libs='-L/usr/lib
-lgsl -lgslcblas -lm'
| make
| make install
| In Fedora the gsl libraries do not come with a .pc file. It is
necessary to add the --with-gsl-libs flag during configuration.
