#!/usr/bin/env python

from pylab import *

exact = load ('exact.dat');
intrp = load ('interp.dat');
noint = load ('nointerp.dat');

plot (exact[:,0], exact[:,5],'k',\
      noint[:,0], noint[:,5],'b',\
      intrp[:,0], intrp[:,5],'r')

legend (('Exact', 'Tricubic B-spline', 'Interp. Tricubic B-spline'),loc='upper left');
title ('Laplacian comparison');
savefig ('Laplacians.ps');

figure(2);

plot (noint[:,0], noint[:,5]-exact[:,5], 'b',\
      noint[:,0], intrp[:,5]-exact[:,5], 'r')

legend (('Error in B-spline', 'Error in interpolating B-spline'));
title ('Laplacian error comparison');
savefig ('LaplacianError.ps');

show()
