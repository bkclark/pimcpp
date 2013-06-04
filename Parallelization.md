---
layout: base
title: Parallelization
---

One of the favorable characteristics of Monte Carlo simulation is that
parallelization can be trivial. In the simplest case, "cloning,"
\<amsmath\>M\</amsmath\> independent runs are performed with identical
simulation conditions except that the random number generator has a
unique seed for each clone. Assuming the initial configuration is taken
from an equilibrated distribution and there are not problems ensuring
ergodicity, \<amsmath\>M\</amsmath\> clones will give a factor of
\<amsmath\>M\</amsmath\> more statistical samples. Leaving aside some
technicalities, this means observable quantities will be better
converged by a factor of \<amsmath\>\\sqrt{M}\</amsmath\>.

In the case of Path Integral Monte Carlo, when a relatively low
simulation temperature requires a large number of time slices
\<amsmath\>P\</amsmath\>, another common strategy is time-slice
parallelization. Because the potential action is diagonal in imaginary
time (i.e. only interactions among particles at the same time slice take
place) and the kinetic action is computed between nearest-neighbor time
slices, any time slice only needs information about its two neighbors;
computations on all other time slices can occur independently. To this
end, computations on sets of time slices are split up among the
available processors, with the only overhead being the need to
communicate information about time slices at the boundaries to
neighboring processors. For \<amsmath\>P\</amsmath\> roughly an order of
magnitude larger than \<amsmath\>M\</amsmath\>, this overhead is not
prohibitive. A technical elaboration: for many time slices
\<amsmath\>P\</amsmath\>, it will be necessary to use a so-called
multiple time slice move algorithm like
[Bisection](Moves#Bisection_Block). Here, the kinetic action is computed
not just between nearest-neighbor slices, but between
\<amsmath\>2\^n\</amsmath\>-fold nearest-neighbor slices. In this
realistic case, the need for communication across larger sets of
time-slices is necessary. Hence the rough requirement of \<amsmath\>P
\\approx 10\\cdot M\</amsmath\>.
