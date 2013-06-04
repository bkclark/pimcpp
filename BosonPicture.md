---
layout: base
title: BosonPicture
---

Below we see a cartoon of three closed loops (i.e particles) in 2
dimensions. The numbers represent different points in imaginary time
(i.e \<amsmath\>\\tau\</amsmath\>,\<amsmath\>2\\tau\</amsmath\>, and
\<amsmath\>3\\tau\</amsmath\>.) For distinguishable particles in path
integral monte carlo, you sample over all the positions (of all the
beads) of these closed loops.
![Image:1perm.png](1perm.png "Image:1perm.png")

For Bosons, you sample over everything above as well as all the
positions (of all the beads) of the small loops permuted onto each other
to form bigger loops:

Here we have two particles permuting onto each other to form a bigger
loop:\
 ![Image:2perm.png](2perm.png "Image:2perm.png")\
 \
 \
 \
 Here we have all three particles permuting onto each other:\
 ![Image:3perm.png](3perm.png "Image:3perm.png")

This picture helps make it clear when statistics are important in a
quantum system (as opposed to just zero point motion). Whenever the
particles are close enough such that they can permute onto each other
without stretching their connection (springs) out too much, then
statistics becomes important. The length scale for a spring's
displacement from equilibrium is \<amsmath\>\\sqrt{4 \\lambda
\\tau}\</amsmath\>, which should be compared with the average
interparticle separation, \<amsmath\>r\_s\</amsmath\>, which is defined
by \<amsmath\>\\frac{4}{3}\\pi r\_s\^3 = \\frac{V}{N}\</amsmath\>, where
the total system volume and particle number are used.
