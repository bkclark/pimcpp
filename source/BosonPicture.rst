Below we see a cartoon of three closed loops (i.e particles) in 2
dimensions. The numbers represent different points in imaginary time
(i.e :math:`\tau`,\ :math:`2\tau`, and :math:`3\tau`.) For
distinguishable particles in path integral monte carlo, you sample over
all the positions (of all the beads) of these closed loops.
Image:1perm.png

For Bosons, you sample over everything above as well as all the
positions (of all the beads) of the small loops permuted onto each other
to form bigger loops:

| Here we have two particles permuting onto each other to form a bigger
loop:
|  Image:2perm.png
| Here we have all three particles permuting onto each other:
|  Image:3perm.png

This picture helps make it clear when statistics are important in a
quantum system (as opposed to just zero point motion). Whenever the
particles are close enough such that they can permute onto each other
without stretching their connection (springs) out too much, then
statistics becomes important. The length scale for a spring's
displacement from equilibrium is :math:`\sqrt{4 \lambda \tau}`, which
should be compared with the average interparticle separation,
:math:`r_s`, which is defined by
:math:`\frac{4}{3}\pi r_s^3 = \frac{V}{N}`, where the total system
volume and particle number are used.
