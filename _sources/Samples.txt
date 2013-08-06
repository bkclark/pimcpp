In path integral Monte Carlo (more generically Markov chain Monte Carlo
using the Metropolis move), the probability of accepting a move is
:math:`A_{a\rightarrow b} = \text{min}\left[1, \frac{T_{b\rightarrow a}}{T_{a\rightarrow b}}e^{-\tau \Delta\mathcal S}\right]`
where :math:`\Delta \mathcal S` is the change in action and
:math:`T_{a\rightarrow b}` is the probability of proposing a move to
position b given that you've started in position a.

To understand why the bisection move with free particles (i.e. using
only the kinetic action :math:`\mathcal T`) accepts 100% of the time, we
should understand what :math:`exp^{-\Delta \mathcal T}` is and what
:math:`\frac{T_{b\rightarrow a}}{T_{a\rightarrow b}}` is. To do this,
let's look at the simplest example of a bisection move.

| Consider a single, free particle at three consecutive time slices
along a path, :math:`r_1`, :math:`r_2`, and :math:`r_3`. In this move,
we will keep :math:`r_1` and :math:`r_3` fixed, and will update the
position of :math:`r_2`. We wish to do this optimally, so that it
samples the right distribution given by

.. raw:: html

   <center>

:math:`
\pi(r_2) = \rho_0(r_1, r_2; \tau) \rho_0(r_2,r_3; \tau),
`

.. raw:: html

   </center>

| 
| where :math:` \rho_0` is the free-particle density matrix given by

.. raw:: html

   <center>

:math:`
\rho_0(r,r';\tau) = (4\pi\lambda\tau)^{-\frac{3}{2}} \exp \left[\frac{-(r-r')^2}{4\lambda\tau}\right]
`

.. raw:: html

   </center>

| 
| Writing this product explicitly, we have

.. raw:: html

   <center>

:math:`
\pi(r_2) = (4\pi\lambda\tau)^{-3} \exp \left[-\frac{(r_1-r_2)^2 + (r_3-r_2)^2}{4\lambda\tau}\right]
`

.. raw:: html

   </center>

| 
| Expanding the squares and grouping the :math:`r_2`-dependent terms, we
have

.. raw:: html

   <center>

:math:`
\pi(r_2) = (4\pi\lambda\tau)^{-3} \exp\left[-\frac{(r_1-r_3)^2}{4\lambda\tau}\right]
\exp\left[-\frac{(r_2 - \bar{r})^2}{2\lambda\tau}\right],
`

.. raw:: html

   </center>

where :math:` \bar{r} = (r_1 + r_3)/2 `. We have separated out the
:math:`r_2`-dependence and it is a simple gaussian centered at the
midpoint of the line segment from :math:`r_1` to :math:`r_3`. This is
probably the origin of the name "bisection move". This is special case
of a more general result that the product of two Gaussian distributions
is always another Gaussian distribution.

There exist simple methods to sample a point from a Gaussian
distrubution. We use these methods, then, to randomly choose a value for
:math:`r_2` according to :math:`\pi(r_2)` above.

Image:bisection.png

The new (old) action for this system would be
:math:`exp[-(l_1^2+l_2^2)/4\lambda\tau]` (when we have only the kinetic
action)

:math:`T_{a\rightarrow b}` (the probability of selecting our new
midpoint (the middle bead)) is :math:`exp[-(d_1^2+d_2^2)/4\lambda\tau]`
which when simplified geometrically ends up being
:math:`exp[-(l_1^2+l_2^2-2m^2)/4\lambda\tau]`

When we take the ratio of the new to old, the :math:`m^2` term cancels
out of the probability and we get the whole acceptance probability to be
1.

It is important to sample the kinetic action exactly even when a
potential term is added because failure to do so will cause very bad
acceptance ratios.
