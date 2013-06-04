---
layout: base
title: Samples
---

In path integral Monte Carlo (more generically Markov chain Monte Carlo
using the Metropolis move), the probability of accepting a move is
\\(A\_{a\rightarrow b} = \text{min}\left[1,   \frac{T\_{b\rightarrow a}}{T\_{a\rightarrow b}}e^{-\tau   \Delta\mathcal S}\right]\\)
where \\(\Delta \mathcal   S\\) is the change in action and
\\(T\_{a\rightarrow   b}\\) is the probability of proposing a move to
position b given that you've started in position a.

To understand why the bisection move with free particles (i.e. using
only the kinetic action \\(\mathcal T\\)) accepts 100% of the time, we
should understand what \\(exp^{-\Delta   \mathcal T}\\) is and what
\\(\frac{T\_{b\rightarrow   a}}{T\_{a\rightarrow b}}\\) is. To do this, let's
look at the simplest example of a bisection move.

Consider a single, free particle at three consecutive time slices along
a path, \\(r\_1\\), \\(r\_2\\), and \\(r\_3\\). In this move, we will keep \\(r\_1\\) and
\\(r\_3\\) fixed, and will update the position of \\(r\_2\\). We wish to do this
optimally, so that it samples the right distribution given by\

\\(\pi(r\_2) = \rho\_0(r\_1, r\_2; \tau) \rho\_0(r\_2,r\_3;   \tau),\\)

\
 where \\(\rho\_0\\) is the free-particle density matrix given by\

\\(\rho\_0(r,r';\tau) = (4\pi\lambda\tau)^{-\frac{3}{2}} \exp   \left[\frac{-(r-r')^2}{4\lambda\tau}\right]\\)

\
 Writing this product explicitly, we have\

\\(\pi(r\_2) = (4\pi\lambda\tau)^{-3} \exp   \left[-\frac{(r\_1-r\_2)^2 +   (r\_3-r\_2)^2}{4\lambda\tau}\right]\\)

\
 Expanding the squares and grouping the \\(r\_2\\)-dependent terms, we have\

\\(\pi(r\_2) = (4\pi\lambda\tau)^{-3}   \exp\left[-\frac{(r\_1-r\_3)^2}{4\lambda\tau}\right]   \exp\left[-\frac{(r\_2 - \bar{r})^2}{2\lambda\tau}\right],\\)

where \\(\bar{r} = (r\_1 + r\_3)/2\\) . We have separated out the
\\(r\_2\\)-dependence and it is a simple gaussian centered at the midpoint of
the line segment from \\(r\_1\\) to \\(r\_3\\). This is probably the origin of the
name "bisection move". This is special case of a more general result
that the product of two Gaussian distributions is always another
Gaussian distribution.

There exist simple methods to sample a point from a Gaussian
distrubution. We use these methods, then, to randomly choose a value for
\\(r\_2\\) according to \\(\pi(r\_2)\\) above.

![Image:bisection.png](bisection.png "Image:bisection.png")

The new (old) action for this system would be
\\(exp[-(l\_1^2+l\_2^2)/4\lambda\tau]\\) (when we have only the kinetic
action)

\\(T\_{a\rightarrow b}\\) (the probability of selecting our new midpoint (the
middle bead)) is \\(exp[-(d\_1^2+d\_2^2)/4\lambda\tau]\\) which when
simplified geometrically ends up being
\\(exp[-(l\_1^2+l\_2^2-2m^2)/4\lambda\tau]\\)

When we take the ratio of the new to old, the \\(m^2\\) term cancels out of
the probability and we get the whole acceptance probability to be 1.

It is important to sample the kinetic action exactly even when a
potential term is added because failure to do so will cause very bad
acceptance ratios.
