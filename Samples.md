---
layout: base
title: Samples
---

In path integral Monte Carlo (more generically Markov chain Monte Carlo
using the Metropolis move), the probability of accepting a move is
\<amsmath\>A\_{a\\rightarrow b} = \\text{min}\\left[1,]\</amsmath\>
where \<amsmath\>\\Delta \\mathcal S\</amsmath\> is the change in action
and \<amsmath\>T\_{a\\rightarrow b}\</amsmath\> is the probability of
proposing a move to position b given that you've started in position a.

To understand why the bisection move with free particles (i.e. using
only the kinetic action \<amsmath\>\\mathcal T\</amsmath\>) accepts 100%
of the time, we should understand what \<amsmath\>exp\^{-\\Delta
\\mathcal T}\</amsmath\> is and what \<amsmath\>\\frac{T\_{b\\rightarrow
a}}{T\_{a\\rightarrow b}}\</amsmath\> is. To do this, let's look at the
simplest example of a bisection move.

Consider a single, free particle at three consecutive time slices along
a path, \<amsmath\>r\_1\</amsmath\>, \<amsmath\>r\_2\</amsmath\>, and
\<amsmath\>r\_3\</amsmath\>. In this move, we will keep
\<amsmath\>r\_1\</amsmath\> and \<amsmath\>r\_3\</amsmath\> fixed, and
will update the position of \<amsmath\>r\_2\</amsmath\>. We wish to do
this optimally, so that it samples the right distribution given by\

\<amsmath\> \\pi(r\_2) = \\rho\_0(r\_1, r\_2; \\tau) \\rho\_0(r\_2,r\_3;
\\tau), \</amsmath\>

\
 where \<amsmath\> \\rho\_0\</amsmath\> is the free-particle density
matrix given by\

\<amsmath\> \\rho\_0(r,r';\\tau) = (4\\pi\\lambda\\tau)\^{-\\frac{3}{2}}
\\exp \\left[\\frac{-(r-r')\^2}{4\\lambda\\tau}\\right] \</amsmath\>

\
 Writing this product explicitly, we have\

\<amsmath\> \\pi(r\_2) = (4\\pi\\lambda\\tau)\^{-3} \\exp
\\left[-\\frac{(r\_1-r\_2)\^2] \</amsmath\>

\
 Expanding the squares and grouping the
\<amsmath\>r\_2\</amsmath\>-dependent terms, we have\

\<amsmath\> \\pi(r\_2) = (4\\pi\\lambda\\tau)\^{-3}
\\exp\\left[-\\frac{(r\_1-r\_3)\^2}{4\\lambda\\tau}\\right]
\\exp\\left[-\\frac{(r\_2], \</amsmath\>

where \<amsmath\> \\bar{r} = (r\_1 + r\_3)/2 \</amsmath\>. We have
separated out the \<amsmath\>r\_2\</amsmath\>-dependence and it is a
simple gaussian centered at the midpoint of the line segment from
\<amsmath\>r\_1\</amsmath\> to \<amsmath\>r\_3\</amsmath\>. This is
probably the origin of the name "bisection move". This is special case
of a more general result that the product of two Gaussian distributions
is always another Gaussian distribution.

There exist simple methods to sample a point from a Gaussian
distrubution. We use these methods, then, to randomly choose a value for
\<amsmath\>r\_2\</amsmath\> according to
\<amsmath\>\\pi(r\_2)\</amsmath\> above.

![Image:bisection.png](bisection.png "Image:bisection.png")

The new (old) action for this system would be
\<amsmath\>exp[-(l\_1\^2+l\_2\^2)/4\\lambda\\tau]\</amsmath\> (when we
have only the kinetic action)

\<amsmath\>T\_{a\\rightarrow b}\</amsmath\> (the probability of
selecting our new midpoint (the middle bead)) is
\<amsmath\>exp[-(d\_1\^2+d\_2\^2)/4\\lambda\\tau]\</amsmath\> which when
simplified geometrically ends up being
\<amsmath\>exp[-(l\_1\^2+l\_2\^2-2m\^2)/4\\lambda\\tau]\</amsmath\>

When we take the ratio of the new to old, the
\<amsmath\>m\^2\</amsmath\> term cancels out of the probability and we
get the whole acceptance probability to be 1.

It is important to sample the kinetic action exactly even when a
potential term is added because failure to do so will cause very bad
acceptance ratios.
