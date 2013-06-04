---
layout: base
title: Particle-in-a-Box
---

For a particle in a box with periodic boundary conditions, the quantized
energy levels in a cube of side \<amsmath\>L\</amsmath\> are\
 \<amsmath\>E\_{ijk} = \\frac{4 \\hbar\^2 \\pi\^2 (i\^2 + j\^2 +
k\^2)}{2mL} = \\lambda \\frac{2 \\pi\^2(i\^2 + j\^2 +
k\^2)}{L}\</amsmath\> where
\<amsmath\>\\lambda=\\frac{\\hbar\^2}{2m}\</amsmath\> and
\<amsmath\>i,j,k = 0,1,2,3...\</amsmath\> but
\<amsmath\>i=j=k=0\</amsmath\> is not allowed.\

The expectation value of the energy is\
 \<amsmath\>\\frac{1}{Z}\\sum\_{ijk} \\left\< \\Psi\_{ijk}|\\hat H
\\exp[-\\beta]|\\Psi\_{ijk} \\right\> = \\frac{1}{Z}\\sum\_{ijk}
E\_{ijk}\\exp[-\\beta] \</amsmath\> where \<amsmath\>Z=\\sum\_{ijk}
\\left\< \\Psi\_{ijk}|\\exp[-\\beta]|\\Psi\_{ijk} \\right\>
\</amsmath\>.

This sum will converge rapidly to \<amsmath\>2.51\~K\</amsmath\> and can
be computed numerically with a few lines of Python code, for example.
