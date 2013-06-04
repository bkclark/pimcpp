---
layout: base
title: Particle-in-a-Box
---

For a particle in a box with periodic boundary conditions, the quantized
energy levels in a cube of side \\(L\\) are\

\\(E\_{ijk} = \frac{4 \hbar^2 \pi^2 (i^2 + j^2 + k^2)}{2mL} =   \lambda \frac{2 \pi^2(i^2 + j^2 + k^2)}{L}\\)
where \\(\lambda=\frac{\hbar^2}{2m}\\) and \\(i,j,k =   0,1,2,3...\\) but
\\(i=j=k=0\\) is not allowed.\

The expectation value of the energy is\

\\(\frac{1}{Z}\sum\_{ijk} \left< \Psi\_{ijk}|\hat H   \exp[-\beta \hat H]|\Psi\_{ijk} \right> = \frac{1}{Z}\sum\_{ijk}   E\_{ijk}\exp[-\beta E\_{ijk}]\\)
where
\\(Z=\sum\_{ijk}   \left< \Psi\_{ijk}|\exp[-\beta E\_{ijk}]|\Psi\_{ijk}   \right>\\)
.

This sum will converge rapidly to \\(2.51~K\\) and can be computed
numerically with a few lines of Python code, for example.
