| For a particle in a box with periodic boundary conditions, the
quantized energy levels in a cube of side :math:`L` are
| :math:`E_{ijk} = \frac{4 \hbar^2 \pi^2 (i^2 + j^2 + k^2)}{2mL} = \lambda \frac{2 \pi^2(i^2 + j^2 + k^2)}{L}`
where :math:`\lambda=\frac{\hbar^2}{2m}` and :math:`i,j,k = 0,1,2,3...`
but :math:`i=j=k=0` is not allowed.
|  The expectation value of the energy is
| :math:`\frac{1}{Z}\sum_{ijk} \left< \Psi_{ijk}|\hat H \exp[-\beta \hat H]|\Psi_{ijk} \right> = \frac{1}{Z}\sum_{ijk} E_{ijk}\exp[-\beta E_{ijk}] `
where
:math:`Z=\sum_{ijk} \left< \Psi_{ijk}|\exp[-\beta E_{ijk}]|\Psi_{ijk} \right> `.

This sum will converge rapidly to :math:`2.51~K` and can be computed
numerically with a few lines of Python code, for example.
