For a free particle, the value of the density matrix can be shown to be
:math:`C \exp\left(\frac{-d^2}{4\lambda\tau}\right)` where :math:`d` is
the distance between two particles at adjacent slices. In a periodic
box, this is no longer the case. Instead, the density matrix should be
periodic i.e :math:`\sum_i C \exp[- d(r(0);r_i(\tau))]^2`\ where the
index :math:`i` is the sum over all the periodic boxes. In the situation
where the box length :math:`l \gg \sqrt{4 \lambda \tau}` the former
approximation is a good one. When this is not the case, then the kinetic
density matrix will be incorrect and you will get wrong answers.

To remove this issue, there is a NumImage (Section: Actions) variable in
the input. By setting this variable to 1, it will sum over the images in
the nearest neighbor boxes (all 26 of them). Of course, doing so
significantly slows down the calculation of the kinetic action in the
code and should only be used when you are doing a calculation where it
is critical. To see this effect in action, start making the free
particle box iteratively smaller and compare against the analytical
answer of a free particle in box. Eventually you will find a
disagreement. Then increase NumImages and see the disagreement vanish.
Complete this exercise by lowering the temperature (increasing
:math:`\beta`) so that the system is confined to its ground state.
