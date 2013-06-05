To have confidence in the accuracy of our simulation results, we should
make sure our simulation is well-converged with respect to the relevant
physical paramaters. For bosons, these will typically be

-  box size (finite size effects): You can simply do larger systems and
   check that observables don't change. Sophisticated methods exist for
   correcting the finite size contributions to the energy but this is
   beyond the scope of this tutorial.
-  time step error: tau should be made smaller until observable
   quantities do not depend on its value. (Note: If you are going to try
   to do this, you should make sure you are working in a bigger box that
   ensures the pair action is 0 for any distance larger than half the
   box length, in contrast to our pathological example. Otherwise, the
   energy is likely to be erratic.)
-  off diagonal terms in the density matrix: In many cases, when we fit
   or calculate the density matrix, we only keep so many terms "off the
   diagonal". Certain observables (especially energy) can be sensitive
   to this.
-  k-points: This is only really applicable when doing simulations with
   long range potentials and beyond the scope of this tutorial.

