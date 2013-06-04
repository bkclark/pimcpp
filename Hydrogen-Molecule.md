---
layout: base
title: Hydrogen-Molecule
---

In this tutorial, we will work toward the calculation of a hydrogen
molecule. A key new feature you will explore in this simulation is
PIMC++ ability to deal with long range potentials (Because we are
working with a molecule and most the major difficulties involved in long
range potentials involve periodicity, there is some triviality in
dealing with this. Nonetheless, the techniques you learn here will
immediately generalize to actual periodic system and will allow us to
learn on a real physical system without worrying about fermions yet).

Again, we will start with our free particle in a box and work from
there.

The two key steps in this tutorial will involve changing the input file
for the two species and then adding the required information to get the
coulumb part of the action to work correctly.

Let us start by adding the two relevant species in our simulation,
electrons and protons. In our simulation, we will allow our proton to be
infinitely massive because the quantum effects will be irrelevant on the
energy scales we will be dealing with. Therefore, we will add a Species
Proton (underneath the System section). Because it is infinitely
massive, the proton will be given a \\(\lambda=\frac{\hbar^2}{2m}\\) of 0.
PIMC++ knows that particles with \\(\lambda=0.0\\) are classical and may
treat them specially. We will also add a Species electron (also
underneath the System section). It will have \\(\lambda=0.5\\). PIMC++ has
no units built in and we must simply use consistent units throughout.
For a hydrogen molecule it is convenient to work in Bohr or something.
