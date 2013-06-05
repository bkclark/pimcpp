**What is it?**

PIMC++ is a code designed to perform fully-correlated simulations of
quantum systems in continuous space at finite temperature using Path
Integral Monte Carlo. Written in object-oriented C++, it is designed in
a modular way to facilitate easy addition of new methods, observables,
moves and techniques.

**Why does it exist:**

Computational quantum mechanics involves the use of complicated
algorithms and sophisticated codes to calculate properties of quantum
systems. In many situations, people end up reinventing the wheel by
reimplementing algorithms or techniques that already are in existing
codes. PIMC++ is designed to alleviate this problem. It is designed to
be sufficiently general to be of use for a wide variety of Path Integral
calculations, optimized for high performance calculations (both in terms
of speed and parallelization) to allow for quick results, and written in
an object oriented way to facillitate fast prototyping and addition of
new algorithms and techniques. The goal is for this to become a
community-wide code that speeds up the rate of research and that PIMC++
will grow from trhough the contributions of the community. A similar
vision for ground state calculations is manifested in the code QMCPack.

**A note about community involvement:**

PIMC++ has been designed with the goal of being a community-wide code.
Not only does this mean that it is open to users throughout the
community but also that community members can contribute to the
development of the code. This is useful to the general quantum Monte
Carlo community (by developing and maintaining a single code have a wide
variety of research-level features in it) and it is useful to you as a
student as an excellent learning opportunity.

Summer school students are invited to become involved in this process.
To add something to the code, write it, debug it, and contribute it! The
simplest things to add are additional observables for the system and are
a good place to start. Here is a list of possible things you might be
interested in trying your hand at (i.e. exercise for the hard-working
summer school student): (Times are approximate times to add observable
if (you are experienced, you are starting with no knowledge of the code)

-  Lindenmann ratio observable time to complete: (1 hour, 1 day) (easy)
-  Specific Heat observable: (?, ?) (hard)
-  Another observable of your interest

There are many areas in which PIMC++ will benefit from contributions,
including things we haven't even thought of yet! So feel free to talk
with us if you have ideas, applications, or if you would like to be
involved but feel like you need more information first.

Getting started (user)
----------------------

Getting and Installing PIMC++

-  Downloading: Get the latest release from
   `github <https://github.com/etano/pimcpp>`__
-  Installing
-  `User Guide <User Guide>`__
-  `Developer's Guide <Developer's Guide>`__
-  Tutorials

   -  `Particle in a Box <Particle in a Box>`__
   -  `Liquid helium <Bulk helium>`__
   -  `Hydrogen Molecule <H2>`__

-  `PIMC++ applications <PIMC++ applications>`__

