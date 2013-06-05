Actions
=======

PIMC++ supports a number of Actions that can be turned on and off
throughout the simulation.

PIMC++ samples the stationary distribution

.. math::

  \mathcal{Z} = \int dR_{0}\dots dR_{M} \exp{(-\sum_{m=0}^{M}\mathcal{S_{m}})}

where :math:`\mathcal{Z}` is the :math:`N`-body density matrix, :math:`M` is the number of time slice divisions, 
and :math:`\mathcal{S_{m}}` is the action on slice :math:`m`. A typical simulation of PIMC++ uses a Kinetic
Action, a Potential Action (usually a pair action of some sort) and
possibly a Nodal Action. Other actions may be used through the code that
are not being sampled in the stationary distribution but are being used
in some other manner. For example, an action might be used at a higher
level in a multistage move to approximate the higher levels.

A typical run of PIMC++ requires certain paramaters set. These
paramaters are set in the Actions section.

Example input:

::

  Section (Actions)
  {
    int NumImages = 1;
    int MaxLevels = 5;
    Array<string,1> PairActionFiles(1) = ["sq.PairAction"];

    Section (Action)
    {
      string Name = "ShortRange";
      string Type = "ShortRange";
    }

    Section (Action)
    {
      string Name = "DavidLongRange";
      string Type = "DavidLongRange";
      double kCutoff = 4.86077676958;
    }
  }

Kinetic Action
--------------

The kinetic action formally is

.. math::

  \mathcal{K_{m}} = \frac{ND}{2}\log{4\pi\lambda\tau} + \sum_{n} [-\frac{R_{m}-R_{m+1}+nL}{4\lambda\tau}]

where :math:`n` is a vector of integer values and :math:`L` is the box length.

NumImages
^^^^^^^^^

This paramater defines whether the sum is only over images in the first box (=0), the
nearest 27 boxes (=1), etc. Typically this can be 0 unless :math:`\lambda\tau` is on the order of the box size :math:`L`.

Potential Action
----------------

The potential action formally is

.. math::

  \mathcal{U_{m}} \equiv \mathcal{S_{m}} - \mathcal{K_{m}}

and can have a variety of formulations. In the simplest formulation, the "primitive approximation",

.. math::

  \mathcal{U_{m}} = \frac{\tau}{2} [V(R_{m})+V(R_{m+1})]

where :math:`V` is the potential. However, this approximation comes with errors which go like :math:`\tau^{2}`.

A more precise approach is to solve the 2-body problem exactly, and then compose the full N-body action as a product of 2-body ones.
This is called the "pair product" approximation, and it is what is used most often in PIMC++.

Array<string,1> PairActionFiles(numPairActions)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specifies the PairActionFiles. The full
path must be given to these files. The PairAction files that can be used
in PIMC++ include those from squarer++ and those from David Ceperely's
squarer. To use the latter requires writing a PairAction file that wraps
the dm files from this squarer. See below for a description to
accomplish this.

Example PairAction file:

::

  Section (Fits)
  {
    string Type="DavidFit";
    int NumOffDiagonalTerms = 3;
    Section (Particle1)
    {
      string Name = "e";
      double lambda =0.0625;
      int Ndim = 3;
    }
    Section (Particle2)
    {
      string Name = "e";
      double lambda =0.0625;
      int Ndim = 3;
    }
    string Daviddmfile = "sq.h5";
  }

.. todo::

  Explain the Squarer process in detail

MaxLevels
^^^^^^^^^

Specifies the total number of levels that
are read in from the density matrix file. This must be less then (or
equal to) the total number of levels in the density matrix. It also must
be more then the number of levels that you are using in bisection (which
is why more then one level is read in). The highest level that is read
in is the largest value of tau being used.
