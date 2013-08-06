Simple Harmonic Oscillator
==========================

In many circumstances, a harmonic trap is imposed on a system. Here we provide a simple input for free bosons in a harmonic trap of :math:`\omega = 1`.

The temperature of the system is :math:`T = 1/\beta = 1/(\tau*NumTimeSlices) = 1/2 K`. At this temperature, the energy of the system should be :math:`\sim 3.68 K` [#analytical]_.

Input
-----

The input is very similar to that of the :doc:`Particle-in-a-Box` tutorial with the following exceptions:

First the cell is not periodic. This is implied by `IsPeriodic(3) = [false,false,false];`. Because of this the Box size does not matter, so we just default it to `1.0`.

Also, we start with `2` particles with bosonic statistics. This is reflected both in the `Species` section as well as the `PermuteType` and `Gamma` of the `BisectionBlock` move.

Finally, we have added the action `HarmonicPotential` in the `Actions` section with a oscillator frequency of :math:`\omega = 1`.

::

 Section (Output)
 {
   string OutFileBase = "SHO";
 }

 Section (Parallel)
 {
   int ProcsPerClone = 1;
 }

 Section (System)
 {
   double tau = 0.1;
   int NumTimeSlices = 20;
   Array<double,1> Box(3) = [1.0,1.0,1.0];
   Array<bool, 1> IsPeriodic(3) = [false,false,false];
   Section (Particles)
   {
     Section (Species)
     {
       string Name = "e";
       string Type = "e";
       double lambda = 0.5;
       string Statistics = "BOSON";
       int NumParticles = 2;
       int NumDim = 3;
       string InitPaths = "SHO";
     }
   }
 }

 Section (Actions)
 {
   int NumImages = 0;
   int MaxLevels = 3;
   Array<string,1> PairActionFiles(1) = ["zero.PairAction"];

   Section (Action)
   {
     string Name = "HarmonicPotential";
     string Type = "HarmonicPotential";
     string Species = "e";
     double omega = 1.0;
   }
 }

 Section (Observables)
 {
   Section (Observable)
   {
     string Type = "Energy";
     string Name = "Energy";
     string Description = "Total Energy";
     int Frequency = 1;
     double HistStart = 0.0;
     double HistEnd = 1.0;
     int HistPoints = 2;
   }
 }

 Section (Moves)
 {

   Section (Move)
   {
     string Type = "BisectionBlock";
     string Name = "BisectionBlock";
     string PermuteType = "TABLE";
     Array<double,1> Gamma(4) = [1.0,1.0,1.0,1.0];
     Array<string,1> SamplingActions(1) = ["HarmonicPotential"];
     double epsilon = 1e-10;
     string Species = "e";
     int NumLevels = 3;
     int StepsPerBlock = 2;
   }

   Section (Move)
   {
     string Type = "ShiftMove";
     string Name = "Shift";
   }

 }

 Section (Algorithm)
 {

   Section (Loop){
     int Steps = 10;

     Section (Loop){
       int Steps = 100000;
       bool Equilibrate = false;

       Section (Move) {string Name = "BisectionBlock";}
       Section (Observe) {string Name = "Energy";}
       Section (Move) {string Name = "Shift";}
     }
     Section (WriteData){}
   }

 }


Running
-------


Analysis
--------


.. rubric:: Footnotes

.. [#analytical] Free particles in a harmonic trap represents one of the few systems that can be solved analytically through path integral integration. For a script which gives the exact result, please see pimcpp/Scripts/SHO.py
