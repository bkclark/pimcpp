Particle in a Box
=================

In this tutorial, we calculate the energy of a single particle in a box with periodic boundary conditions. Although not of physical interest, this calculation will introduce you to the PIMC++ software suite, including its input format and output analysis and processing tools.

Specifically, we will calculate the energy of a single particle in a box of size :math:`10 A` at a temperature of :math:`1.666 K`. We are working in a system of units in which energy is in Kelvin, length is in Angstroms, and :math:`\hbar=1`. In these units, the free particle mass will be :math:`0.5 K^{-1} A^{-2}`, which corresponds to :math:`4.034*10^{-26} kg`.

**A note about units:** PIMC++ has no intrinsic units and any system of units can be used by specifying consistent values in the input file. :math:`\beta = M\tau = (k_B T)^{-1} `, where :math:`M` is the number of time slices. :math:`\beta` has units of inverse energy and establishes the temperature scale. For the tutorial, input files will use units of Kelvin and Angstrom with :math:`\hbar = 1`.

In this simple case, it will be possible to compare the simulation with the analytical result. If you do this, you should get an energy of :math:`2.499 K` [#analytical]_.

Input
-----

For more information on each section of the input file, please refer to the documentation for :doc:`Input`.

::

 Section (Output)
 {
   string OutFileBase = "SingleParticle";
 }

 Section (Parallel)
 {
   int ProcsPerClone = 1;
 }

 Section (System)
 {
   double tau = 0.05;
   int NumTimeSlices = 12;
   Array<double,1> Box(3) = [10.0,10.0,10.0];
   Array<bool, 1> IsPeriodic(3) = [true,true,true];
   Section (Particles)
   {
     Section (Species)
     {
       string Name = "free";
       string Type = "free";
       double lambda = 1.0;
       string Statistics = "BOLTZMANNON";
       int NumParticles = 1;
       int NumDim = 3;
       string InitPaths = "BCC";
     }
   }
 }

 Section (Actions)
 {
   int NumImages = 0;
   int MaxLevels = 2;
   Array<string,1> PairActionFiles(1) = ["zero.PairAction"];
 }

 Section (Observables)
 {
   Section (Observable)
   {
     string Type = "Energy";
     string Name = "Energy";
     string Description = "Total Energy";
     int Frequency = 1;
   }
 }

 Section (Moves)
 {

   Section (Move)
   {
     string Type = "BisectionBlock";
     string Name = "BisectionBlock";
     string PermuteType = "NONE";
     string Species = "free";
     int NumLevels = 2;
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
     int Steps = 100;

     Section (Loop){
       int Steps = 100000;

       Section (Move) {string Name = "BisectionBlock";}
       Section (Observe) {string Name = "Energy";}
       Section (Move) {string Name = "Shift";}
     }
     Section (WriteData){}
   }

 }


Running
-------

PIMC++ can be run through the following line

::

    pimc++ SingleParticle.in

In the "Observables" section of the input file, there is variable called "OutFileBase". This specifies the prefix of the output file(s). Because "OutFileBase" (Section: Observables) was set to "SingleParticle" in this case, the output is named "SingleParticle.0.h5". For a serial run, the output filename will be the prefix plus "0.h5". For a parallel run, each processor generates an output file by appending its processor number and ".h5" to the prefix. The output files are written in a portable, hierarchical file format known as `HDF5 <http://hdf.ncsa.uiuc.edu/HDF5/>`__.


Analysis
--------

In the subdirectory pimcpp/Scripts, you will find some analysis tools for parsing output data from PIMC++. These include Analysis.py, which will give averages and errors for every observable, and more specific scripts, such as Energy.py which will output the average energy. Here, we only need the latter to compare to the analytic result.

::

    python Energy.py SingleParticle.0.h5

.. rubric:: Footnotes

.. [#analytical] For a particle in a box with periodic boundary conditions, the quantized energy levels in a cube of side :math:`L` are :math:`E_{ijk} = \frac{4 \hbar^2 \pi^2 (i^2 + j^2 + k^2)}{2mL} = \lambda \frac{2 \pi^2(i^2 + j^2 + k^2)}{L}` where :math:`\lambda=\frac{\hbar^2}{2m}` and :math:`i,j,k = 0,1,2,3...` but :math:`i=j=k=0` is not allowed. The expectation value of the energy is :math:`\frac{1}{Z}\sum_{ijk} \left< \Psi_{ijk}|\hat H \exp[-\beta \hat H]|\Psi_{ijk} \right> = \frac{1}{Z}\sum_{ijk} E_{ijk}\exp[-\beta E_{ijk}]` where :math:`Z=\sum_{ijk} \left< \Psi_{ijk}|\exp[-\beta E_{ijk}]|\Psi_{ijk} \right>`. This sum will converge rapidly to :math:`2.499 K` and can be computed numerically with a few lines of Python code, for example see FG.py in pimcpp/Scripts.
