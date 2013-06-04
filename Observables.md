---
layout: base
title: Observables
---

To add an observable to your calculation, there are three steps:

1.  Define the observable.
2.  Tell the Algorithms section when to evaluate the observable.
3.  Tell the Algorithms section when to write out the observable. Note:
    ALL observables (except for PathDumps) are written out at the same
    time.

To begin with, one needs to define the observable in the observables
section and then one needs to call the observable in the algorithm.

Table of Contents

-   -   [Defining the observable](#Defining_the_observable)
    -   [Inserting the observable into the
        algorithm](#Inserting_the_observable_into_the_algorithm)
    -   [Telling the observables to write to the output
        file](#Telling_the_observables_to_write_to_the_output_file)
    -   [Paramaters](#Paramaters)

-   [List of Observables](#List_of_Observables)
    -   [Autocorrelation](#Autocorrelation)
    -   [Energy](#Energy)
    -   [Forces](#Forces)
    -   [Grids](#Grids)
    -   [Pair Correlation Function](#Pair_Correlation_Function)
    -   [n(r)](#nr)
    -   [Particle Average Location](#Particle_Average_Location)
    -   [PathDump](#PathDump)
    -   [Permutation Counting](#Permutation_Counting)
    -   [Pressure](#Pressure)
    -   [Structure Factor](#Structure_Factor)
    -   [Superfluid Fraction](#Superfluid_Fraction)
    -   [Time Analysis](#Time_Analysis)
    -   [Vacancy Location](#Vacancy_Location)
    -   [Weight](#Weight)
    -   [Winding Number](#Winding_Number)

Defining the observable
----------------------------------------------------------------------------------------------------------------

All observables are defined in their respective section *Observable*
within the overall section *Observables*

     Section(Observables)
     {
       Section (Observable)
       {
         string Type="Energy";
         string Name="He-He_Energy";
       }
       Section (Observable)
       {
         ...
       }
     }

Inserting the observable into the algorithm
--------------------------------------------------------------------------------------------------------------------------------------------------------

A sample algorithm looks like:

     Section (Algorithm)
     {
       Section (Loop){ //Start Accumulating observables
         int Steps=3000000000;
         Section (Loop){
           int Steps=100;
           Section (Move) {string Name="BisectionMoveForHelium";}
           Section (Loop){
             int Steps=10;
             Section (Move) {string Name="Displace";}
           }
           Section (Move) {string Name="CenterOfMass";}
           Section (Observe) {string Name = "Vacancy"; }
           Section (Observe) {string Name = "PathDump"; }
           Section (Observe) {string Name = "StructureFactor"; }
           Section (Observe) {string Name = "TimeAnalysis"; }
           Section (Observe) {string Name = "CycleCount"; }
           Section (Move) {string Name = "Shift"; }
         }
         Section (WriteData){}
       }
     }

To add your move to the algorithm, you must place in a line like:

      Section (Observe) {string Name="Vacancy";}

where the name is the name you have chosen for your observable.

Telling the observables to write to the output file
------------------------------------------------------------------------------------------------------------------------------------------------------------------------

All observables (except PathDump) are written when the following line is
reached in the algorithm:

      Section (WriteData){}

Paramaters
--------------------------------------------------------------------------------------

All observables have (at minimum) the variables *Type*, *Name*, and
*Frequency* which must be specified. The variable *Type* tells the code
which observable to use (i.e. energy, paircorrelation, etc). You must
use a string picked from the list of observables. The variable *Name* is
a name you choose for the observable to be used later on when you call
it in the algorithm. For example, your Name might be
"EnergyForTheUpElectrons". The variable *Frequency* tells the observable
how often it accumulates (NOT writes) its running average. Every
Frequency number of iterations through the algorithm, the observable is
accumulated to the running average.

All observables also have a number of other paramaters that need not be
specified. They are:

1.  Description: A description of the observable that gets placed in the
    output file solely for the sake of the user to read.

Many observables have other paramaters. For example, the pair
correlation function has an observable that indicates what two species
it works on. Inside its observable block, there might be

     Section (Observable)
       {
         string Type="PairCorrelation";
         string Name="He-He3_PC";
         string Species1="He";
         string Species2="He3";
       }

In cases where you are just editing the input to an observable instead
of adding one, you will need only to change the paramaters of the
observable.

* * * * *

List of Observables
========================================================================================================

Note: All observables include a paramater Name that is user specified.

Autocorrelation
------------------------------------------------------------------------------------------------

This is a specialized observable to compute the dipole-dipole
autocorrelation time for molecules such as water. It will be generalized
for other systems. However, at present, using this observable is not
recommended for any system other than rigid five-site water models.

*Example input:*

      Section (Observable)
        {
            string Type = "AutoCorr";
            string Name = "AutoCorrelation";
            int Frequency=1250;
            int dumpFrequency=50; // #steps = freq*dumpFreq
            int numSlots=50;      // should ensure dumpFrequency >= numSlots

            Section (Grid)
              {
                string type = "Linear";
              }
          }

**Paramaters:**

**Type:** AutoCorr

**Frequency:** As for all observables, Frequency gives the interval of
Monte Carlo steps between measurements of the instantaneous dipole
moments of the molecules.

**dumpFrequency:** Specifies the number of measurements to accumulate
(not the number of MC steps) before the autocorrelation time is
computed. So, the number of MC steps over which the autocorrelation time
is computed will be Frequency\*dumpFrequency.

**numSlots:** Specifies the number of bins in which dipole moments are
accumulated. So the maximum possible autocorrelation time (in MC steps)
is given by Frequency\*numSlots. It is practically important to set
dumpFrequency\>numSlots

Energy
------------------------------------------------------------------------------

The energy observable calculates the energy by evaluating the average
beta derivative of the action.

Example input: Input Paramaters: Type: Energy

ComputeEnergies: This Array of strings is an optional parameter used to
specify additional Action objects that can compute a system energy but
are not computed by default. If the strings specified are recognized as
valid Action objects, they will be computed whenever the Energy
observable is called and will be included in the sum of all energies
reported as Total in the [output]. Currently, the supported Action
objects that can be specified are:

[ST2WaterClass]

[QMCSamplingClass]

[IonIonActionsClass]

[CEIMCActionClass]

Forces
------------------------------------------------------------------------------

This observable generates a distribution of the force on a given
species. It has not been widely tested and should be regarded as under
development. *Example input:*

      Section (Observable)
          {
              string Type = "Forces";
              string Name = "ForceOnNa";
              int Frequency = 10;
              string Species = "Na";
           }

*Input Parameters:*

**Type:** Forces

**Species:** Specify the species on which to compute forces

Grids
----------------------------------------------------------------------------

Grid is a separate Section that is contained within a number of
observables, but it is not an observable *per se*.

*Example input:*

      Section (Grid)
            {
              string Type = "Linear";
              double start = 0.0;
              int NumPoints = 100;
            }

*Input Paramaters:*

-   **Type:** Linear
-   **start:** An optional parameter specifying the start of the grid
    (default is 0.0)
-   **end:** An optional parameter specifying the end of the grid
    (default is the box size for a periodic system)
-   **NumPoints:** A required parameter specifying the number of grid
    points (i.e. the resolution) of the grid

Pair Correlation Function
--------------------------------------------------------------------------------------------------------------------

The pair correlation function calculates \<amsmath\>g(r)\</amsmath\>

*Example Input:*

      Section (Observable)
          {
            string Type = "PairCorrelation";
            string Name = "HeHePC";
            string Species1 = "He";
            string Species2 = "He";
            string Description="Helium-Helium Pair Correlation";
            int Frequency=2;
            Section (Grid)
              {
                string Type = "Linear";
                double start = 0.0;
                int NumPoints = 100;
              }
          }

*Input Paramaters:*

-   **Type:** PairCorrelation
-   **Species1:** One species to be used in pair correlation function
-   **Species2:** Second species to be used in pair correlation
    function.
-   **Description:** Description of the observable
-   **Grid Section:** See section in observables on defining
    [Grids](#Grids)
-   **Frequency:** Observable accumulated every Frequency number of
    times it is called in the algorithm.
-   **Name:** User specified and consistent throughout input file.

n(r)
------------------------------------------------------------------------

Particle Average Location
--------------------------------------------------------------------------------------------------------------------

PathDump
----------------------------------------------------------------------------------

Every Frequency steps, the particle paths and the permutations are
written to the output file in .h5 format. If you wish to restart from an
old configuration, it is important that you have included a PathDump in
your simulation (otherwise there will be no data from which to restart)

**Example Input:**

      Section (Observable)
      {   
          string Type="PathDump";
          string Name="PathDump";
          int Frequency=2;
      }

*Input Parameters:*

-   **Type:** PathDump
-   **Name:** User defined that is consistent throughout the algorithm
-   **Frequency:** DIFFERENT then all other observables. Every frequency
    number of times pathdump is called in the algorithm it WRITES to a
    file. All other observables, accumulate every frequency number of
    calls but only write at a WriteBlock. PathDump does nothing at a
    WriteBlock.
-   **AllClones:** boolean variable. If True, then all the parallel
    clones will dump the path. If False, then only the 0'th clone will
    dump the path.

Permutation Counting
----------------------------------------------------------------------------------------------------------

Pressure
----------------------------------------------------------------------------------

*Example input:*

      Section (Observable)
          {
            string Type = "Pressure";
            string Name = "Pressure";
            string Description="Total Pressure";
            int Frequency=2;
            double Prefactor=138.065;
          }

*Input Paramaters:* **Type:** Pressure **Prefactor:** This is a number
that the pressure can be multiplied by. Typically this is chosen to
change the pressure into more sane units. the number 138.065 changes the
pressure into bars.

*Output Paramaters:* **ShortRangePressure:** Specifies the component of
the pressure that comes from the shrot range ...

Structure Factor
--------------------------------------------------------------------------------------------------

This observable computes the Structure Factor
\<amsmath\>\\frac{1}{\\sqrt{N\_a
N\_b}}\\rho\_{k,a}\\rho\_{-k,b}\</amsmath\> between species A and B

Example input:

    Section (Observable)
         {
           string Type="StructureFactor";
           string Name="HeliumStructureFactor";
           double kCutoff=5.0;
           int Frequency=1;
           string Species1="Helium4";
           string Species2="Helium4";

         }

*Input Paramaters:*

-   **Type:** StructureFactor
-   **Name:** Anything user defined that is consistent throughout the
    input
-   **kCutoff:** maximum k that should be included in the structure
    factor. If the long range action is being used, this MUST be the
    same cutoff as used in the long range action. If you wish to include
    other k-vectors, use the AdditionalkVecs variable (see below)
-   **Frequency:** Observable accumulated every Frequency number of
    times it is called in the algorithm
-   **Species1** and **Species2**: The names of the respective species
    used in calculating the structure factor.
-   **Array\<double,2\>\</double,2\>
    AdditionalkVecs(number\_of\_k\_vecs,NDIM):** This is a list of
    additional kvectors that are to be included beyond those that are
    defined by the kCutoff.

*Output Paramaters:*

Superfluid Fraction
--------------------------------------------------------------------------------------------------------

The superfluid fraction calculates
\<amsmath\>\\frac{\\rho\_s}{\\rho}=\\frac{\\left\<W\^2\\right\>}{2\\lambda\\beta
N}\</amsmath\>

*Example Input*:

    Section (Observable)
    &#123;
       string Type&#61;&quot;SuperfluidFraction&quot;&#59;
       string Name&#61;&quot;SuperfluidFraction&quot;&#59;
       Array&lt;string,1&gt; SpeciesList(1)&#61;&#91;&quot;He&quot;&#93;&#59;
       int Frequency&#61;1&#59; 
    &#125;

*Input Paramaters:*

-   **Type:** SuperfluidFraction
-   **Name:** User specified and consistent throughout input file.
-   **Array\<string,1\>\</string,1\> SpeciesList(num\_of\_species):**
    List of species to calculate the superfluid fraction on (currently
    this doesn't work. It calculates it on the entire system independent
    of this variable)
-   **Frequency:** Observable accumulated every Frequency number of
    times it is called in the algorithm

*Output Paramaters:*

Time Analysis
--------------------------------------------------------------------------------------------

The time analysis measures how long spend in different sections of the
code. (note: if you suspend you're code, this observable will give
extremely erroneous results (the time measurements are depending on the
wallclock time)

*Example Input:*

    Section (Observable)
    &amp;&#35;123&#59;
      string Type&amp;&#35;61&#59;&amp;quot&#59;TimeAnalysis&amp;quot&#59;&amp;&#35;59&#59;
      string Name&amp;&#35;61&#59;&amp;quot&#59;TimeAnalysis&amp;quot&#59;&amp;&#35;59&#59;
      int Frequency&amp;&#35;61&#59;1&amp;&#35;59&#59;
    &amp;&#35;125&#59;

*Input Paramaters:*

-   **Type:** TimeAnalysis
-   **Name:** User specified and consistent throughout input file.
-   **Frequency:** Needs to be specified but not actually used

Vacancy Location
--------------------------------------------------------------------------------------------------

Weight
------------------------------------------------------------------------------

Winding Number
----------------------------------------------------------------------------------------------

*Example input:*

      Section (Observable)
        {
          string Type="WindingNumber";
          string Name="WindingNumber";
          Array&lt;string,1&gt;&lt;/string,1&gt; SpeciesList(1) = ["He"];
          string Description="Winding Number";
          int Frequency=2;
          int dumpFrequency=20;
          double kCutoff=2.55;
        }

*Input Paramaters:*

**Type:** WindingNumber

*Output Paramaters:*

     Specifies the component of the pressure that comes from the shrot range ...
