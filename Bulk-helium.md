---
layout: base
title: Bulk-helium
---

**Liquid Helium Tutorial**

In this tutorial, we will estimate the transition temperature between
normal liquid Helium and superfluid liquid Helium at saturated vapor
pressure (SVP). To estimate the superfluid transition temperature, we
need to scan over a range of temperatures to identify the critical
point.

**Warning:** Due to time constraints in this tutorial, our calculation
will use a small number of particles (10) and a correspondingly small
simulation box. Consequently, we will be violating the rule that
short-range pair interactions must be zero beyond half the box length
(they should be truncated and tapered to ensure a smooth potential
energy surface). In other words, our simulation box is smaller than the
length scale of the particle pair interactions, which introduces
discontinuities into the potential energy surface. Thus, this system is
too small to give accurate results and is intended as a pedagogical
example only.

You are encouraged to prepare your own input file for these simulations
by converting the single particle input file we just used. By doing so,
you will learn a lot about the capabilities and the options of the code.
Instructions for doing this are contained
[here](WritingHeliumInputfile). However, this process will be
time-consuming for first-time users, so we will provide the relevant
input files below so you can get started with the runs right away. If
you choose to move on and use the prepared input files for your runs,
you should come back and review [the detailed
instructions](WritingHeliumInputfile) at the end.

Table of Contents

-   [Understanding quantum effects](#Understanding\_quantum\_effects)
-   [Superfluid Transition](#Superfluid\_Transition)
    -   [Cloning Parallelization](#Cloning\_Parallelization)
    -   [Time slice Parallelization](#Time\_slice\_Parallelization)

-   [Superfluid Transition Temperature Calculation (larger
    system)](#Superfluid\_Transition\_Temperature\_Calculation\_larger\_system)

Understanding quantum effects
----------------------------------------------------------------------------------------------------------------------------

In this section, we will explore the nature of quantum effects and how
they are manifested in PIMC simulations. As we lower the temperature of
Helium-4, we will see a transition from particles that are effectively
classical, to particles where quantum zero point motion becomes
important, and finally to particles where bosonic exchange statistics
become important.

We will start with the input file contained here:
[Helium.in](Helium\_File\_2). If you generated your own input, compare it
against this example to check for typos. You also need
He4.95.OffDiag.dm, He4.95.Sampling.in and He4.95.PairAction; copies of
these files are in a location that the setup\\_PIMC script will find so
you don't need to worry about them if you are running on a lab machine.

Now, we want to start with liquid helium in the high temperature regime.
We will start by examining a temperature of 64K. For this temperature we
will use tau=0.0015625 with NumTimeSlices=10. To prepare this run, open
the input file for editing and change tau (Section: System) and
NumTimeSlices (Section: System) appropriately. Run pimc++ with

    pimc++ Helium.in

This takes some time. You may want to look over how to write the [file
input](WritingHeliumInputfile) file, generate the [pair density
matrix](Constructing%20a%20density%20matrix), or how you can[explicitly
integrate over permutations](BosonPicture).

Now that the simulation has completed, we can analyze its results. To do
this, run the report generator (Report.py) on the h5 file that it
produced and examine the output with your web browser. As before, you
can use firefox or konqueror to view this output. It should look very
similar to
[this.](http://esler.physics.uiuc.edu/tutorial/liquidHelium/Helium64/)

Notice that not only is there information about the energy but we also
see the structure factor, superfluid fraction, and pair correlation
function. Looking at the structure factor, we can see that the system we
are simulating is a liquid (the structure factor of a solid would have a
noticeable peak.) The value of \\(\frac{\rho\_s}{\rho}\\) is 0 indicating
that the system is not superfluid. From the pair correlation, we can
notice a couple things. To begin with, we clearly see that Helium has a
hard core (there are effectively no particles within a distance of 2).
We also can tell that the well of the inter-particle potential is
approximately where the peak of the pair correlation function is shown.

The analysis should give a summary of all the observables except the
PathDump. The PathDump observable dumps complete configurations of the
path so that one can look at visual images of them. This will allow us
to get physical intuition about how quantum effects become important. On
the lab machines, we will use the pathvis++ tool set up by the
\\(setup\_pimc++\\) script. Type

    pathvis++ Helium.0.h5

We now see a screen like this:\

**^4^He at 64 K**
![Image:Helium.64.png](Helium.64.png "Image:Helium.64.png")

The simulation box can be rotated by left-clicking and dragging the
mouse. You may zoom with the scroll wheel or by middle-clicking and
dragging left (zoom out) or right (zoom in). At the bottom of the screen
is a slider that allows us to change through frames of each successive
Monte Carlo step that was dumped to the output file. You should notice
here that the particle paths are are much smaller than the
inter-particle spacing. This compactness essentially indicates that at
all points through imaginary time, the particle is in nearly the same
place. In this regime, the helium atoms would be reasonably well
approximated by treating them as classical particles.

Let's now lower the temperature and see what happens when quantum
effects become important. We will lower the temperature to approximately
4.2 K. Change tau (Section: System) to 0.0125 and find an appropriate
value to change NumTimeSlices (Section: System) to in order to get the
desired temperature. (Caveat: Typically, we wouldn't change tau
significantly during a study because the time step error would change
but we are doing it for the tutorial to keep the number of time slices
small for different temperatures). Now, run your simulation and again
view the simulation using pathvis++. You should find paths much like
this:\

**^4^He at 4.2 K**
![Image:Helium\\_4.2.png](Helium\_4.2.png "Image:Helium\_4.2.png")

Notice that these paths are significantly larger than the paths at 64 K.
Their non-negligible extent is a manifestation of zero point motion. In
this regime, the system would be poorly represented by a classical
potential (which "point" in space along the path would you even use to
represent the classical particle?). There are many simulations in which
zero point motion is important (for example ambient water), even though
exchange statistics aren't critical. You can tell that statistics aren't
relevant in this simulation because the paths aren't permuting onto each
other. The inter-particle separation is sufficiently large that a
permutation would be very unfavorable in terms of the energetics of the
springs. Notice that in the simulations at both 6 K and 64 K, the
winding number is 0.

Now, let's run a simulation where the bosonic nature of Helium is
important. Set the temperature to 2.0K (let tau be 0.025). Run this
simulation and view the results using pathvis++. As you push the slider
and advance through Monte Carlo steps we see the paths in our system
expand and begin to permute onto each other. The existence of permuted
paths indicates the importance of boson exchange statistics in this
simulation. At some point, you will see a series of permuting paths wind
across the entire box (to see this, it is sometimes helpful to toggle to
Wrap mode with the button at top. This forces the paths to stay within
the box whereas the other mode tries to enforce continuity of the
paths). This wrapping around the box is an indication that the system is
superfluid. See mathematically how [superfluidity manifests itself
within Path Integral Monte Carlo](Path%20Integral%20Superfluidity).
Winding paths are indicated in yellow. Of course, the jumping around of
these paths sometimes obfuscates their global structure which makes it
hard to see what's going on. If we click "Smooth", these paths will be
smoothed by truncating their Fourier spectrum. The "Detail" slider
controls how many Fourier components are retained, i.e. you can adjust
this slider to tune the desired amount of detail in the path. Finally,
once you have a configuration and perspective you like, you can click on
tubes to add depth cues to the paths. At this point, you should have a
visualization similar to\

**^4^He at 1 K**
![Image:Helium\\_1K.png](Helium\_1K.png "Image:Helium\_1K.png")

In summary, we see that as the temperature decreases, quantum effects
become increasingly important, specifically zero point motion and then
quantum statistics.

**Summary**

[![64 K](Helium.64.png "64 K")]( "64 K")

64 K

[![4.2 K](Helium\_4.2.png "4.2 K")]( "4.2 K")

4.2 K

[![1 K](Helium\_1K.png "1 K")]( "1 K")

1 K

Before we move on, we should make sure we are using simulation
parameters (time step and system size) that lead to a well-converged
result (i.e. introduce negligible systematic errors). We will also want
to tune the algorithm (i.e., you may notice that you are spending too
much time computing observables and not enough time actually performing
the Monte Carlo moves of the simulation).

1.  [Checking
    Convergence](Checking%20Convergence%20(Liquid%20Helium%20Tutorial))
2.  [Algorithm
    Tuning](Algorithm%20Tuning%20(Liquid%20Helium%20Tutorial))

\
 Of course, as we change other parameters, we would want to re-optimize
these settings. However, an initial adjustment is adequate for our
purposes.

Superfluid Transition
------------------------------------------------------------------------------------------------------------

Let us now proceed to calculate the superfluid transition temperature of
liquid helium. We can start to get a rough estimate of this value by
measuring the superfluid fraction at a variety of different
temperatures. Let us do a series of calculations from \\(1.0 - 3.0 K\\). The
temperature is set as \\(1/\beta\\) where \\(\beta\\) is set by the product of
the input variables tau (Section: System) and NumTimeSlices (Section:
System). Let's make sure our time step \\(\tau\\) is 0.025.

There is a only a discrete set of temperatures (in our relevant
temperature range) we can choose because \\(\beta\\) must be an integral
multiple of tau, whose values are constrained by the density matrix
files that were generated. Set up 5 different simulations in this
temperature range (by changing the NumTimeSlices (Section: System)
variable and see if you can locate the Superfluid transition (remember,
you will also need to change the name for the OutFileBase (Section:
Observables) if you want it to write different output files). Once these
runs have completed, use the python analysis script Report.py to examine
them

    Report.py theFileName

and then plot the Superfluid Fraction as a function of temperature. (5
minute runs should be sufficient to get a rough outline of the graph
especially at higher temperature). A reasonably converged graph should
look like this (the red points are from the PIMC++ simulation. The green
lines are experimental results from Donnelly, 1967):\

![Image:SuperfluidGraph\\_10ptcl\\_longerRun.png](SuperfluidGraph\_10ptcl\_longerRun.png "Image:SuperfluidGraph\_10ptcl\_longerRun.png")\
 An example of what you may see after 5 minutes is
[here](SuperfluidGraph\_10ptcl)

Although we know that this transition must be sharp in the thermodynamic
limit (\\(N \rightarrow \infty, N/V\\) constant), it appears very smooth and
not at any definitive temperature in our graph. At best we can only
qualitatively locate the transition temperature. This is because 10
particles is very far from the thermodynamic limit. (Remember that phase
transitions are only formally defined in the thermodynamic limit) We
could obviously mitigate these finite size effects by going to a system
of more than 10 particles (and we will do so in a moment). (An even
better method is to do a variety of system sizes and then use a finite
size scaling analysis to pinpoint the exact transition temperature, but
an explanation of how this works is beyond the scope of the tutorial.)
Before we embark on these larger systems, though, it behooves us
understand how to use parallelization with PIMC++.

Because we want to calculate the superfluid fraction as quickly as
possible, we will want to take advantage of different modes of
parallelization within PIMC++. The code supports two modes of
parallelism: cloning and time slice parallization.

### Cloning Parallelization

**We will skip this section of the tutorial, since the computer lab is
not set up with a parallel run environment.** You can skip to the next
[Section](Bulk%20helium#Superfluid%20Transition%20Temperature%20Calculation%20(larger%20system)).
You can read on if you are interested in parallelization in PIMC++.
However, the lab is not set up to run PIMC++ in parallel.

The first mode of parallelism is to run many independent copies of the
simulation. Because the error bars of our simulation depend on the
number of independent samples we collect, Monte Carlo is "trivially
parallel" in this sense. PIMC++ automatically takes advantage of this
form of parallelism as its *default* behavior. If you run the code
(using mpirun or the equivalent for your MPI implementation) on more
than one processor, it will run independently on each processor,
creating output files labeled basename.PROCNUMBER.h5. Let's see this now
by running PIMC++ on multiple processors. For a preconfigured system
(these machines aren't), you can run\

    runPIMC Helium.in -np 4

\

Once the run is complete, you should see Helium.0.h5, Helium.1.h5,
Helium.2.h5 and Helium.3.h5. Running the python analysis script

     Report.py Helium 

will analyze all 4 of these files simultaneously. You should see that
the simulation has error bars that are approximately half as large as
they were earlier.

### Time slice Parallelization

**We will also skip this section, since the computer lab is not set up
with a parallel environment.**

In discussing time slice parallelization, we will refer to two important
quantities. \\(NP\\) refers to the number of processors allocated to PIMC++
at runtime. For example, if we had used the generic MPI command %mpirun
-np 4 pimc++ myInput.in or the Tungsten script %runPIMC nyInput.in -np
4, NP would be 4. The other important quantity is ProcsPerClone,
specified in the PIMC++ input file. For the Cloning mode of
parallelization discussed in the previous section, ProcsPerClone=1.

Unfortunately, for very large systems, the trivial (or natural!)
parallelization of Monte Carlo may not be enough. This is because the
amount of wall clock time it takes for the system to come into
equilibrium may become very large. In this case, we can also parallelize
over the time slices. In this case, the path of M slices is divided
between n processors (i.e. processor 0 can be in charge of slices 0-20,
processor 1 can be in charge of slices 20-40, processor 2 in charge of
40-60, etc.). The two modes of parallelism can be combined.

Effective use of time slice parallelism can make one's wallclock running
time largely independent of temperature. We can utilize this mode of
parallelism by editing the Section:Parallel of our input file. The
relevant variable, ProcsPerClone, sets how many processors are used for
each independent copy of the simulation. The time slices will then be
distributed evenly over these processors. To choose the correct number
of processors, we should keep the following considerations in mind. To
begin with, it is a rule that each processor must have at least \\(2^l\\)
slices where l is the maximum level that is used in any Bisection Move
throughout the simulation. On the other hand, you don't want to have too
many slices per processor or we won't gain as much a speedup as desired.
Also, the total number of processors must be a multiple of
ProcsPerClone. If \\(NP \neq   ProcsPerClone\\), PIMC++ will run
\\(\frac{NP}{ProcsPerClone}\\) clones, each of which has ProcsPerClone
processors assigned to it, thereby using both methods of parallelism in
the same simulation. When \\(NP=ProcsPerClone\\), a run will produce only
one output file that looks as if the simulation had been run in series.

Let's now experiment with these modes of parallelism. We will do a
system of 10 time slices split over 2 processors so each processor has 5
slices (practically you would more likely do 500 slices spread over 25
processors so each processor had 20 slices). Set *ProcsPerClone*
(Section:Parallel) equal to 2 instead of 1 and the *NumTimeSlices*
(Section: System) variable to 10. Now, run your job on 2 processors. On
Tungsten,

    runPIMC Helium.in -np 2

You'll notice that it produces a single output file. Run the analysis on
this output file and examine it. If you had run this same temperature
earlier, you should notice your error bars are approximately a factor of
\\(\sqrt{2}\\) smaller than they were previously.

Superfluid Transition Temperature Calculation (larger system)
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

We will now try to pin down the superfluid transition temperature more
accurately by doing a larger system. In the winter school, we will do
this in a group setting (if you are following this tutorial yourself,
acquiring all the data will involve performing several longer runs). We
will do a system of 32 particles and again try to find the transition
temperature. With 32 particles, the size of the box is now
\\(11.361397~\AA\\) for the same saturated vapor pressure density. Everyone
should make sure they make the changes:

-   *NumParticles* (Section:System/Particles/Species) to 32 and
-   *Box* (Section:System) to [11.361397,]

Now, choose a temperature in the range from 1-3K (preferably something
different from your neighbor) and adjust *NumTimeSlices* (Section:
System) to reflect that.

Now, everyone should run their simulation using

    pimc++ Helium.in

As it is progressing, run your analysis script (Report.py) to make sure
it looks like it is producing reasonable data. We will collect these
numbers on the board and try to better locate the transition
temperature.

Continue now to [H2
molecule](http://cms.mcc.uiuc.edu/pimcpp/index.php/H2)
