---
layout: base
title: System
---

There are a number of different paramaters used to define the physical
system on which the simulation is performing. They include things such
as the number of time slices, the time step, information about the
species, the box, etc.

They are all included in the section **System**

*Example input:*

       }
     
     }

At the level of **System** (i.e. not at a lower level in the hierarchy)
the parameters are:

**NumTimeSlices**: Specifies the number of time slices for the system.
The beta of the system is defined as tau\*NumTimeSlices

**tau**: Specifies the time step (in the trotter breakup) for the
simulation.

**Array\<bool,1\>\</bool,1\> IsPeriodic(3)**: Specifies which dimensions
the system is periodic in. Currently an entirely periodic box is well
tested and a box with no periodicity should largely work. A box that is
periodic in only some dimensions is not well implemented right this
second but could be without too much difficulty.

**Array\<double,1\>\</double,1\> Box(3)**: Specifies the size of the box
for the system. Must be the number of dimensions of the entire system.

The **System** section, contains a section **Particles** that defines
parameters associated with the particles.

The **Particles** section defines a **Species** section for each species
in the simulation.

Each **Species** in the simulation is set up in a section **Species**
inside the section **Particles *which is inside the section*** *System*

Example input:

     Section (Species)
        {
          string Name="He";
          string Type="He";
          double lambda=6.059615;
          string Statistics="BOSON";
          int NumParticles=48;
          int NumDim=3;
          string InitPaths="BCC";
        }

The paramaters that are required to define a species are:

*Paramaters:* **Type:** Specifies the type of the species. This is used
to look up its Pair Action and must correspond to the corresponding name
inside the PairAction file.

**Name:** This is a user-supplied name that will be used to reference
the species. lambda: Equal to \\hbar\^2/2m.

**Statistics:** string variable indicating the statistics of the
particle. Choice of "Boson", "Fermion" or "Boltzmannon". Currently using
the "Fermion" flag turns on some of the nodal actions. None of the other
flags currently do anything. Maybe we should check if the Boson flag is
on before allowing permutations to happen.

**NumParticles:** How many particles of this type are being used in the
simulation.

**NumDim:** The dimensionality of this particle. Lower dimensions then
the number of dimensions of the simulation will allow things such as
strings and planes as particles but this currently does not work.

**InitPaths:** string used to specify how to initialize the particles.
Choices are:

1.  **BCC:** Sets up the particles in a BCC lattice
2.  **FIXED:** Allows you to specifiy the intial conditions of all the
    particles. All time slices will start at the same initial condition.
    Using this paramater requires the use of the paramater
3.  **Array \<double,2\>\</double,2\> Positions(numParticles,NDIM):**
    Specifies the particle locations.
4.  **FILE:** Allows you to specify a file from which the last
    configuration in the PathDump of the file will be read. Using this
    flag requires the user of the paramater
    1.  **string File="fileName":** Specifies the BASE (no .0.h5) of the
        filename to be read. If cloning is being used it will attempt to
        read fileName.n.h5 for the n'th clone.


