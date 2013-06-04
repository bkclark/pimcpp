---
layout: base
title: Parallel
---

The parallel section of the of the input section is used to specify
paramaters associated with how the system works in parallel.

A typical parallel input section looks like Example input:

     Section (Parallel)
     {
       int ProcsPerClone=1;
     }

*Input Paramaters*:

**ProcsPerClone**: This sets the number of processors over which the
time slices will be parallelized (i.e if ProcsPerClone is 4 and there
are 200 time slices then each processor will have 50 time slice). Your
ProcsPerClone must divide your total number of processors. If you have
more then ProcsPerClone processors, the rest of the processors will be
used to do cloning (i.e. identical runs with different random numbers).
These will be stored as basename.procnum.h5
