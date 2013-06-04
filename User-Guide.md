---
layout: entry
title: User-Guide
---

The general workflow of Path Integral Monte Carlo involves the
following:

You propose a trial move, check to see if it is accepted or rejected
based on its respective action, and repeat. Ocassionally you accumulate
observables, These steps all operate on a set of paths which consist of
a number of particles of a variety of species.

PIMC++ is highly configurable allowing you to select which set of moves
will be made on your particles at what time, which actions are to be
used for these moves, and then when and which observables you will
accumulate. It also allows you to configure properties of your
particles.

In order to start using PIMC++, you might take a look at the following
choices:

-   InputExamples: Generic input file I can use to get started
-   Input File Explanation
    -   <Observables>: Adding an observable to the system and a list of
        observables, their descriptions and parameters.
    -   <Moves>: Adding a move to the system and a list of moves, their
        description, and parameters.
    -   <Actions>: Actions that exist and how you set paramaters for
        them
    -   <System>: Changing the nature of the system (adding species,
        etc.)
    -   <Algorithm>: How to change the order in which moves/observables
        and writing of data is called
    -   <Parallel>
    -   [Input Consistency](Input%20Consistency)

-   Analyzing the Output using the Report generator
-   Visualizing the paths of your system with PathVis

