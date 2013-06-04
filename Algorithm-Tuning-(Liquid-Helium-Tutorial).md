---
layout: base
title: Algorithm-Tuning-(Liquid-Helium-Tutorial)
---

The speed of a Path Integral simulation is dependent on a number of
factors. We would like to maximize the rate of diffusion of the system
through configuration space (or equivalently minimize the number of
seconds it takes to get a statistically independent sample from our
Markov Chain). PIMC++ has a number of features built in to help you
measure the efficiency of your simulation. Moreover, beyond a
qualitatively different choice of move algorithms, there are a number of
quantitative knobs that we can tune to speed up the simulation.

Ways to monitor the speed of your simulation (open up your HTML analysis
files and find where these are in the output):

-   Time analysis observable (must be turned on in your input file):
    This tells you what percent of time you spend on each move and
    observable in your system. If you find that you're spending too much
    time in an observable, turn down the *Frequency*
    (Section:Observables/Observable) at which it is being observed. For
    moves, a (provably?) good heuristic is to spend the same percentage
    of time in each move (bounds your regret against the "perfect"
    algorithm). Also, this gives you an idea of what moves you might
    optimize.
-   Moves Acceptance ratios: Tells you the percentage that your move is
    being accepted. Low acceptance can be indicative of a move being not
    useful in your simulation (a good heurestic for low is less then
    10%). Also, it can indicate that you need to do fewer levels in a
    multi-stage move. Information on the individual stages and
    individual permutations are also supplied to allow fine-tuning of
    your *Gamma* (Section: Moves/Move (Bisection Block)) parameters and
    number of levels.
-   Move center of mass diffusion: Indicates how far each move has
    helped the center of mass of the path diffuse. Sometimes a better
    indicator of the efficacy of a move. A low acceptance ratio may
    still be acceptable if the accepted moves push around the path
    significantly. Less useful in situations where the path center of
    mass is poorly defined (think winding paths). To use this metric as
    a fair comparison between different moves, you should divide this by
    the number by time spent in the BisectionBlock move (calculable in
    the Time Analysis observable)

There are a number of useful knobs to turn that have the potential to
make a significant difference in speed. They include:

-   Bisection Level: The higher this is, the more the path moves at once
    but the lower your acceptance ratio
-   Frequency: The less often you observe, the faster it is (but the
    worse your statistics may be)
-   StepsPerBlock (Section: Moves/Move (BisectionBlock): In the
    bisectionblock, it does StepsPerBlock moves at a single set of time
    slices before continuing in the algorithm. For Boltzmannons there is
    no efficiency improvement by having a large number here (in fact, if
    it's much more then the number of particles it might be marginally
    inefficient since you'll be ignoring other slices. For Bosons, there
    is a large upfront cost that is amortized over. Consequently higher
    numbers here help efficiency significantly. Heuristically, this
    should be between \\(N\\) and \\(N^2\\).

