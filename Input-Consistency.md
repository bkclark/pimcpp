---
layout: base
title: Input-Consistency
---

There are certain paramaters in the input that must be consistent for
the code to work. Here is a list of some of thos necessary input
consistencies:

If you choose Boson as your species type then your BisectionBlock needs
to have Permutations turned on (currently at the only chose of TABLE).

The total number of levels read in the action must be greater then the
number of levels used in the BisectionBlock.

The total number of time slices must be more then 2\^Levels used in the
bisection block.

If OpenLoops is on, you must be using an approximate action at the
higher level in your bisection block.
