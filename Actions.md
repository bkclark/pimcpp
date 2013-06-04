---
layout: base
title: Actions
---

PIMC++ supports a number of Actions that can be turned on and off
throughout the simulation.

PIMC++ samples the stationary distribution

where S is the action.

     A typical simulation of PIMC++ uses a Kinetic Action, a Potential Action (usually a pair action of some sort) and possibly a nodal action.

Other actions may be used through the code that are not being sampled in
the stationary distribution but are being used in some other manner. For
example, an action might be used at a higher level in a multistage move
to approximate the higher levels.

A typical run of PIMC++ requires certain paramaters set. These
paramaters are set in the Actions section (At some point I think we
would like to have hiearchal sections for kinetic action, potential
action, etc.)

Kinetic Action Paramaters: NumImages: The kinetic action formally is

where k,k' index's the k,k' image of particle i,j. This paramater
defines whether the sum is only over images in the first box (=0), the
nearest 27 boxes (=1), etc. Typically this can be 0 unless

is on the order of the box size. Potential Action Paramaters:
Array\<string,1\>\</string,1\> PairActionFiles(numPairActions):
Specifies the PairActionFiles. The full path must be given to these
files. The PairAction files that can be used in PIMC++ include those
from squarer++ and those from David Ceperely's squarer. To use the
latter requires writing a PairAction file that wraps the dm files from
this squarer. See below for a description to accomplish this. MaxLevels:
Specifies the total number of levels that are read in from the density
matrix file. This must be less then (or equal to) the total number of
levels in teh density matrix. It also must be more then the number of
levels that you are using in bisection (which is why more then onle
level is read in). The highest level that is read in is the value of tau
tbeing used. Wrapping David Ceperley's squarer dm files in a PairAction
files such that they can be used.

(details to come)
