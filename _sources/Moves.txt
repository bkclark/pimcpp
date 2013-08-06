Moves
=====

To add a move to your calculation, there are two steps. One needs to
define the move and then one needs to tell PIMC++ when to use the move
in the algorithm.

(Note: Moves sometimes write data like acceptance ratios. This data is
written when the algorithm reaches the WriteData section)

Defining the move
-----------------

All moves are defined in their respective section **Move** in the
overall section **Moves**

::

 Section(Moves)
 {
   Section (Move)
   {
     string Type="BisectionBlock";
     string Name="HeliumBisection";
     string PermuteType="NONE";
   }
   Section (Move)
   {
     ...
   }
 }

All moves have (at minimum) the variables *Type* and *Name* must be
specified. The variable *Type* tells the code which move to use (i.e.
bisection block, displace move, etc). The string to use may be found
below in the list of moves. The variable *Name* is a name you choose for
the move to be used later on when you call it in the algorithm. For
example, your Name might be "BisectionMoveForHelium".

Inserting the move into the algorithm
-------------------------------------

A sample algorithm looks like:

::

 Section (Algorithm)
 {
   Section (Loop){ //Start Accumulating observables
     int Steps=3000000000;
     Section (Loop){
       int Steps=100;
       Section (Move) {string Name="BisectionMoveForHelium";}
       Section (Loop){
         int Steps=10;
         Section (Move) {string Name="Displace";}
       }
       Section (Move) {string Name="CenterOfMass";}
       Section (Observe) {string Name = "Vacancy"; }
       Section (Observe) {string Name = "PathDump"; }
       Section (Observe) {string Name = "StructureFactor"; }
       Section (Observe) {string Name = "TimeAnalysis"; }
       Section (Observe) {string Name = "CycleCount"; }
       Section (Move) {string Name = "Shift"; }
     }
     Section (WriteData){}
   }
 } 

To add your move to the algorithm, you must place in a line like:
Section (Move) {string Name="BisectionMoveForHelium";} where the name is
the name you have chosen for your move.

List of Moves
-------------

Note: All Moves include a paramater *Name* that is user specified.

Bisection Block
^^^^^^^^^^^^^^^

The bisection block is the primary move in PIMC++. It implements the
Bisection move as documented in David Ceperley rev mod phys article
(somebody please link to this).

The bisection block has a number of sophisticated switches that turn
certain actions on and off. I am currently documenting it's most typical
use.

*Example input:*

::

 Section (Move) {
     string Type="BisectionBlock";
     string Name="BisectionBlock";
     string PermuteType="TABLE";
     string Species="He";
     Array<double,1> Gamma(4) = [1.0, 1.0, 1.0, 1.0];
     double epsilon=1e-5;
     int NumLevels=2;
     int StepsPerBlock=100;
   }

*Input Paramaters:*

-  **Type:** BisectionBlock
-  **Species:** Specifies a string for the Name (NOT Type) of the
   species BisectionBlock acts on (this species must show up as the Name
   of a species specified earlier in the input file).
-  **NumLevels:** Number of levels in which to perform a bisection.
-  **StepsPerBlock:** A BisectionBlock chooses a set of time slices (of
   size 2^level) in which to work. Then it makes moves on these time
   slices StepsPerBlock times. This is done to amortize the cost of
   building tables for the permutations. If you are using permutations
   in the code it is suggested that you select StepsPerBlock to be
   between the number of particles and number of particles squared.
-  **PermuteType:** Specifies how permutations are done on the species.
   The current choices are:

   -  NONE: No permutations are done
   -  TABLE: Permutations up to size four are done. These permutations
      are chosen by building a table of all possible four particle
      permutations and choosing a permutation proportional to (check
      that) where d is the distance between the first and last time
      slice the Bisection Block is working on. If TABLE is chosen, there
      are two additional input paramaters that are to be used:

      -  **epsilon:** A double that sets a minimum unnormalized
         probability for any permutations for permutations which are
         pushed onto the table to be allowed to have. Permutations with
         an unnormalized probability less then this are not included.
         Because of how things are implemented this leads to a small
         bias.
      -  **Array Gamma(4):** An array of size four that specifies a
         multiplicative factor to increase the unnormalized probability
         of permutations of size 1,2,3, and 4 respectively. These
         multiplicative factors must all be above 1 (we should check
         this in the code!)

*Output Paramaters*

-  **AcceptRatio:** Contains a global acceptance ratio of the entire
   system as well as an acceptance ratio of each bisection stage. The
   remaining stages are different levels of the bisection move where the
   last stage is the lowest level. Each stage specifies the percentage
   of moves that got to that level and then succeeded in getting to the
   next stage. (So the acceptance ratio of the total move should be the
   product of all the stages).
-  **CenterOfMassDrift:** This specifies how far the squared center of
   mass has been moved by this move. By dividing it by the time spent in
   this move you can calculate your average diffusion rate.

Center of Mass Move
^^^^^^^^^^^^^^^^^^^

When PIMC++ is initialized, it calculates the center of mass of all the
particles. Then after a move has been accepted it calculates how that
center of mass is changed. This is done to resolve problems concerning
calculating the center of mass in a periodic box. This is done
independent of whether or not this move is defined.

If this move is defined and called in the algorithm, the entire system
is moved so that the center of mass displacement is now 0. The move
always accepts (if you are running a simulation that has an action (such
as an external potential) that changes when the entire system is moved
by a constant DO NOT USE THIS MOVE as it is currently implemented)

*Example input:*

::

 Section (Move){
     string Type="CenterOfMass";
     string Name="CenterOfMass";
   }

*Input Paramaters:*

- **Type:** CenterOfMass

*Output Paramaters:*

- **AcceptRatio:** Although appearing in the output file, this is irrelevant because everything is always accepted.

Displace Move
^^^^^^^^^^^^^

The displace move chooses a particle among all particles of a certain
species and displaces that particle a distance d where d is chosen from
a gaussian with width sigma. The displace move does not attempt to move
a particle if it is part of a permutation. Note: There is a bug in the
displace move that well cause the code to loop forever if all particles
are in a permutation. This needs to be fixed. Note: The displace move
doesn't push the kinetic action back onto its action list and
consequently will never calculate it. If for some unknown reason your
system changes the kinetic action when particles are displaced, this
will obviously cause a bug.

*Example input:*

::

 Section (Move) {
     string Type="Displace";
     string Name="Displace";
     double Sigma=0.5;
     Array<string,1> ActiveSpecies(1)=["He"];
     int NumToMove=1;
   }

*Input Paramaters:* 

- **Type:** Displace
- **Sigma:** Width of gaussian from which displacement is chosen
- **Array ActiveSpecies(1):** List of species which displacement acts upon.
- **int NumToMove:** Number of particles in which to move at one time.

*Output Paramaters:*

- **AcceptRatio:** Outputs the percent of these moves that have been accepted. Due to technical details of the implementation also outputs the acceptance ratio of a single stage that is the same value.

Shift Move
^^^^^^^^^^

All code in PIMC++ follows the invariant that the initial and final time
slice in memory are not to be altered. Notice, without time slice
parallelization these two time slices will be permutations of one
another. With time slice parallelization, the last time slice of
processor i will be a permutation of the first time slice of processor
i+1 (modulo the number of processors). Keeping this invariant is crucial
to the correct functioning of the code.

A Shift Move shifts the data by a number of time slices (randomly chosen
in some range) to allow time slices that are on the border between
processors (or the last time slice on a processor if there is only one
processor) to be moved.

It is very important to use a shift move occassionally or certain time
slices will never be updated.

Shift Moves are always accepted.

*Example input:*

::

 Section (Move)
     {
       string Type="ShiftMove";
       string Name="Shift";
     }

*Input Paramaters:*

- **Type:** ShiftMove

*Output Paramaters:*

- **AcceptRatio:** Although appearing in the output file, this is irrelevant because everything is always accepted
