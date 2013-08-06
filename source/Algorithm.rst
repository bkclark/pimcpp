Algorithm
=========

The algorithm in the input specifies the order in which
moves/observables/etc. are called as well as when the data is written.

A typical example of an algorithm is

::

 Section (Algorithm)
 {
   Section (Loop) //Start Accumulating observables
   {
     int Steps=100;
     Section (Loop){
       int Steps=10000;
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

One basic aspect of the algorithm is the Loop. The only important
paramater is

**Steps**: Number of steps that the loop should execute.

The moves and observables are defined by including a Section Move and
Observe respectively that includes a single parameter **Name** that
defines the name as used in the definition of the Move and Observable
earlier in the input file.

To write out the data, a Section WriteData is included.
