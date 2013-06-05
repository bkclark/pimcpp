::

    double tau=.05;
    Section (Parallel)
    {
      int ProcsPerClone = 1;
    }
    Section (System)
    {
      int NumTimeSlices=12;
      Array<double,1> Box(3)=[7.709894,7.709894,7.709894];
      Array<bool, 1 > IsPeriodic(3)=[true,true,true];
      Section (Particles)
        {
          Section (Species)
        {
          string Name="Helium4";
          string Type="He";
          double lambda=6.059615;
          string Statistics="BOSON";
          int NumParticles=10;
           int NumDim=3;
           string InitPaths="BCC";
         }
       } 
      }
     Section (Action)
     {
       int NumImages=0;
       int MaxLevels = 2;
       Array<string,1> PairActionFiles(1) = ["He4.95.OffDiag.PairAction"];

     }

     Section (Observables)
     {
       string OutFileBase = "Helium";
       Section (Observable)
         {
           string Type = "Energy";
           string Name = "Energy";
            string Description="Total Energy";
           int Frequency=1;
         }
     }   


     Section (Moves){
       Section (Move) {
         string Type="BisectionBlock";
         string Name="BisectionBlock";
         string PermuteType="TABLE";
         Array<double,1> Gamma(4)=[1.0,10.0,100.0,400.0];
         double epsilon=1e-5;
         string Species="Helium4";
         int NumLevels=2;
         int StepsPerBlock=2;
       }
       Section (Move)
         {
           string Type="ShiftMove";
           string Name="Shift";
         }
      }  

     Section (Algorithm)
     {
       Section (Loop){
         int Steps=1000;
         Section (Loop){
           int Steps=50;
           Section (Move) {string Name="BisectionBlock";}
           Section (Observe) {string Name = "Energy"; }
           Section (Move) {string Name = "Shift"; }
         }
         Section (WriteData){}
       }
     }

