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
           string Name="He";
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
       Array<string,1> PairActionFiles(1) = ["He4_cut.PairAction"];

     }

     Section (Observables)
     {
       string OutFileBase = "Helium_05";
       Section (Observable)
         {
           string Type = "Energy";
           string Name = "Energy";
            string Description="Total Energy";
           int Frequency=10;
         }

    Section (Observable)
         {
           string Type="StructureFactor";
           string Name="HeliumStructureFactor";
           double kCutoff=5.0;
           int Frequency=25;
           string Species1="He";
           string Species2="He";

         }
        Section (Observable)
          {
            string Type="SuperfluidFraction";
            string Name="SuperfluidFraction";
            Array<string,1> SpeciesList(1)=["He"];
            int Frequency=1; 
          }

       Section (Observable)
         {   
           string Type="PathDump";
           string Name="PathDump";
           int Frequency=100;
         }
       
       Section (Observable)
         {
           string Type = "PairCorrelation";
           string Name = "HeHePC";
           string Species1 = "He";
           string Species2 = "He";
           string Description="Helium-Helium Pair Correlation";
           int Frequency=2;
           Section (Grid)
             {
               string Type = "Linear";
               double start = 0.0;
               int NumPoints = 100;
             }
         }
       
       Section (Observable)
         {
           string Type="TimeAnalysis";
           string Name="TimeAnalysis";
           int Frequency=1;
         }
     }   


     Section (Moves){
       Section (Move) {
         string Type="BisectionBlock";
         string Name="BisectionBlock";
         string PermuteType="TABLE";
         string Species="He";
          int NumLevels=2;
         int StepsPerBlock=10;
         Array<double,1> Gamma(4)=[1.0,10.0,10.0,10.0];
         double epsilon=1e-5;
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
           int Steps=5000;
           Section (Move) {string Name="BisectionBlock";}
           Section (Observe) {string Name = "Energy"; }
    //       Section (Observe) {string Name = "HeliumStructureFactor"; }
           Section (Observe) {string Name = "PathDump"; }
    //       Section (Observe) {string Name = "HeHePC"; }
           Section (Observe) {string Name = "TimeAnalysis"; }
           Section (Observe) {string Name = "SuperfluidFraction"; }
           Section (Move) {string Name = "Shift"; }
         }
         Section (WriteData){}
       }
     }

