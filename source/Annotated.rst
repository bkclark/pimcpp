::

    // ANNOTATED TUTORIAL INPUT FILE
    // This file is still valid input; the parser will ignore all comments

    double tau=.05; // 1/kT -- sets the time step and units

    // Specify the parallelization strategy
    Section (Parallel)
    {
      int ProcsPerClone = 1;
      // For ProcsPerClone > 1, each clone will use 
      // ProcsPerClone processors for time-slice parallelization.
      // Otherwise, each processor will run an independent simulation with a different RNG seed.
    }

    // Description of physical system to be simulated
    Section (System)
    {
      int NumTimeSlices=12; // Inverse temperature beta = tau * NumTimeSlices
      Array<double,1> Box(3)=[10.0,10.0,10.0];  // dimensions of the simulation box -- make sure units are consistent
      Array<bool, 1 > IsPeriodic(3)=[true,true,true];  // Periodic boundary conditions for each dimension
      Section (Particles)
      {
        Section (Species)
        {
          string Name="free";  // arbitrary, but must be unique and used consistently
          string Type="free";  // H, O, Al, etc.
          double lambda=1.0;   // hbar^2/(2*m) in appropriate units
          string Statistics="BOLTZMANNON"; // other options are "BOSON" or "FERMION"
          int NumParticles=1;
          int NumDim=3; // dimensionality
          string InitPaths="BCC"; // See documentation on initialization options
        }
      } 
    }

    // Defines how to compute the action (tau * Energy) difference
    // used to accept or reject moves
    Section (Action)
    {
      int NumImages=0;
      int MaxLevels = 2;
      Array<string,1> PairActionFiles(1) = ["zero.PairAction"];
      // These files define the action as a function of distance
      // between two particles of the given species.
      // Files must be provided for all relevant combinations.
      // The file zero.PairAction defines the (null) free-particle action
    }

    // Set up measurement of quantities of interest
    Section (Observables)
    {

      string OutFileBase = "SingleParticle"; // base name for output file

      // Basic energy observable  
      Section (Observable)
      {
        string Type = "Energy"; // specifies energy observable
        string Name = "Energy"; // arbitrary but must be unique and consistent
        string Description="Total Energy";
        int Frequency=1; // Accumulates every Frequency times it is called
       }
    }   

    // Define moves to update system configurations
    Section (Moves){

      // See documentation for more information about this method
      Section (Move) {
        string Type="BisectionBlock";
        string Name="BisectionBlock";
        string PermuteType="NONE";
        string Species="free";
        int NumLevels=2;  
        int StepsPerBlock=2;
      }

      // This is a pseudo-move that is needed for technical reasons
      // when multiple time slices are used.
      // See documentation for more information.
      Section (Move)
      {
        string Type="ShiftMove";
        string Name="Shift";
      }
    }  

    // Define the sequence and frequency of events
    // during the simulation
    Section (Algorithm)
    {
      Section (Loop)
      {
        int Steps=1000; // Do these events 1000 times
        Section (Loop)
        {
          int Steps=50; // Do these events 50 times
          Section (Move) {string Name="BisectionBlock";} // Make move; accept or reject
          Section (Observe) {string Name = "Energy"; }   // Accumulate energy of current configuration
          Section (Move) {string Name = "Shift"; }       // do Shift pseudo-move
        }
        Section (WriteData){} // Dump accumulated observable data to output file
      }
    }

