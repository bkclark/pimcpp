/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "PIMCClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"
#include <sstream>
#include <fstream>
#include "Blitz.h"
#include "IO/FileExpand.h"
#include "QMCWrapper.h"


bool PIMCClass::Read(IOSectionClass &in)
{
  // tells whether to run or be a dummy
  bool doPIMCRun = false;

  // Read the parallelization strategy
  PathData.Read (in);

  // this is set to true in PathDataClass::Read when not built with qmcpack
  if(PathData.IAmQMCManager){
    doPIMCRun = true;
    // Read in the system information and allocate the path
    assert(in.OpenSection("System"));
    PathData.Path.Read(in);
    in.CloseSection();

#ifdef USE_QMC
    PathData.AssignPtclSetStrings();
#endif

    // if (PathData.Path.ExistsCoupling){
    //   int myProc=PathData.InterComm.MyProc();
    //   PathData.Path.ExistsCoupling=(double)(myProc)/100;
    // }

    // Read in the action information
    cerr <<PathData.Path.Communicator.MyProc()<<" Reading Actions"<<endl;
    assert(in.OpenSection("Action"));
    PathData.Actions.Read(in);
    in.CloseSection();

    // Now actually initialize the paths
    cerr <<PathData.Path.Communicator.MyProc()<<" Initializing Paths"<<endl;
    assert(in.OpenSection("System"));
    PathData.Path.InitPaths(in);



    in.CloseSection();

    if (PathData.Path.UseCorrelatedSampling())
      PathData.Path.SetIonConfig(0);

    // Init Actions caches
    cerr <<PathData.Path.Communicator.MyProc()<<" Initializing Actions Caches"<<endl;
    PathData.Actions.Init();

    // Read in the Observables
    cerr <<PathData.Path.Communicator.MyProc()<< " Reading Observables"<<endl;
    assert(in.OpenSection("Observables"));
    ReadObservables(in);
    in.CloseSection();

    // Check for root processor
    bool iAmRootProc = (PathData.Path.Communicator.MyProc()==0);

    // Create Actions section in output file
    if (iAmRootProc)
      OutFile.NewSection("Actions");

    // Append Long Range Action
    if (PathData.Actions.HaveLongRange()) {
      cerr << "Initializing Long Range" << endl;
      assert (in.OpenSection ("Action"));
      PathData.Actions.LongRange.Init (in, OutFile);
      if (PathData.Actions.UseRPA)
        PathData.Actions.LongRangeRPA.Init(in);
      in.CloseSection();
    }
    if (iAmRootProc) {
      PathData.Actions.WriteInfo(OutFile);
      OutFile.CloseSection(); // "Actions"
    }

    // Read in the Moves
    cerr <<PathData.Path.Communicator.MyProc()<<" Reading Moves"<<endl;
    assert(in.OpenSection("Moves"));
    ReadMoves(in);
    in.CloseSection();

    // Read in the Algorithm
    cerr <<PathData.Path.Communicator.MyProc()<<" Reading Algorithm"<<endl;
    assert(in.OpenSection("Algorithm"));
    ReadAlgorithm(in);
    in.CloseSection();
  } else {
    QMCWrapper = new QMCWrapperClass(PathData);
  }

  return doPIMCRun;

}


void PIMCClass::ReadObservables(IOSectionClass &in)
{
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot= myProc==0;
  if (iAmRoot) {
    string outFileBase;
    assert(in.ReadVar("OutFileBase",outFileBase));
    int fileStart;
    if (!in.ReadVar("FileStart", fileStart))
      fileStart = 0;
    // Allow for tilde-expansion in these files
    outFileBase = ExpandFileName (outFileBase);
    ostringstream cloneNum;
    cloneNum << (PathData.GetCloneNum() + fileStart);
    bool restart;
    if (!in.ReadVar("Restart",restart))
      restart=false;
    if (restart){
      stringstream tempStream;
      int counter=0;
      tempStream<<outFileBase<<"."<<counter<<"."<<(PathData.GetCloneNum()+fileStart)<<".h5";
      //cerr<<"Checking for "<<tempStream.str();
      while (fileExists(tempStream.str())){
	counter++;
	tempStream.str("");
	tempStream<<outFileBase<<"."<<counter<<"."<<(PathData.GetCloneNum()+fileStart)<<".h5";
	//cerr<<"Checking for "<<tempStream.str();
      }
      ostringstream counterNum;
      counterNum<<counter;
      OutFileName=outFileBase+"."+counterNum.str()+"."+cloneNum.str()+".h5";
    }
    else{
      OutFileName = 
	outFileBase+ "." + cloneNum.str() + ".h5";
    }
    OutFile.NewFile(OutFileName);
    /// This is needed so that all of the decendents of the root
    /// LoopClass object have a real output file that they can flush.
    /// In partcular, WriteData will flush the file.
    Algorithm.SetOutfile(OutFile);

    /////////////////////////////////////
    // Write input file to output file //
    /////////////////////////////////////
    ifstream infile;
    infile.open(in.GetFileName().c_str());
    infile.seekg(0,ios::end);
    int length=infile.tellg();
    infile.seekg(0,ios::beg);
    char *buffer=new char[length+1];
    infile.read(buffer,length);
    buffer[length] = '\0';
    infile.close(); 
    string fileCopy(buffer);   
    delete buffer;
    OutFile.WriteVar("InputFile",fileCopy);


    OutFile.NewSection("RunInfo");
    RunInfo.Write(OutFile);
    OutFile.CloseSection();
    OutFile.NewSection("System");
    WriteSystemInfo();
    OutFile.CloseSection(); // "System" 
    OutFile.NewSection ("Observables");
    Array<double,1> weights;
    if (in.ReadVar("Weights", weights)) {
      double myWeight = weights(PathData.GetCloneNum());
      OutFile.WriteVar("Weight", myWeight);
    }
  }
  int numOfObservables=in.CountSections("Observable");
  
  for (int counter=0;counter<numOfObservables;counter++){
    in.OpenSection("Observable",counter);
    string observeType, observeName;
    assert(in.ReadVar("Type",observeType));
    assert(in.ReadVar("Name",observeName));
    if (iAmRoot)
      OutFile.NewSection(observeType);
    ObservableClass* tempObs;
    if (observeType=="PairCorrelation") 
	tempObs = new PairCorrelationClass(PathData,OutFile);
    else if (observeType=="nofr")
      tempObs=new nofrClass(PathData,OutFile);
    else if (observeType=="PlaneDensity")
      tempObs=new PlaneDensityClass(PathData,OutFile);
    else if (observeType=="ParticleAverageLoc")
      tempObs= new ParticleAverageLocClass(PathData,OutFile);
    else if (observeType=="DropletSuperfluidity")
      tempObs = new SuperfluiDrop(PathData,OutFile);
    else if (observeType=="SpecificHeatA")
      tempObs = new SpecificHeatAClass(PathData,OutFile);
    else if (observeType=="SpecificHeat")
      tempObs = new SpecificHeatClass(PathData,OutFile);
    else if (observeType=="SuperfluidFractionPerLayer")
      tempObs = new SuperfluidFractionPerLayerClass(PathData,OutFile); 
   else if (observeType=="SuperfluidFraction")
      tempObs = new SuperfluidFractionClass(PathData,OutFile);
    else if (observeType=="Sign")
      tempObs = new SignClass(PathData,OutFile);
    //    else if (observeType=="Vacancy")
    //      tempObs = new VacancyLocClass(PathData,OutFile);
    //    else if (observeType=="VacancyNear")
    //      tempObs = new VacancyLocNearbyClass(PathData,OutFile);
    //    else if (observeType=="VacancyNear")
    //      tempObs = new VacancyLoc2Class(PathData,OutFile);
    //    else if (observeType=="VacancyDensity")
    //      tempObs = new VacancyDensityClass(PathData,OutFile);
    //    else if (observeType=="Conductivity")
    //      tempObs = new ConductivityClass(PathData,OutFile);
    else if (observeType=="Coupling")
      tempObs = new CouplingClass(PathData,OutFile);
    else if (observeType=="Energy")
      tempObs = new EnergyClass(PathData,OutFile);
    else if (observeType=="Hexatic")
      tempObs = new HexaticClass(PathData,OutFile);
    else if (observeType=="HBond")
      tempObs = new HbondClass(PathData,OutFile);
    else if (observeType=="MeanSqDiffusion" || observeType=="MSD")
      tempObs = new ObsDiffusionClass(PathData,OutFile);
    else if (observeType=="DistanceToOpen")
      tempObs = new HeadLocClass(PathData,OutFile);
    else if (observeType=="PhiK")
      tempObs = new PhiKClass(PathData,OutFile);
    else if (observeType=="Pressure")
      tempObs = new PressureClass(PathData,OutFile);
    //    else if (observeType=="VacancyLocation")
    //      tempObs = new VacancyLocClass(PathData,OutFile);
    else if (observeType=="TimeAnalysis")
      tempObs = new MCTimeClass(PathData,OutFile,Moves,Observables,
				PathData.Actions.ActionList);
    else if (observeType=="TimeLindenman")
      tempObs= new TimeLindenmanClass(PathData,OutFile);
    else if (observeType=="TimeHexatic")
      tempObs= new TimeHexaticClass(PathData,OutFile);
    else if (observeType=="Angular")
      tempObs = new AngularClass(PathData,OutFile);
    else if (observeType=="PathDump")
      tempObs = new PathDumpClass(PathData,OutFile);
    else if (observeType=="WindingNumber")
      tempObs = new WindingNumberClass(PathData,OutFile);
    else if (observeType=="Vacancy")
      //      tempObs = new VacancyLocClass(PathData,OutFile);
      //    else if (observeType=="CycleCount")
      tempObs = new PermutationCountClass(PathData,OutFile);
    else if (observeType=="StructureFactor")
      tempObs = new StructureFactorClass(PathData,OutFile);
    else if (observeType=="DynamicStructureFactor")
      tempObs = new DynamicStructureFactorClass(PathData,OutFile);
    else if (observeType=="Weight")
      tempObs = new WeightClass(PathData,OutFile);
    else if (observeType=="Forces")
      tempObs = new ForcesClass(PathData,OutFile);
    else if (observeType=="PairCorrelationReweighting")
      tempObs = new PairCorrelationReweightingClass(PathData,OutFile);
    //  else if ( "OpenOrientation")
    //   tempObs = new OpenOrientationClass(PathData,OutFile);
    else {
      perr << "We do not recognize the observable " << observeType << endl;
      abort();
    }
    tempObs->Name = observeName;
    tempObs->Read(in);
    Observables.push_back(tempObs);
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();//Observable
  }
  if (iAmRoot)
    OutFile.CloseSection(); // "Observables"
}




void PIMCClass::ReadMoves(IOSectionClass &in)
{

  int numOfMoves=in.CountSections("Move");
  int steps;
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot = (myProc == 0);
  MoveClass* move;  
  if (iAmRoot)
    OutFile.NewSection("Moves");
  for (int counter=0;counter<numOfMoves;counter++){
    in.OpenSection("Move",counter);
    string moveType, moveName;
    assert(in.ReadVar("Type",moveType));
    assert(in.ReadVar("Name",moveName));
    if (iAmRoot)
      OutFile.NewSection(moveType);
    if (moveType=="ShiftMove")
      move = new ShiftMoveClass(PathData, OutFile);
    else if (moveType=="PrintMove")
      move = new PrintMoveClass(PathData, OutFile);
    else if (moveType=="BisectionBlock")
      move = new BisectionBlockClass(PathData,OutFile);
    else if (moveType=="SwapMove")
      move = new SwapMoveClass(PathData,OutFile);
    else if (moveType=="CorrelatedBisectionBlock")
      move = new CorrelatedBisectionBlockClass(PathData,OutFile);
    else if (moveType=="CouplingMove")
      move = new CouplingMoveClass(PathData,OutFile);
    else if (moveType=="CenterOfMass")
      move = new CenterOfMassMoveClass(PathData,OutFile);
    else if (moveType=="ReadPath")
      move = new ReadPathClass(PathData,OutFile);
    else if (moveType=="BisectionSphereBlock")
      move = new BisectionSphereBlockClass(PathData,OutFile);
    else if (moveType=="CenterDroplet")
      move = new CenterDropletClass(PathData,OutFile);
    else if (moveType=="GrowWorm")
      move = new WormGrowMoveClass(PathData,OutFile);
//     else if (moveType=="CloseWorm")
//       move = new WormCloseMoveClass(PathData,OutFile);
//     else if (moveType=="RemoveWorm")
//       move = new WormRemoveMoveClass(PathData,OutFile);
    else if (moveType=="OpenEnd")
      move = new OpenEndMoveClass(PathData,OutFile);
    else if (moveType=="RefSlice")
      move = new RefSliceMoveClass(PathData,OutFile);
    else if (moveType=="Displace")
      move = new DisplaceMoveClass(PathData,OutFile);
    else if (moveType=="DisplaceFast")
      move = new DisplaceFastMoveClass(PathData,OutFile);
    //     else if (moveType=="WaterRotate")
//     else if (moveType=="VariationalDisplace")
//        move = new VariationalDisplaceMoveClass(PathData,OutFile);
		/// This is just here for debugging; should be deleted in the future
    else if (moveType=="WaterMove"){
			cerr << "ERROR: 'WaterMove' IS OBSOLETE: USE 'MoleculeMove' NOW" << endl;
			assert(0);
		}
    else if (moveType=="MoleculeMove")
      move = new MoleculeMoveStageManagerClass(PathData, OutFile);
    else if (moveType=="PreSampleMoleculeMove")
      move = new PreSamplingClass(PathData, OutFile);
    else if (moveType=="SPSMoleculeMove")
      move = new SPSClass(PathData, OutFile);

    else if (moveType=="LocalFlip")
      move =new LocalFlip(PathData,OutFile);

    else if (moveType=="LocalFlip")
      move = new LocalFlip(PathData, OutFile);
    else if (moveType=="GlobalFlip")
      move = new GlobalFlip(PathData, OutFile);

    else{
      perr<<"This type of move is not recognized: "<< moveType <<endl;
      abort();
    }
    move->Name = moveName;
    move->Read(in);
    Moves.push_back(move);
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();
  }
  if (iAmRoot) {
    OutFile.CloseSection (); // "Moves"
    OutFile.FlushFile();
  }
  //cerr<<"I have finished reading the moves"<<endl;
}



void PIMCClass::ReadAlgorithm(IOSectionClass &in)
{
  int maxWallTime;
  if (in.ReadVar("MaxWallTime", maxWallTime)) {
    PathData.SetMaxWallTime(maxWallTime);
    int hours = maxWallTime/3600;
    int minutes = (maxWallTime-3600*hours)/60;
    int seconds = maxWallTime%60;
    perr << "Maximum wall time is " << hours 
	 << ((hours != 1) ? " hours, " : " hour, ") << minutes
	 << ((minutes != 1) ? " minutes, and " : " minute, and ") << seconds 
	 << ((seconds != 1) ? " seconds.\n" : " second.\n");
  }
  //cerr<<"Calling algorithm read"<<endl;
  Algorithm.Read(in,1);
  
}


void PIMCClass::Run()
{
  cerr <<PathData.Path.Communicator.MyProc()<< " Simulation started." << endl;
  Algorithm.DoEvent();
  cerr<<PathData.Path.Communicator.MyProc()<<" PIMC++ has completed"<<endl;
  //  Array<MoveClass*,1> Moves;
//   for (int counter=0;counter<Moves.size();counter++){
//     cout<<"My name is "<<((MoveClass*)Moves(counter))->Name<<endl;
//     cout<<"My acceptance ratio is "<<((MoveClass*)Moves(counter))->AcceptanceRatio()<<endl;
//   }
  
}

void PIMCClass::Dummy()
{
	while(true)
		QMCWrapper->QMCDummy(PathData);
}

void PIMCClass::WriteSystemInfo()
{
  dVec box = PathData.Path.GetBox();
  Array<double,1> boxArray(3);
  boxArray(0) = box[0];
  boxArray(1) = box[1];
#if NDIM==3
  boxArray(2) = box[2];
#endif 
  OutFile.WriteVar ("Box", boxArray);
  OutFile.WriteVar("tau",PathData.Path.tau);
  OutFile.WriteVar("NumTimeSlices",PathData.Path.TotalNumSlices);
  OutFile.WriteVar("seed",PathData.Seed);
  for (int speciesIndex=0; speciesIndex < PathData.Path.NumSpecies(); 
       speciesIndex++) {
    SpeciesClass &species = PathData.Path.Species(speciesIndex);
    OutFile.NewSection("Species");
    OutFile.WriteVar ("Name", species.Name);
    OutFile.WriteVar ("NumParticles", species.NumParticles);
    OutFile.WriteVar ("lambda", species.lambda);
    ParticleType type = species.GetParticleType();
    if (type == FERMION)
      OutFile.WriteVar ("ParticleType", "Fermion");
    if (type == BOSON)
      OutFile.WriteVar ("ParticleType", "Boson");
    if (type == BOLTZMANNON)
      OutFile.WriteVar ("ParticleType", "Boltzmannon");
    OutFile.WriteVar("NEquilibrate", PathData.Path.NEquilibrate);
    OutFile.CloseSection(); //"Species"
  }
}
