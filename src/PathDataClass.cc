/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
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
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "PathDataClass.h"

void PathDataClass::Next(int &slice, int &ptcl)
{
  if (Join==slice)
    ptcl=Path.Permutation(ptcl);
  slice=(slice+1) % Path.NumTimeSlices();
}


void 
PathDataClass::WormInfo(int &headSlice, int &headPtcl,
		    int &tailSlice, int &tailPtcl,
		    int &numEmpty, int &wormSize)
{
  tailSlice=-1;tailPtcl=-1;headSlice=-1;headPtcl=-1;
  numEmpty=0;
  wormSize=0;
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
      if (SliceFullandPreviousSliceEmpty(slice,ptcl)){
	headSlice=slice;
	headPtcl=ptcl;
      }
      if (SliceFullandNextSliceEmpty(slice,ptcl)){
	tailSlice=slice;
	tailPtcl=ptcl;
      }
      if (Path.ParticleExist(slice,ptcl)==0 && slice!=Path.NumTimeSlices()-1)
	numEmpty++;
    }
  int currSlice=headSlice;
  int currPtcl=headPtcl;
  while (currSlice!=tailSlice || currPtcl!=tailPtcl){
    if (currSlice!=Path.NumTimeSlices()-1)
      wormSize++;
    Next(currSlice,currPtcl);
  } 
}
void
PathDataClass::MoveTailToSliceZero()
{
  int tailPtcl,tailSlice;
  int lastSlice=NumTimeSlices()-1;
  MoveJoin(lastSlice);
  FindTail(tailSlice,tailPtcl);
  int needToShift=tailSlice;
  ShiftData(-needToShift);
  Join=lastSlice-needToShift;


}


///Will only work in serial
void PathDataClass::MoveLinkToEnd(int linkToMove)
{
  //  cerr<<"Being called "<<linkToMove<<endl;
  
  int endTimeSlice=Path.NumTimeSlices()-1; //las slice on processor
  int needToShift=endTimeSlice-linkToMove;
  MoveJoin(linkToMove);
  if (needToShift!=0){
    ShiftData(needToShift);
    Join=linkToMove+needToShift;
    assert(Join==Path.NumTimeSlices()-1);
  }
  return;


}


///Will only work in serial
void PathDataClass::MoveOpenLinkToEnd()
{
  //  cerr<<"Calling moveopenlinktoend"<<endl;
  int openLink=Path.OpenLink;
  int endTimeSlice=Path.NumTimeSlices()-1; //las slice on processor
  int needToShift=endTimeSlice-openLink;
  //  cerr<<openLink<<" "<<endTimeSlice<<" "<<needToShift<<endl;
    MoveJoin(openLink);
  if (needToShift!=0){
    ShiftData(needToShift);
    Join=Path.OpenLink;
  }
  //  cerr<<"Join: "<<Join<<" "<<Path.OpenLink<<" "<<Path.NumTimeSlices()-1<<endl;
  assert(Join==Path.NumTimeSlices()-1);
  return;


}


/// Calculates the centroid position of each path
void PathDataClass::GetCentroids(Array<TinyVector<double,NDIM>,1> &CentPos)
{
  int N = Path.NumParticles();
  int M = Path.NumTimeSlices()-1;

  //Move the join to the end so we don't have to worry about permutations
  MoveJoin(M);

  // Loop through particles, making sure to visit each only once
  int initSlice = 0;
  bool usedPtcl[N];
  Array<TinyVector<double,NDIM>,1> tmpCentPos(N);
  for (int ptcl = 0; ptcl < N; ptcl++) {
    usedPtcl[ptcl] = false;
    tmpCentPos(ptcl) = 0.0;
  }
  for (int ptcl = 0; ptcl < N; ptcl++) {
    if (!usedPtcl[ptcl]) { // Only visit each particle once

      // Calculate centroid position
      int currSlice = initSlice;
      int currPtcl = ptcl;
      int nextSlice = -1;
      int nextPtcl = -1;
      dVec centPos = 0.0;
      int numSlices = 0;
      int numPtcls = 0;
      while (nextSlice != initSlice || nextPtcl != ptcl) {
        nextSlice = (currSlice + 1) % M;
        if (currSlice == Join) {
          nextPtcl = Path.Permutation(currPtcl);
          numPtcls++;
        } else
          nextPtcl = currPtcl;
        usedPtcl[nextPtcl] = true;
        dVec slicePos = Path(currSlice,currPtcl);
        //Path.PutInBox(slicePos); // HACK: DO I NEED THIS???
        centPos = centPos + slicePos;
        currSlice = nextSlice;
        currPtcl = nextPtcl;
        numSlices++;
      }
      tmpCentPos(ptcl) = centPos;
    }
  }

  // Gather up all the position data to get the centroids
  Path.Communicator.AllSum(tmpCentPos,CentPos);

  // Calculate the variance in each direction
  CentPos = CentPos/Path.TotalNumSlices;

}


///Worm Moves////////
bool PathDataClass::SliceFullandNextSliceEmpty(int slice,int ptcl)
{
  int nextSlice=(slice+1) % NumTimeSlices();
  int nextPtcl;
  if (Join==slice)
    nextPtcl=Path.Permutation(ptcl);
  else
    nextPtcl=ptcl;
  return (Path.ParticleExist(slice,ptcl)==1.0 && Path.ParticleExist(nextSlice,nextPtcl)==0.0);
}

bool PathDataClass::SliceFullandPreviousSliceEmpty(int slice,int ptcl)
{
  int prevSlice=((slice-1)+Path.NumTimeSlices() ) % Path.NumTimeSlices();
  int prevPtcl=0;
  if (Join==prevSlice){
    while (Path.Permutation(prevPtcl)!=ptcl)
      prevPtcl++;
  }
  else
    prevPtcl=ptcl;

  return (Path.ParticleExist(slice,ptcl)==1.0 && Path.ParticleExist(prevSlice,prevPtcl)==0.0);

}
  

void PathDataClass::FindHead(int &headSlice,int &headPtcl)
{
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++)
      if (SliceFullandPreviousSliceEmpty(slice,ptcl)){
	headSlice=slice;
	headPtcl=ptcl;
	return;
      }
}

void PathDataClass::FindTail(int &tailSlice,int &tailPtcl)
{
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++)
      if (SliceFullandNextSliceEmpty(slice,ptcl)){
	tailSlice=slice;
	tailPtcl=ptcl;
	return;
      }
}






void PathDataClass::Read (IOSectionClass &in)
{

#ifdef USE_QMC
  useDefaultStrings = true;
  int M = MetaWorldComm.NumProcs();
  // cerr << MetaWorldComm.MyProc() << ": PathDataClass.cc: N is " << M << endl;
  string QMCFilename;
  int ceimcProcs = 1;
  RUN_QMC = false;
  in.ReadVar("CEIMC_MODE", RUN_QMC);
  if(RUN_QMC){
    cerr << "CEIMC_MODE: Reading CEIMC Section..." << endl;
    assert (in.OpenSection ("CEIMC"));
    assert (in.ReadVar("QMCFile", QMCFilename));
    in.ReadVar("QMCPerCEIMC", ceimcProcs);
    assert(M%ceimcProcs == 0);
    if(!in.ReadVar("Timestep", dt))
      dt = 0.05;
    if(!in.ReadVar("Walkers", walkers))
      walkers = 10;
    if(!in.ReadVar("Steps", steps))
      steps = 100;
    if(!in.ReadVar("Blocks", blocks))
      blocks = 50;
    if(!in.ReadVar("Chains", chains))
      chains = 51;
    assert(in.ReadVar("Correlated",correlated));
    if(correlated){
      cerr << "Using correlated samping to compute energy differences." << endl;
      assert(in.ReadVar("QMCMethod",QMCMethod));
      if(QMCMethod == "VMC")
        cerr << "Using VMCMultiple driver." << endl;
      else if(QMCMethod == "RQMC")
        cerr << "Using RQMCMultiple driver." << endl;
      else {
        cerr << "Method " << QMCMethod << " not recognized.  Supported options are VMC or RQMC." << endl;
        assert(0);
      }
    }
    else
      cerr << "NOT using correlated sampling for energy differences." << endl;

    if(in.ReadVar("ParticleSet0", ptclSet0) && in.ReadVar("ParticleSet1", ptclSet1))
      useDefaultStrings = false;

    in.CloseSection();
    int argc = 1;
    char* argv[argc];
    argv[0]  = (char*)(QMCFilename.c_str());
    OHMMS::Controller->initialize(argc,argv);
    OhmmsInfo Welcome(argc,argv,OHMMS::Controller->mycontext());
    //qmcplusplus::QMCInterface* qmc = new qmcplusplus::QMCInterface(argc,argv);
    qmc = new qmcplusplus::QMCInterface(argc,argv);
    //qmcplusplus::QMCMain *qmcmain = new qmcplusplus::QMCMain(argc,argv);
    if(qmc->parse(argv[0])) {
      qmc->initialize(MetaWorldComm.MyProc(), MetaWorldComm.NumProcs());
    }
  }

  int managers = M/ceimcProcs;
  IAmQMCManager = false;
  Array<int, 1> WorldMembers(managers);
  int mgrIndex = 0;
  // Assemble list of managers, subset to be WorldComm
  for(int m=0; m<M; m++){
    if(m%ceimcProcs == 0){
      WorldMembers(mgrIndex) = m;
      mgrIndex++;
      if(MetaWorldComm.MyProc() == m)
        IAmQMCManager = true;
    }
  }
  int myCEIMCNum = MetaWorldComm.MyProc()/ceimcProcs;
  MetaWorldComm.Split(myCEIMCNum,QMCComm);
  cerr << MetaWorldComm.MyProc() << ": # of managers is " << managers << ".  IAmQMCManager is " << IAmQMCManager << " and I belong to CEIMC family " << myCEIMCNum << " of which I am proc " << QMCComm.MyProc() << endl;
  MetaWorldComm.Subset(WorldMembers, WorldComm);

  if(IAmQMCManager){
  int N = WorldComm.NumProcs();
    cerr << WorldComm.MyProc() << ": WorldComm size is " << N << endl;
  int procsPerClone = 1;
  if (N > 1) {
    assert (in.OpenSection ("Parallel"));
    assert (in.ReadVar("ProcsPerClone", procsPerClone));
    in.CloseSection();
  }

  ////  cerr << "  Set up Inter- and Intra- " << endl;
  // Setup Inter- and IntraComms
  assert ((N % procsPerClone) == 0);
  NumClones = N / procsPerClone;
  MyCloneNum = WorldComm.MyProc()/procsPerClone;
  // Create IntraComm
  //// cerr << "  Going to initialize IntraComm with MyCloneNum " << MyCloneNum << endl;
  WorldComm.Split(MyCloneNum, IntraComm);
  /// cerr << "  initialized IntraComm" << endl;  
  // cerr << "  skipped IntraComm" << endl;  
  Array<int,1> ranks (NumClones);
  for (int clone=0; clone<NumClones; clone++)
    ranks(clone) = clone*procsPerClone;
  /// cerr << "  ranks is " << ranks << "; going to creat Intercomm." << endl;
  WorldComm.Subset (ranks, InterComm);
  /// cerr << "  PIMC: initialized InterComm; ranks is " << ranks << endl;

  int seed;
  bool haveSeed = in.ReadVar ("Seed", seed);
  // Now, set up random number generator

  //  int seed;
  if (in.ReadVar("Seed",Seed)){
    Random.Init (Seed, NumClones);
  }
  else {
    Seed=Random.InitWithRandomSeed(NumClones);
  }
  //    Random.Init (314159, numClones);

  Path.MyClone=WorldComm.MyProc()/procsPerClone;
  }

#else
  // Has no function when PIMC++ is not built with qmcpack
  // but needs to be true to continue Path.Read (see PIMCClass.cc)
  // THIS IS THE DEFAULT!
  IAmQMCManager = true;
  int N = WorldComm.NumProcs();
  int procsPerClone = 1;
  bool sameSeed = false;
  if (N > 1) {
    assert (in.OpenSection("Parallel"));
    assert (in.ReadVar("ProcsPerClone",procsPerClone));
    in.ReadVar("SameLocalSeed",sameSeed);
    in.CloseSection();
  }

  // Setup Inter- and IntraComms
  assert ((N % procsPerClone) == 0);
  NumClones = N/procsPerClone;
  MyCloneNum = WorldComm.MyProc()/procsPerClone;
  WorldComm.Split(MyCloneNum, IntraComm);
  Array<int,1> ranks (NumClones);
  for (int clone=0; clone<NumClones; clone++)
    ranks(clone) = clone*procsPerClone;
  WorldComm.Subset(ranks, InterComm);

  // Now, set up random number generator
  int seed;
  bool haveSeed = in.ReadVar ("Seed", seed);
  if (in.ReadVar("Seed",Seed)){
    cout<<"RANDOM SEED: "<<Seed<<endl;
    Random.Init (Seed, NumClones, sameSeed);
  }
  else {
    Seed=Random.InitWithRandomSeed(NumClones);
  }

  // Assign Path clone variables
  Path.MyClone = MyCloneNum;
  Path.NumClones = NumClones;

  // Create CloneStr
  stringstream tempCloneStr;
  tempCloneStr << Path.MyClone << " " << Path.Communicator.MyProc();
  Path.CloneStr = tempCloneStr.str();

  if (WorldComm.MyProc() == 0) {
    cout <<"# Procs: "<<N<< ", # Clones: "<<NumClones<<", Procs/Clones: "
         <<procsPerClone<<", Root Host: "<<Path.Communicator.MyHost()<<endl;
  }

#endif

  // this is the NEW molecule stuff!!
  // initialize MoleculeHelperClass
  if(in.OpenSection("Molecules")){
    Mol.Read(in);
    in.CloseSection();
  }

  moveClock = 0;
}

#if USE_QMC
// Function for qmcpack functionality
void PathDataClass::AssignPtclSetStrings()
{
	if(!useDefaultStrings){
		assert(ptclSet0.size() == Path.NumSpecies());
		cerr << "  Got the following qmcpack particleset labels:\n";
	}
	else{
		ptclSet0.resize(Path.NumSpecies());
		ptclSet1.resize(Path.NumSpecies());
		for(int s=0; s< ptclSet0.size(); s++){
			ptclSet0(s) = Path.Species(s).Name;
			ptclSet1(s) = Path.Species(s).Name + '1';
		}
		cerr << "  Didn't read in particleset names; generated defaults:" << endl;
	}

	cerr << "SPECIES		SET 0		SET 1" << endl; 
	for(int s=0; s<ptclSet0.size(); s++) cerr << Path.Species(s).Name << "  " << ptclSet0(s) << "  " << ptclSet1(s) << endl;
}
#endif

void PathDataClass::MoveRefSlice (int absSlice)
{
  int first, last;
  int myProc   = Path.Communicator.MyProc();
  int numProcs = Path.Communicator.NumProcs();
  Path.SliceRange(numProcs-1, first, last);
  /// Min slices is the minimum number of slices owned by a processor
  /// It is the maximum we can shift the path at a time.
  int minSlices = last - first-1;
  Path.SliceRange(myProc, first, last);
  if (absSlice < Path.GetRefSlice()) { // Shift to left
    while ((Path.GetRefSlice()-absSlice) > minSlices) {
      MoveJoin(minSlices);
      ShiftData (-minSlices);
      Join = 0;
    }
    int shift = absSlice - Path.GetRefSlice();
    MoveJoin(-shift);
    ShiftData (shift);
    Join = 0;
  }
  else if (absSlice > Path.GetRefSlice()) { // Shift to right
    while ((absSlice-Path.GetRefSlice()) > minSlices) {
      MoveJoin(0);
      Path.ShiftData(minSlices);
      Join = minSlices;
    }
    int shift = absSlice - Path.GetRefSlice();
    MoveJoin(0);
    ShiftData(shift);
    Join = shift;
  }
  MoveJoin (Path.NumTimeSlices()-1);
  assert (Path.GetRefSlice() == absSlice);
}


#include <sys/time.h>

int
PathDataClass::GetWallTime()
{
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return (int)tv.tv_sec;
}


bool
PathDataClass::ExceededWallTime()
{
  if (MaxWallTime == -1)
    return false;
  return ((GetWallTime()-StartWallTime) > MaxWallTime);
}

void
PathDataClass::SetMaxWallTime(int maxWallTime)
{
  MaxWallTime = maxWallTime;
}

