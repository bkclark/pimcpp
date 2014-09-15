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

#include "PathClass.h"
#include "Actions/ActionsClass.h"
#include "IO/FileExpand.h"
#include <unistd.h>

/// This function initializes the paths depending on how they are
/// specified to be initialized in the input file.
void PathClass::InitPaths (IOSectionClass &in)
{
  SetMode(NEWMODE);
  assert(in.OpenSection ("Particles"));
  for (int speciesIndex=0; speciesIndex<NumSpecies(); speciesIndex++) {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(in.OpenSection("Species", speciesIndex));
    string InitPaths;
    in.ReadVar ("InitPaths", InitPaths);
    if (InitPaths == "LANGEVIN") {
      InitLangevin(in, species);
    } else if (InitPaths == "CUBIC") {
      InitCubic(in, species);
    } else if(InitPaths=="SPHERE") {
      InitSphere(in, species);
    } else if (InitPaths == "BCC") {
      InitBCC(in, species);
    } else if (InitPaths == "DIAG") {
      InitDiag(in, species);
    } else if (InitPaths == "LINE") {
      InitLine(in, species);
    } else if (InitPaths == "HYDROGENBCC") {
      InitHydrogenBCC(in, species);
    } else if (InitPaths=="ALLFIXED"){
      InitAllFixed(in, species);
    } else if (InitPaths == "FIXED") {
      InitFixed(in, species);
    } else if (InitPaths == "ADDVACANCIES") {
      InitAddVacancies(in, species);
    } else if (InitPaths == "ALLPATHS") {
      InitAllPaths(in, species);
    } else if (InitPaths=="RESTART") {
      InitRestart(in, species);
    } else if (InitPaths == "LEVIFLIGHT") {
      InitLeviFlight(in, species, speciesIndex);
    } else if (InitPaths == "RANDOM") {
      InitRandom(in, species);
    } else {
      perr << "Unrecognize initialization strategy " << InitPaths << endl;
      abort();
    }
    InitRealSlices();
    in.CloseSection(); // Species
  }
  in.CloseSection(); // "Particles"

  // Init open paths
  NowOpen = false;
  NowOpen.AcceptCopy();
  if (OpenPaths) {
    string openSpeciesName;
    assert(in.ReadVar("OpenSpecies",openSpeciesName));
    OpenSpeciesNum = SpeciesNum(openSpeciesName);
    InitOpenPaths();
  }

  // Init O(N)
  if (OrderN) {
    double cutoff;
    assert(in.ReadVar("CutoffDistance",cutoff));
    Cell.CutoffDistance=cutoff;
    Array<int,1>  numGrid;
    assert(in.ReadVar("CellGridSize",numGrid));
    assert(numGrid.size()==NDIM);
    numGrid(0)=25;
    numGrid(1)=25;
#if NDIM==3
    numGrid(2)=30;
#endif
    Cell.Init(Box,numGrid);
    for (int slice=0;slice<NumTimeSlices();slice++)
      Cell.BinParticles(slice);
  }

  // Calculate current sign
  SignWeight = GetSign();

  //Everything needs to be accepted
  Array<int,1> allParticles(NumParticles());
  for (int i=0; i<NumParticles(); i++)
    allParticles(i) = i;
  AcceptCopy(0, TotalNumSlices, allParticles);
  Path.AcceptCopy();
  BroadcastRefPath();
  RefPath.AcceptCopy();

  // Print particle coordinates
  //for (int i=0; i<NumParticles(); i++)
  //  for (int j=0; j<NumTimeSlices(); j++)
  //    cout << i << " " << j << " " << Path(j,i) << endl;
}

void PathClass::InitRestart(IOSectionClass &in, SpeciesClass &species)
{
  // Get restart h5 file
  string fileName;
  assert(in.ReadVar("File",fileName));

  // Decided whether or not to replicate
  bool replicate = false;
  in.ReadVar ("Replicate", replicate);

  int myProc = Communicator.MyProc();

  /// Decide whether to restart from several clones or just one (the 0th clone).
  bool parallelFileRead;
  IOSectionClass inFile;
  if (!in.ReadVar("ParallelFileRead",parallelFileRead))
    parallelFileRead = true;
  int tmpMyClone;
  if (!parallelFileRead)
    tmpMyClone = 0;
  else
    tmpMyClone = MyClone;

  /// If NumClones > MaxClones, reuse earlier .h5 files to restart
  int MaxClones;
  if (in.ReadVar("MaxClones",MaxClones)) {
    if (tmpMyClone > MaxClones-1)
      tmpMyClone = tmpMyClone % MaxClones;
  }

  /// Get previous clone filename
  stringstream oss;
  int counter = 0;
  oss<<fileName<<"."<<counter<<"."<<tmpMyClone<<".h5";
  while (fileExists(oss.str())) {
    counter++;
    oss.str("");
    oss<<fileName<<"."<<counter<<"."<<tmpMyClone<<".h5";
  }
  counter--;
  if (counter < 0) {
    cerr << CloneStr << " Previous clones for " << fileName << " don't exist! Aborting. " << endl;
    abort();
  }

  /// Try to find previous configuration to start from.
  /// Start with previous clone, then work way backwards.
  bool FoundPathDump = false;
  while(!FoundPathDump) {
    oss.str("");
    oss<<fileName<<"."<<counter<<"."<<tmpMyClone<<".h5";

    /// Check for Path Dump
    string fullFileName = oss.str();
    assert (inFile.OpenFile(fullFileName.c_str()));
    inFile.OpenSection("Observables");
    inFile.OpenSection("PathDump");
    IOVarBase *tmpPermVar = inFile.GetVarPtr("Permutation");
    /// Check if returning NULL pointer
    if (!tmpPermVar) {
      cerr << CloneStr << " WARNING: Path Dump not found in " << fullFileName << "! Trying previous clone..." << endl;
      counter--;
    } else {
      FoundPathDump = true;
    }
    inFile.CloseSection();
    inFile.CloseSection();
    inFile.CloseFile();
  }
  /// Broadcast so everyone agrees on the answer.
  Communicator.Broadcast(0,counter);

  /// Reopen the file
  oss.str("");
  oss<<fileName<<"."<<counter<<"."<<tmpMyClone<<".h5";
  string fullFileName = oss.str();
  assert (inFile.OpenFile(fullFileName.c_str()));

  /// Get the Box
  inFile.OpenSection("System");
  Array<double,1> oldBox;
  inFile.ReadVar("Box",oldBox);
  inFile.CloseSection();

  /// Get Permutations
  /// Only read in on processor 0 and then broadcast
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<int,1> oldPermutation;
  int extent0; int extent1; int extent2;
  if (myProc == 0) {
    IOVarBase *permutationVar = inFile.GetVarPtr("Permutation");
    int numDumps = permutationVar->GetExtent(0);
    permutationVar->Read(oldPermutation,numDumps-1,Range::all());
    extent0 = oldPermutation.extent(0);
  }
  Communicator.Broadcast(0,extent0);
  int numPerms = extent0;
  if (myProc != 0) {
    oldPermutation.resize(extent0);
  }
  Communicator.Broadcast(0,oldPermutation);
  /// Put permutation on proper last slice and accept
  SetMode(NEWMODE);
  if (myProc == Communicator.NumProcs()-1){
    for (int ptcl=0;ptcl<NumParticles();ptcl++)
      Permutation(ptcl) = oldPermutation(ptcl);
    Permutation.AcceptCopy();
  }

  /// Get Paths
  /// Only read in on processor 0 and then broadcast
  Array<double,3> oldPaths;
  IOVarBase *pathVar = inFile.GetVarPtr("Path");
  /// Check if returning NULL pointer
  if (!pathVar) {
    cerr << CloneStr << " Path Dump not found! Aborting " << endl;
    abort();
  }
  int numDumps=pathVar->GetExtent(0);
  if (myProc == 0){
    pathVar->Read(oldPaths,numDumps-1,Range::all(),Range::all(),Range::all());
    extent0 = oldPaths.extent(0);
    extent1 = oldPaths.extent(1);
    extent2 = oldPaths.extent(2);
  }
  Communicator.Broadcast(0,extent0);
  Communicator.Broadcast(0,extent1);
  Communicator.Broadcast(0,extent2);
  if (myProc != 0){
    oldPaths.resize(extent0,extent1,extent2);
  }
  Communicator.Broadcast(0,oldPaths);

  if (myProc == 0)
    cout<<CloneStr<<" Restarting from "<<oss.str()<<", Path Dumps: "<<numDumps<<" Permutations: "<<numPerms<<endl;

  /// Assign positions to beads
  int myFirstSlice,myLastSlice;
  SliceRange (myProc, myFirstSlice, myLastSlice);
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    int endSlice=min(myLastSlice,TotalNumSlices-1);
    for (int slice=0; slice<TotalNumSlices; slice++) {
      int relSlice = slice-myFirstSlice;
      if (myFirstSlice<=slice && slice<=myLastSlice){
        dVec pos;
        pos = 0.0;
        for (int dim=0; dim<NDIM; dim++)
          if (replicate) {
            pos(dim) = oldPaths(ptcl,0,dim)*(Box[dim]/oldBox(dim));
          } else if (slice>=oldPaths.extent(1)) {
            pos(dim) = oldPaths(ptcl,oldPaths.extent(1)-1,dim)*(Box[dim]/oldBox(dim));
          } else {
            pos(dim) = oldPaths(ptcl,slice,dim);//*(Box[dim]/oldBox(dim));
          }
        Path(relSlice,ptcl) = pos;
      }
    }
    /// If you are the last processors you must make sure the last
    /// slice is the same as the first slice on the first processors.
    /// You want to end up on the bead you permute onto.
    if (myProc==Communicator.NumProcs()-1) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        pos(dim) = oldPaths(Permutation(ptcl),0,dim);//*(Box[dim]/oldBox(dim));
      Path(NumTimeSlices()-1,ptcl) = pos;
    }
  }
  if (myProc==Communicator.NumProcs()-1)
    MoveJoin(NumTimeSlices()-1,1);

  /// Close things
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();

   ofstream outfile;
   outfile.open("positions.dat");
   for (int slice=0;slice<NumTimeSlices();slice++)
     for (int ptcl=0;ptcl<NumParticles();ptcl++)
       outfile<<slice<<" "<<ptcl<<"  "<<Path(slice,ptcl)[0]<<" "<<Path(slice,ptcl)[1]<<" "<<Path(slice,ptcl)[2]<<endl;
   outfile.close();
}

void PathClass::InitCubic(IOSectionClass &in, SpeciesClass &species)
{
  int num = species.NumParticles;
  bool isCubic = 0;
#if NDIM==2
  isCubic = (Box[0]==Box[1]);
#endif
#if NDIM==3
  isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
  if (!isCubic) {
    perr << isCubic << " A cubic box is current required for cubic initilization\n";
    abort();
  }
  int numPerDim = (int) ceil (pow((double)num, 1.0/NDIM)-1.0e-6);
  double delta = Box[0] / numPerDim;
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    int ip = (ptcl-species.FirstPtcl);
    int ix, iy, iz;
#if NDIM==2
    ix = ip/(numPerDim);
    iy = ip-(ix*numPerDim);
#endif
#if NDIM==3
    ix = ip/(numPerDim*numPerDim);
    iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
    iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
#endif
    dVec r;
    r[0] = ix*delta-0.5*Box[0];
    r[1] = iy*delta-0.5*Box[1];
#if NDIM==3
    r[2] = iz*delta-0.5*Box[2];
#endif
    for (int slice=0; slice<NumTimeSlices(); slice++)
      Path(slice,ptcl) = r;
    dVec disp = r;
  }
}

void PathClass::InitSphere(IOSectionClass &in, SpeciesClass &species)
{
  dVec r0,r;
  species = *SpeciesArray(0);
  double SphereRadius;
  assert(in.ReadVar("SphereRadius",SphereRadius));
  double sigma = sqrt(2*species.lambda*tau);
  for (int ptcl=0; ptcl<NumParticles(); ptcl++){
    Random.LocalGaussianVec(1,r0);
    for (int slice=0; slice<NumTimeSlices(); slice++){
      Random.LocalGaussianVec(sigma,r);
      r = r0 + r;
      SetPos(slice,ptcl,r*SphereRadius*(1./sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])));
    }
  }
}

void PathClass::InitBCC(IOSectionClass &in, SpeciesClass &species)
{
  int num = species.NumParticles;
  bool isCubic = 0;
#if NDIM==2
  isCubic = (Box[0]==Box[1]);
#endif
#if NDIM==3
  isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
  if (!isCubic) {
    perr << "A cubic box is current required for cubic initilization\n";
    abort();
  }
  int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/NDIM)-1.0e-6);
  double delta = Box[0] / numPerDim;
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    int ip = (ptcl-species.FirstPtcl)/2;
    int ix, iy, iz;
#if NDIM==2
    ix = ip/(numPerDim);
    iy = ip-(ix*numPerDim);
#endif
#if NDIM==3
    ix = ip/(numPerDim*numPerDim);
    iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
    iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
#endif
    dVec r;
    r[0] = ix*delta-0.5*Box[0];
    r[1] = iy*delta-0.5*Box[1];
#if NDIM==3
    r[2] = iz*delta-0.5*Box[2];
#endif
    if (ptcl % 2)
      r += 0.5*delta;
    for (int slice=0; slice<NumTimeSlices(); slice++)
      Path(slice,ptcl) = r;
  }
}

void PathClass::InitDiag(IOSectionClass &in, SpeciesClass &species)
{
  // Diagonal paths
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    dVec r;
    r[0] = ptcl;
    r[1] = ptcl;
#if NDIM==3
    r[2] = ptcl;
#endif
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec temp = r*(1. + 0.1/sqrt(3.0*ptcl*ptcl));
      Path(slice,ptcl) = temp;
    }
  }
}

void PathClass::InitLine(IOSectionClass &in, SpeciesClass &species)
{
  // Forming a straight line
  int num = species.NumParticles;
  double delta = 1.0;
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    dVec r;
    r[0] = ptcl*delta;
    r[1] = ptcl*delta;
#if NDIM==3
    r[2] = ptcl*delta;
#endif
    for (int slice=0; slice<NumTimeSlices(); slice++)
      Path(slice,ptcl) = r;
  }
}

void PathClass::InitHydrogenBCC(IOSectionClass &in, SpeciesClass &species)
{
  int num = species.NumParticles;
  bool isCubic = 0;
#if NDIM==2
  isCubic = (Box[0]==Box[1]);
#endif
#if NDIM==3
  isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
  if (!isCubic) {
    perr << "A cubic box is current required for cubic initilization\n";
    abort();
  }
  int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/NDIM)-1.0e-6);
  double delta = Box[0] / numPerDim;
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    int ip = (ptcl-species.FirstPtcl)/2;
    int ix, iy, iz;
#if NDIM==2
    ix = ip/(numPerDim);
    iy = ip-(ix*numPerDim);
#endif
#if NDIM==3
    ix = ip/(numPerDim*numPerDim);
    iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
    iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
#endif
    dVec r;
    r[0] = ix*delta-0.5*Box[0];
    r[1] = iy*delta-0.5*Box[1];
#if NDIM==3
    r[2] = iz*delta-0.5*Box[2];
#endif
    if (ptcl % 2)
      r += 0.5*delta;
    if ( species.isIon ) {
      r[0] += 0.5;
    }
    for (int slice=0; slice<NumTimeSlices(); slice++)
      Path(slice,ptcl) = r;
  }

}

void PathClass::InitAllFixed(IOSectionClass &in, SpeciesClass &species)
{
  int myFirstSlice, myLastSlice, myProc;
  myProc = Communicator.MyProc();
  SliceRange (myProc, myFirstSlice, myLastSlice);
  Array<double,3> Positions;
  assert(in.ReadVar("Positions",Positions));
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<species.NumDim; dim++)
        pos(dim) = Positions(slice,ptcl-species.FirstPtcl,dim);
      Path(slice,ptcl) = pos;
    }

    ///If you are the last processors you must make sure the last
    ///slice is the same as the first slice on the first
    ///processors. The  join should be at the
    if (myProc==Communicator.NumProcs()-1){
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        pos(dim) = Positions(0,ptcl-species.FirstPtcl,dim);
      Path(NumTimeSlices()-1,ptcl) = pos;
    }
  }
}

void PathClass::InitFixed(IOSectionClass &in, SpeciesClass &species)
{
  Array<double,2> Positions;
  assert (in.ReadVar ("Positions", Positions));
  assert (Positions.rows() == species.NumParticles);
  assert (Positions.cols() == species.NumDim);
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<species.NumDim; dim++)
        pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
      Path(slice,ptcl) = pos;
    }
  }
}

void PathClass::InitAddVacancies(IOSectionClass &in, SpeciesClass &species)
{
  Array<double,2> Positions;
  assert (in.ReadVar ("Positions", Positions));
  Positions.resizeAndPreserve(species.NumParticles,Positions.extent(1));
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++){
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<species.NumDim; dim++)
        pos(dim) = Positions(ptcl-species.FirstPtcl,dim)*ScaleBox;
      Path(slice,ptcl) = pos;
    }
  }
}

void PathClass::InitAllPaths(IOSectionClass &in, SpeciesClass &species)
{
  Array<double,3> Positions;
  assert (in.ReadVar ("Positions", Positions));
  perr<<"NumTimeSlices: "<<TotalNumSlices<<" "<<Positions.extent(0);
  assert (Positions.extent(0) == species.NumParticles);
  assert (Positions.extent(1) == TotalNumSlices);
  assert (Positions.extent(2) == species.NumDim);
  for (int ptcl=species.FirstPtcl; 
       ptcl<=species.LastPtcl; ptcl++){
    for (int slice=0; slice<TotalNumSlices; slice++) {
      perr<<ptcl;
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<species.NumDim; dim++)
        pos(dim) = Positions(ptcl-species.FirstPtcl,slice,dim);
      Path(slice,ptcl) = pos;
    }
    int slice=NumTimeSlices()-1;
    dVec pos;
    pos = 0.0;
    for (int dim=0; dim<species.NumDim; dim++)
      pos(dim) = Positions(ptcl-species.FirstPtcl,0,dim);
    Path(slice,ptcl) = pos;
  }
}

void PathClass::InitLeviFlight(IOSectionClass &in, SpeciesClass &species, int speciesIndex)
{
  Array<double,2> Positions(species.NumParticles,species.NumDim);
  double sigmaFactor = 1.0;
  if(!in.ReadVar ("Positions", Positions)) {
    InitBCC(in,species);
    for (int ptcl=0; ptcl<species.NumParticles; ptcl++) 
      for (int dim=0; dim<NDIM; dim++)
        Positions(ptcl,dim) = Path(0,ptcl)[dim];
  }
  if (in.ReadVar ("SigmaFactor", sigmaFactor, 1.0))
    perr << "SigmaFactor = " << sigmaFactor << endl;
  assert (Positions.rows() == species.NumParticles);
  assert (Positions.cols() == species.NumDim);
  Array<dVec,1> R0(species.NumParticles);
  for (int ptcl=0; ptcl<species.NumParticles; ptcl++) 
    for (int dim=0; dim<NDIM; dim++)
      R0(ptcl)[dim] = Positions(ptcl,dim);

  bool haveNodeAction = Actions.NodalActions(speciesIndex)!=NULL;
  bool isNodeAvoiding = false;
  in.ReadVar ("NodeAvoiding", isNodeAvoiding);
  bool isPhaseAvoiding = false;
  in.ReadVar ("PhaseAvoiding", isPhaseAvoiding);
  if (isNodeAvoiding||haveNodeAction)
    NodeAvoidingLeviFlight(speciesIndex, R0);
  else if (isPhaseAvoiding)
    PhaseAvoidingLeviFlight(speciesIndex, R0, sigmaFactor);
  else {
    int myStart, myEnd;
    SliceRange (Communicator.MyProc(), myStart, myEnd);
    Array<double,2> Positions;
    assert (in.ReadVar ("Positions", Positions));
    assert (Positions.rows() == species.NumParticles);
    assert (Positions.cols() == species.NumDim);
    Array<dVec,1> flight(TotalNumSlices+1);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      for (int dim=0; dim<NDIM; dim++) {
        flight(0)[dim] = Positions(ptcl-species.FirstPtcl,dim);
        flight(TotalNumSlices)[dim] = Positions(ptcl-species.FirstPtcl,dim);
      }
      LeviFlight (flight, species.lambda);
//    for (int i=0; i<MyNumSlices; i++)
//      Path(i, ptcl) = flight(i-RefSlice);
      for (int slice=myStart; slice<=myEnd; slice++)
        Path(slice-myStart, ptcl) = flight(slice);
    }
  }
}

void PathClass::InitRandom(IOSectionClass &in, SpeciesClass &species)
{
  double radius;
  if(!in.ReadVar("Radius", radius))
    radius = Box[0]/species.NumParticles;
  int N = species.NumParticles;
  Array<dVec,1> R(N);
  for (int i=0; i<N; i++) {
    bool overlap = true;
    int tries = 0;
    while ((overlap) && (tries < 1000)) {
      tries++;
      for (int dim=0; dim<NDIM; dim++) 
        R(i)[dim] = Box[dim]*(Random.World()-0.5);
      overlap = false;
      for (int j=0; j<i; j++) {
        dVec disp = R(j)-R(i);
        PutInBox(disp);
        overlap = overlap || (dot(disp,disp) < (4.0*radius*radius));
      }
    }
    if (tries == 1000) {
      cerr << "Exceed 1000 tries in InitRandom.  "
           << "Decrease excluded radius.\n";
      abort();
    }
  }
  for (int slice=0; slice<NumTimeSlices(); slice++) {
    for (int ptcl=0; ptcl<N; ptcl++) {
      Path[0](slice, ptcl + species.FirstPtcl) = R(ptcl);
      Path[1](slice, ptcl + species.FirstPtcl) = R(ptcl);
    }
  }
}

/// Constructs a Levi flight beginning in the vec(0) and ending
/// in vec(N-1) if vec has length N.  This is a path which samples
/// exactly the free particle density matrix.
void PathClass::LeviFlight (Array<dVec,1> &vec, double lambda)
{
  // HACK HACK HACK
  //  tau *= 0.000001;

  int N = vec.size();
  for (int slice=1; slice<(N-1); slice++) {
    double di = (double)slice;
    double delta = (double)(N-slice-1);
    dVec center = (1.0/(delta+1.0))*(delta*vec(slice-1) + vec(N-1));
    double taueff = tau*(1.0 - 1.0/delta);
    double sigma = sqrt (2.0*lambda*taueff);
    Random.CommonGaussianVec(sigma, vec(slice));
    vec(slice) += center;
  }
  // DEBUG DEBUG DEBUG DEBUG DEBUG!!!!
  FILE *fout = fopen ("Levi.dat", "w");
  for (int i=0; i<N; i++) {
    for (int dim=0; dim<NDIM; dim++) 
      fprintf (fout, "%1.12e ", vec(i)[dim]);
    fprintf (fout, "\n");
  }
  fclose (fout);
}

void PathClass::NodeAvoidingLeviFlight (int speciesNum, Array<dVec,1> &R0)
{
  SpeciesClass &species = Species(speciesNum);
  double lambda = species.lambda;
  bool haveNodeAction = Actions.NodalActions(speciesNum)!=NULL;

  int myFirstSlice, myLastSlice, myProc;
  myProc = Communicator.MyProc();
  SliceRange (myProc, myFirstSlice, myLastSlice);

  // HACK to get ground state plane wave calculations to happen
  // simultaneously rather than sequentially.
  if (haveNodeAction) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      (*this)(0,ptcl) = R0(ptcl-species.FirstPtcl);
    Actions.NodalActions(speciesNum)->IsPositive(0);
  }

  int numPtcls = species.NumParticles;
  Array<dVec,1> prevSlice (numPtcls), newSlice(numPtcls);
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    prevSlice (ptcl) = R0(ptcl);
    if (myProc == 0)
      (*this)(0, ptcl+species.FirstPtcl) = R0(ptcl);
    RefPath(ptcl+species.FirstPtcl) = R0(ptcl);
  }

  /// Check slice 0 -- important for ground state nodes
  bool swapFirst = false;
  if (myProc == 0)
    if (haveNodeAction)
      if (!Actions.NodalActions(speciesNum)->IsPositive(0)) {
        perr << "Initially negative node action.  Swapping two particles "
             << "for species " << species.Name << ".\n";
        swapFirst = true;
      }
  Communicator.Broadcast (0,swapFirst);

  if (swapFirst) {
    dVec tmp = R0(0);
    R0(0) = R0(1);
    R0(1) = tmp;
    RefPath(species.FirstPtcl)     = R0(0);
    RefPath(species.FirstPtcl+1)   = R0(1);
    prevSlice(0)                   = R0(0);
    prevSlice(1)                   = R0(1);
    if (myProc==0) {
      (*this)(0,species.FirstPtcl)   = R0(0);
      (*this)(0,species.FirstPtcl+1) = R0(1);  
      if (!Actions.NodalActions(speciesNum)->IsPositive(0)) {
        perr << "Still not positive after swap!!!!!!!!!!!!\n";
        /// HACK HACK HACK -- commenting out abort for now to allow
        /// fixed-phase to continue.
        // abort();
      }
    }
  }

  for (int slice=1; slice<TotalNumSlices+1; slice++) {
    int sliceOwner = SliceOwner(slice);
    int relSlice = slice-myFirstSlice;
    double delta = (double)(TotalNumSlices-slice);
    double taueff = tau*(1.0 - 1.0/(delta+1.0));
    double sigma = sqrt (2.0*lambda*taueff);
    bool positive = false;
    int numRejects = 0;
    do {
      // Randomly construct new slice
      for (int ptcl=0; ptcl<numPtcls; ptcl++) {
        dVec center = (1.0/(delta+1.0))*(delta*prevSlice(ptcl) + R0(ptcl));
        Random.CommonGaussianVec(sigma, newSlice(ptcl));
        newSlice(ptcl) += center;
      }
      // Now check the nodal sign if we're a fermion species
      if (!haveNodeAction)
        positive = true;
      else {
        // Now assign to Path
        if (sliceOwner == myProc) {
          for (int ptcl=0; ptcl<numPtcls; ptcl++)
            (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
          positive = Actions.NodalActions(speciesNum)->IsPositive(relSlice);
        }
        // Now broadcast whether or not I'm positive to everyone
        Communicator.Broadcast(sliceOwner, positive);
      }
      if (!positive) {
        numRejects++;
        if ((numRejects%10)==0)
          cerr << "numRejects = " << numRejects << endl;
      }
    } while (!positive);
    // Copy slice into Path if I'm the slice owner.
    if ((slice>=myFirstSlice) && (slice<=myLastSlice)) {
      for (int ptcl=0; ptcl<numPtcls; ptcl++)
        (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
      // Check to make sure we're positive now.
      if (haveNodeAction) {
        if (!Actions.NodalActions(speciesNum)->IsPositive(relSlice)) {
          perr << "Still not postive at slice " << slice << " myProc = " << myProc << "relslice=" << relSlice <<endl;
          abort();
        }
      }
    }
    // continue on to next slice
    prevSlice = newSlice;
  }
  Array<int,1> changedParticles(species.NumParticles);
  if (haveNodeAction) {
    for (int i=0; i<species.NumParticles; i++)
      changedParticles(i) = i+species.FirstPtcl;
    double localAction = Actions.NodalActions(speciesNum)->Action(0, NumTimeSlices()-1, changedParticles, 0);
    Communicator.PrintSync();
    double globalAction = Communicator.AllSum(localAction);
    if (Communicator.MyProc()==0)
      perr << "Nodal Action after Levi flight = " << globalAction << endl;
  }
  Actions.NodalActions(speciesNum) -> AcceptCopy(0, NumTimeSlices()-1);
  AcceptCopy(0, NumTimeSlices()-1, changedParticles);

}

void PathClass::PhaseAvoidingLeviFlight (int speciesNum, Array<dVec,1> &R0, double sigmaFactor)
{
  SpeciesClass &species = Species(speciesNum);
  int numPtcls = species.NumParticles;
  double lambda = species.lambda;
  bool haveNodeAction = Actions.NodalActions(speciesNum)!=NULL;
  double maxAction = 0.1*M_PI*M_PI/(double)numPtcls;


  int myFirstSlice, myLastSlice, myProc;
  myProc = Communicator.MyProc();
  SliceRange (myProc, myFirstSlice, myLastSlice);

  // HACK to get ground state plane wave calculations to happen
  // simultaneously rather than sequentially.
  if (haveNodeAction) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      (*this)(0,ptcl) = R0(ptcl-species.FirstPtcl);
    Actions.NodalActions(speciesNum)->IsPositive(0);
  }

  Array<dVec,1> prevSlice (numPtcls), newSlice(numPtcls);
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    prevSlice (ptcl) = R0(ptcl);
    if (myProc == 0)
      (*this)(0, ptcl+species.FirstPtcl) = R0(ptcl);
    RefPath(ptcl+species.FirstPtcl) = R0(ptcl);
  }
  
  int N = TotalNumSlices+1;
  for (int slice=1; slice<N; slice++) {
    int sliceOwner = SliceOwner(slice);
    int relSlice = slice-myFirstSlice;
    
    double delta = (double)(N-slice-1);
    double taueff = tau*(1.0 - 1.0/(delta+1.0));
    double sigma = sigmaFactor*sqrt (2.0*lambda*taueff);
    bool positive = false;
    
    int numRejects = 0;
    do {
      // Randomly construct new slice
      for (int ptcl=0; ptcl<numPtcls; ptcl++) {
        dVec center = (1.0/(delta+1.0))*(delta*prevSlice(ptcl) + R0(ptcl));
        Random.CommonGaussianVec(sigma, newSlice(ptcl));
        newSlice(ptcl) += center;
      }
      // Now check the nodal sign if we're a fermion species
      if (!haveNodeAction)
        positive = true;
      else {
        // Now assign to Path
        if (sliceOwner == myProc ) {
          for (int ptcl=0; ptcl<numPtcls; ptcl++)
            (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
          double action = Actions.NodalActions(speciesNum)->SingleAction
            (relSlice-1, relSlice, species.Ptcls, 0);
          positive = action < maxAction;
        }
        // Now broadcast whether or not I'm positive to everyone
        Communicator.Broadcast(sliceOwner, positive);
      }
      if (!positive) {
        numRejects++;
        if (numRejects > 200) {
          /// Create a new starting point and start over
          for (int i=0; i<R0.size(); i++)
            for (int dim=0; dim<NDIM; dim++)
              R0(i)[dim] = Box[dim] *(Random.Common()-0.5);
          perr << "Calling recursively to start over.\n";
          PhaseAvoidingLeviFlight (speciesNum, R0, sigmaFactor);
          return;
        }
//      if ((numRejects%50)==0)
//        perr << "numRejects = " << numRejects << endl;
      }
    } while (!positive);
    // Copy slice into Path if I'm the slice owner.
    if ((slice>=myFirstSlice) && (slice<=myLastSlice)) {
      for (int ptcl=0; ptcl<numPtcls; ptcl++) 
        (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
      // Check to make sure we're positive now.
      if (haveNodeAction && (relSlice > 0)) {
        if (Actions.NodalActions(speciesNum)->Action 
            (relSlice-1, relSlice, species.Ptcls, 0) > maxAction) {
          perr << "exceeded maximum action at slice " << slice 
               << " myProc = " << myProc << "relslice=" << relSlice <<endl;
          abort();
        }       
      }
    }
    // continue on to next slice
    prevSlice = newSlice;
  }
  if (haveNodeAction) {
    Array<int,1> changedParticles(species.NumParticles);
    for (int i=0; i<species.NumParticles; i++)
      changedParticles(i) = i+species.FirstPtcl;
    double localAction = 
      Actions.NodalActions(speciesNum)->Action(0, NumTimeSlices()-1, 
                                               changedParticles,0);
    double globalAction = Communicator.AllSum(localAction);
    if (Communicator.MyProc()==0)
      perr << "Phase Action after Levi flight = " << globalAction << endl;
  }
  if (Actions.NodalActions(speciesNum) != NULL)
    Actions.NodalActions(speciesNum)->AcceptCopy(0, NumTimeSlices()-1);  
}

void PathClass::InitLangevin(IOSectionClass &in, SpeciesClass &species)
{
  // Read the name of the previous runs output file
  string langevinName;
  assert (in.ReadVar("LangevinFile", langevinName));
  langevinName = ExpandFileName (langevinName);

  // Read in data from previous runs output file
  IOSectionClass langFile;
  assert (langFile.OpenFile (langevinName));
  assert (langFile.OpenSection("Moves"));
  assert (langFile.OpenSection("Langevin"));

  int first = species.FirstPtcl;
  int last  = species.LastPtcl;
  IOVarBase *Rvar = langFile.GetVarPtr("R");
  assert (Rvar != NULL);
  assert (Rvar->GetRank() == 3);
  int numSteps = Rvar->GetExtent(0);
  Array<double,2> lastPos(last-first+1,NDIM);
  Rvar->Read(lastPos, (numSteps-1), Range::all(), Range::all());
  langFile.CloseFile();

  // Now assign to path
  for (int slice=0; slice<NumTimeSlices(); slice++) 
    for (int ptcl=first; ptcl<=last; ptcl++) 
      for (int dim=0; dim<NDIM; dim++)
        Path(slice,ptcl)[dim] = lastPos(ptcl-first,dim);

  perr << "Successfully read Langevin initial positions."  << endl;
}
