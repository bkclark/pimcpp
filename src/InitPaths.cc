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


void PathClass::Restart(IOSectionClass &in, string fileName, bool replicate, SpeciesClass &species)
{
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
            pos(dim) = oldPaths(ptcl,slice,dim)*(Box[dim]/oldBox(dim));
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
        pos(dim) = oldPaths(Permutation(ptcl),0,dim)*(Box[dim]/oldBox(dim));
      Path(NumTimeSlices()-1,ptcl) = pos;
    }
  }
  if (myProc==Communicator.NumProcs()-1)
    MoveJoin(NumTimeSlices()-1,1);

  /// Close things
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();

//   ofstream outfile;
//   outfile.open("positions.dat");
//   for (int slice=0;slice<NumTimeSlices();slice++)
//     for (int ptcl=0;ptcl<NumParticles();ptcl++)
//       outfile<<slice<<" "<<ptcl<<"  "<<Path(slice,ptcl)[0]<<" "<<Path(slice,ptcl)[1]<<endl;
//   outfile.close();
}


void PathClass::ReadSqueeze(IOSectionClass &in,string fileName,bool replicate)
{
// This was modified on Jan 19 2005 to read in a set of classical (P=1) configs and duplicate them to produce a set of PIMC (P>1) configs.  -jg

///  cerr<<"Read squeezing now"<<endl;
  IOSectionClass inFile;
  stringstream oss;
  ////  oss<<fileName<<"."<<MyClone<<".h5";
  bool parallelFileRead;
  if (!in.ReadVar("ParallelFileRead",parallelFileRead)){
    parallelFileRead=true;
  }
  if (parallelFileRead)
    oss<<fileName<<"."<<MyClone<<".h5";
  else
    oss<<fileName<<"."<<0<<".h5";
  string fullFileName=oss.str();
  ////  cerr<<"THE FULL FILE NAME IS "<<fullFileName<<endl;
  assert (inFile.OpenFile(fullFileName.c_str()));
  //  assert (inFile.OpenFile(fileName.c_str()));
  inFile.OpenSection("System");
  Array<double,1> oldBox;
  inFile.ReadVar("Box",oldBox);
  inFile.CloseSection();
  ///  cerr<<"Read the box"<<endl;
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,3> oldPaths; //(58,2560,2,3);
  Array<int,2> oldPermutation;
  assert(inFile.ReadVar("Permutation",oldPermutation));
  SetMode(NEWMODE);;
  int myProc=Communicator.MyProc();
  if (myProc==Communicator.NumProcs()-1){
    for (int ptcl=0;ptcl<NumParticles();ptcl++)
      Permutation(ptcl) = oldPermutation(oldPermutation.extent(0)-1,ptcl);
    Permutation.AcceptCopy();
  }
    IOVarBase *pathVar = inFile.GetVarPtr("Path");
  int numDumps=pathVar->GetExtent(0);
  cerr<<"Number of path dumps is "<<numDumps<<endl;
  cerr<<"Number of permutations is "<<oldPermutation.extent(0)<<endl;
  pathVar->Read(oldPaths,numDumps-1,Range::all(),Range::all(),Range::all());
  cerr<<"My data is now read"<<endl;
  //  assert(inFile.ReadVar("Path",oldPaths));
  double sliceRatio=(double)oldPaths.extent(1)/(double)TotalNumSlices;
  cerr<<"Sliceratio is "<<sliceRatio<<endl;
  ///  cerr << "My paths are of size"  << oldPaths.extent(0) << " "
  ///       << oldPaths.extent(1)<<" " << oldPaths.extent(2) << endl;

  int myFirstSlice,myLastSlice;
  SliceRange (myProc, myFirstSlice, myLastSlice);  
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    for (int slice=0; slice<TotalNumSlices; slice++) {
      int sliceOwner = SliceOwner(slice);
      int relSlice = slice-myFirstSlice;
      ///      if (myProc==sliceOwner){
      if (myFirstSlice<=slice && slice<=myLastSlice){
        dVec pos;
        pos = 0.0;
        for (int dim=0; dim<NDIM; dim++)
          if (replicate){
            pos(dim) = oldPaths(ptcl,0,dim)*(Box[dim]/oldBox(dim));
          }
        //UNHACK          else if (slice>=oldPaths.extent(2)){
          else if (sliceRatio*slice>=oldPaths.extent(1)){
            pos(dim)=oldPaths(ptcl,oldPaths.extent(1)-1,dim)*(Box[dim]/oldBox(dim));
          }
          else{
            //UNHACK        pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,slice,dim)*(Box[dim]/oldBox(dim));
            pos(dim) = oldPaths(ptcl,(int)trunc(sliceRatio*slice),dim)*(Box[dim]/oldBox(dim));
          }

        Path(relSlice,ptcl) = pos;
      }
    }      
    ///If you are the last processors you must make sure the last
    ///slice is the same as the first slice on the first
    ///processors. The  join should be at the 
    if (myProc==Communicator.NumProcs()-1){
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        if (replicate){//not sure this works for replicate
          pos(dim) = oldPaths(Permutation(ptcl),0,dim)*(Box[dim]/oldBox(dim));
        }
        else{//you want to end up on the person you permute onto
          pos(dim) = oldPaths(Permutation(ptcl),0,dim)*(Box[dim]/oldBox(dim));
        }
    
      Path(NumTimeSlices()-1,ptcl) = pos;
    }
  }
  cerr<<"Finished with particles"<<endl;
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();
  if (myProc==Communicator.NumProcs()-1)
    MoveJoin(NumTimeSlices()-1,1);
}



void 
PathClass::Tile(IOSectionClass &in,string fileName,bool replicate)
{
  //Invariants: There are no permutations.  
  //Invariants: NDIM==2
  cerr<<"starting tile " <<Communicator.MyProc()<<endl;
  IOSectionClass inFile;
  stringstream oss;
  bool parallelFileRead;
  if (!in.ReadVar("ParallelFileRead",parallelFileRead)){
    parallelFileRead=true;
  }
  if (parallelFileRead)
    oss<<fileName<<"."<<MyClone<<".h5";
  else
    oss<<fileName<<"."<<0<<".h5";
  string fullFileName=oss.str();
  assert (inFile.OpenFile(fullFileName.c_str()));
  inFile.OpenSection("System");
  Array<double,1> oldBox;
  inFile.ReadVar("Box",oldBox);
  inFile.CloseSection();
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,3> oldPaths; //(58,2560,2,3);
  Array<int,2> oldPermutation;
  assert(inFile.ReadVar("Permutation",oldPermutation));
  SetMode(NEWMODE);
  int myProc=Communicator.MyProc();
  for (int ptcl=0;ptcl<NumParticles();ptcl++)
    Permutation(ptcl) = ptcl; //oldPermutation(oldPermutation.extent(0)-1,ptcl);
  Permutation.AcceptCopy();
  cerr<<"About to read paths"<<endl;
  IOVarBase *pathVar = inFile.GetVarPtr("Path");
  int numDumps=pathVar->GetExtent(0);
  cerr<<"Number of path dumps is "<<numDumps<<endl;
  pathVar->Read(oldPaths,numDumps-1,Range::all(),Range::all(),Range::all());
  double sliceRatio=(double)oldPaths.extent(1)/(double)TotalNumSlices;
  assert(fabs(sliceRatio-1)<1e-10);
  int oldParticles;
  oldParticles=oldPaths.extent(0);
  int numTiles=NumParticles()/oldParticles;
  cerr<<"numTiles is "<<numTiles<<endl;
  cerr<<"oldParticles is "<<oldParticles<<endl;
  cerr<<"OldBox is "<<oldBox<<endl;
  dVec oldBox_dVec;
  for (int dim=0;dim<NDIM;dim++)
    oldBox_dVec(dim)=oldBox(dim);
  cerr<<"Box is "<<Box<<endl;
  cerr<<"Sliceratio is "<<sliceRatio<<endl;

  int myFirstSlice,myLastSlice;

  Array<dVec,1> first_slice(NumParticles());
  SliceRange (myProc, myFirstSlice, myLastSlice);  
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    dVec prev_pos;
    for (int slice=0; slice<TotalNumSlices; slice++) {
      int relSlice = slice-myFirstSlice;

      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        pos(dim) = oldPaths(ptcl % oldParticles,(int)trunc(sliceRatio*slice),dim); //(Box[dim]/oldBox(dim));
      PutInBox(pos,oldBox_dVec);
      if (ptcl/oldParticles==1 || ptcl/oldParticles==3)
        pos(0)=pos(0)+oldBox(0);
      if (ptcl/oldParticles==2 || ptcl/oldParticles==3)
        pos(1)=pos(1)+oldBox(1);
      if (slice!=0){
        dVec diff=pos-prev_pos;
        PutInBox(diff,oldBox_dVec);
        pos=diff+prev_pos;
      }
      assert(NDIM==2);
      prev_pos=pos;
      if (slice==0)
        first_slice(ptcl)=pos;
      if (myFirstSlice<=slice && slice<=myLastSlice){
        Path(relSlice,ptcl) = pos;
      }
    }
    ///If you are the last processors you must make sure the last
    ///slice is the same as the first slice on the first
    ///processors. The  join should be at the 
    if (myProc==Communicator.NumProcs()-1){
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        pos(dim) = oldPaths(ptcl % oldParticles,0,dim); //*(Box[dim]/oldBox(dim));
      PutInBox(pos,oldBox_dVec);
      if (ptcl/oldParticles==1 || ptcl/oldParticles==3)
        pos(0)=pos(0)+oldBox(0);
      if (ptcl/oldParticles==2 || ptcl/oldParticles==3)
        pos(1)=pos(1)+oldBox(1);
      dVec diff=pos-prev_pos;
      PutInBox(diff,oldBox_dVec);
      pos=diff+prev_pos;
      if (fabs(pos[0]-first_slice(ptcl)[0])>1e-5 || fabs(pos[1]-first_slice(ptcl)[1])>1e-5){
        cerr<<"BROKEN: "<<ptcl<<" "<<pos<<" "<<first_slice(ptcl)<<" "<<prev_pos<<" "<<diff<<" "<<endl;
        cerr<<"BROKEN2: "<<oldPaths(ptcl % oldParticles,0,0)<<" "<<oldPaths(ptcl % oldParticles,0,1)<<" "<<oldPaths(ptcl % oldParticles,0,2)<<endl;
        cerr<<"NUMTIMESLICES "<<NumTimeSlices()<<endl;
        for (int slice=0;slice<TotalNumSlices;slice++)
          cerr<<"A1234 "<<slice<<" "<<oldPaths(ptcl % oldParticles,slice,0)<<" "<<oldPaths(ptcl % oldParticles,slice,1)<<" "<<oldPaths(ptcl % oldParticles,slice,2)<<endl;
      
        assert(pos==first_slice(ptcl));
      }
      Path(NumTimeSlices()-1,ptcl) = pos;
    }
  }
  cerr<<"Almost done initializing"<<endl;
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();
  if (myProc==Communicator.NumProcs()-1)
    MoveJoin(NumTimeSlices()-1,1);
  cerr<<"ending tile " <<Communicator.MyProc()<<endl;
  cerr<<"Done initializing"<<endl;
}


void 
PathClass::ReadOld(string fileName,bool replicate)
{
// This was modified on Jan 19 2005 to read in a set of classical (P=1) configs and duplicate them to produce a set of PIMC (P>1) configs.  -jg
///  cerr<<"Trying to read old"<<endl;
  IOSectionClass inFile;
  assert (inFile.OpenFile(fileName.c_str()));
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,4> oldPaths; //(58,2560,2,3);
  
  assert(inFile.ReadVar("Path",oldPaths));
  perr << "My paths are of size"  << oldPaths.extent(0) << " "
       << oldPaths.extent(1)<<" " << oldPaths.extent(2) << endl;
  
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos; 
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        if (replicate){
          pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,0,dim);
        }
        else{
          pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,slice,dim);
        }
      Path(slice,ptcl) = pos;
    }      
  }
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();
}


/// This function initializes the paths depending on how they are
/// specified to be initialized in the input file.  Currently, the
/// options are :
/// CUBIC:
/// BCC:
/// FIXED:
/// FILE:
/// LEVIFLIGHT
void 
PathClass::InitPaths (IOSectionClass &in)
{
      //cerr<<"Printing Permutation InitPath 0"<<endl;
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
        //cerr<<Permutation(ptcl)<<endl;


  NowOpen=false;
  NowOpen.AcceptCopy();
  //cerr<<"I am in mode "<<GetMode()<<endl;
  //    cerr<<"Printing Permutation InitPath 0p"<<endl;
  //    for (int ptcl=0;ptcl<NumParticles();ptcl++)
        //cerr<<Permutation(ptcl)<<endl;



  SetMode (NEWMODE);
  //  cerr<<"Hello"<<endl;
  //  assert(1==2);
  assert(in.OpenSection ("Particles"));
  //cerr<<"Printing Permutation InitPath 0pp"<<endl;
  //for (int ptcl=0;ptcl<NumParticles();ptcl++)
  //  cerr<<Permutation(ptcl)<<endl;
  
  SetMode(OLDMODE);
  //cerr<<"Printing Permutation InitPath 0ppOLD"<<endl;
  //for (int ptcl=0;ptcl<NumParticles();ptcl++)
  //  cerr<<Permutation(ptcl)<<endl;
  SetMode(NEWMODE);

  int numSpecies = NumSpecies();
//      cerr<<"Printing Permutation InitPath 0ppp"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//      cerr<<Permutation(ptcl)<<endl;

  // Now initialize the Path
//      cerr<<"Printing Permutation InitPath a"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//      cerr<<Permutation(ptcl)<<endl;

  for (int speciesIndex=0; speciesIndex<numSpecies; speciesIndex++)
  {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(in.OpenSection("Species", speciesIndex));
    string InitPaths;
    ///    cerr<<"about to read the string"<<endl;
    in.ReadVar ("InitPaths", InitPaths);
    ///    cerr<<"Read "<<InitPaths<<endl;
    bool replicate = false;
    in.ReadVar ("Replicate", replicate);
//      cerr<<"Printing Permutation InitPath c"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//      cerr<<Permutation(ptcl)<<endl;

    if (InitPaths == "RANDOM") {
      perr << "Don't know how to do RANDOM yet.\n";
      exit(1);
    }
    else if (InitPaths == "LANGEVIN") {
      InitLangevin (in, species);
    }
    else if (InitPaths == "CUBIC") {
      int num = species.NumParticles;
      bool isCubic = 0;
#if NDIM==2   
      isCubic = (Box[0]==Box[1]);
      cerr << isCubic << endl;
#endif
#if NDIM==3    
      isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
      if (!isCubic) {
        perr << isCubic << " A cubic box is current required for cubic initilization\n";
        abort();
      }
      int numPerDim = (int) ceil (pow((double)num, 1.0/3.0)-1.0e-6);
      double delta = Box[0] / numPerDim;
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
        int ix, iy, iz;
        int ip = ptcl-species.FirstPtcl;
        ix = ip/(numPerDim*numPerDim);
        iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
        dVec r;
        r[0] = ix*delta;
        r[1] = iy*delta;
#if NDIM==3    
        iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
        r[2] = iz*delta;
#endif
        for (int slice=0; slice<NumTimeSlices(); slice++) 
          Path(slice,ptcl) = r;
        dVec disp = r;
        cerr << ptcl << " " << r << " " << endl;
      }
    }
    else if(InitPaths=="UniformSphere") {
      dVec r0,r;
      double SphereRadius;
      assert(in.ReadVar("SphereRadius",SphereRadius));
      //      double SphereRadius=31;
      SpeciesClass &species = *SpeciesArray(0);
      double sigma=sqrt(2*species.lambda*tau);
      for (int ptcl =0;ptcl<NumParticles();ptcl++){
        Random.LocalGaussianVec(1,r0);//under common
        ///        cerr<<"in PathClass"<<r0*SphereRadius/sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2])<<endl;
        for (int slice=0;slice <NumTimeSlices();slice++){
          SetPos(slice,ptcl,r0*SphereRadius * (1.0/sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2])));                      
        }
      }
    }
    else if(InitPaths=="UniformPointsSphere") {
      dVec r0,r;
      SpeciesClass &species = *SpeciesArray(0);
      double SphereRadius;
      assert(in.ReadVar("SphereRadius",SphereRadius));
      //      double SphereRadius=31;
      double sigma=sqrt(2*species.lambda*tau);
      for (int ptcl =0;ptcl<NumParticles();ptcl++){
        Random.LocalGaussianVec(1,r0);//under common
        for (int slice=0;slice <NumTimeSlices();slice++){
          Random.LocalGaussianVec(sigma,r);
          r=r0+r;
          SetPos(slice,ptcl,r*SphereRadius  * (1./sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])));
        }
      }
    }
    else if (InitPaths == "BCC") {
      int num = species.NumParticles;
      bool isCubic = 0;
#if NDIM==2   
      isCubic = (Box[0]==Box[1]);
      cerr << isCubic << endl;
#endif
#if NDIM==3
      isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
      if (!isCubic) {
        perr << "A cubic box is current required for cubic initilization\n";
        abort();
      }
      int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/3.0)-1.0e-6);
      double delta = Box[0] / numPerDim;
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
        int ip = (ptcl-species.FirstPtcl)/2;
        int ix, iy, iz;
        ix = ip/(numPerDim*numPerDim);
        iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
#if NDIM==3
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
        //cerr<<"My position is "<<r[0]<<" "<<r[1]<<" "<<r[2]<<endl;
        for (int slice=0; slice<NumTimeSlices(); slice++)
          Path(slice,ptcl) = r;
      }
    }
    else if (InitPaths == "DIAG") {
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
        dVec r;
        r[0] = ptcl;
        r[1] = ptcl;
#if NDIM==3
        r[2] = ptcl;
#endif
        for (int slice=0; slice<NumTimeSlices(); slice++) {
          dVec temp=0.1*r * (1./sqrt(3.0*(double)ptcl*(double)ptcl));
          temp=r+temp;
          Path(slice,ptcl) =temp;
        }
      }
    }
    else if (InitPaths == "BCCe") {
      // BCC lattice with slices forming a circle of radius eps around the lattice site
      int num = species.NumParticles;
      bool isCubic = 0;
#if NDIM==2
      isCubic = (Box[0]==Box[1]);
      cerr << isCubic << endl;
#endif
#if NDIM==3
      isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
#endif
      if (!isCubic) {
        perr << "A cubic box is current required for cubic initilization\n";
        abort();
      }
      int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/3.0)-1.0e-6);
      double delta = Box[0] / numPerDim;
      double eps = 1.0e-2 * delta;
      cout << eps << endl;
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
        int ip = (ptcl-species.FirstPtcl)/2;
        int ix, iy, iz;
        ix = ip/(numPerDim*numPerDim);
        iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
#if NDIM==3
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
        for (int slice=0; slice<NumTimeSlices(); slice++) {
          double dx = eps * cos((slice/(double)NumTimeSlices())*2*M_PI);
          double dy = eps * sin((slice/(double)NumTimeSlices())*2*M_PI);
#if NDIM==3
          double dz = 0.0;
#endif
          dVec dr;
          dr[0] = dx;
          dr[1] = dy;
#if NDIM==3
          dr[2] = dz;
#endif
          Path(slice,ptcl) = r + dr;
        }
      }
    }
    else if (InitPaths=="ALLFIXED"){
      int myFirstSlice, myLastSlice, myProc;
      myProc = Communicator.MyProc();
      SliceRange (myProc, myFirstSlice, myLastSlice);
      Array<double,3> Positions;
      assert(in.ReadVar("Positions",Positions));
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
        for (int slice=0;slice<TotalNumSlices;slice++){

          int sliceOwner = SliceOwner(slice);
          int relSlice = slice-myFirstSlice;
          ///      if (myProc==sliceOwner){
          if (myFirstSlice<=slice && slice<=myLastSlice){
            dVec pos;
            pos = 0.0;
            for (int dim=0; dim<NDIM; dim++)
              pos(dim) = Positions(slice,ptcl,dim);
            Path(relSlice,ptcl) = pos;
          }
        }
        ///If you are the last processors you must make sure the last
        ///slice is the same as the first slice on the first
        ///processors. The  join should be at the
        if (myProc==Communicator.NumProcs()-1){
          dVec pos;
          pos = 0.0;
          for (int dim=0; dim<NDIM; dim++)
            pos(dim) = Positions(0,ptcl,dim);
          Path(NumTimeSlices()-1,ptcl) = pos;
        }
      }
    }
    else if (InitPaths == "FIXED") {
//      cerr<<"Printing Permutation InitPath 1"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//      cerr<<Permutation(ptcl)<<endl;


      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      for (int ptcl=species.FirstPtcl; 
           ptcl<=species.LastPtcl; ptcl++){
        for (int slice=0; slice<NumTimeSlices(); slice++) {
          dVec pos;
          pos = 0.0;
          for (int dim=0; dim<species.NumDim; dim++)
            pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
          Path(slice,ptcl) = pos;
        }      
      }
//      cerr<<"Printing Permutation InitPath 2"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//        cerr<<Permutation(ptcl)<<endl;



//      cerr<<"Printing Permutation InitPath 3"<<endl;
//      for (int ptcl=0;ptcl<NumParticles();ptcl++)
//        cerr<<Permutation(ptcl)<<endl;


      //      PrintRealSlices();
      //      TestRealSlices();
    }
//     else if (InitPaths == "WORM") {
//       NowOpen=true;
//       Array<double,2> Positions;
//       assert (in.ReadVar ("Positions", Positions));
//       assert (Positions.rows() == species.NumParticles);
//       assert (Positions.cols() == species.NumDim);
//       for (int ptcl=species.FirstPtcl; 
//         ptcl<=species.LastPtcl; ptcl++){
//      for (int slice=0; slice<NumTimeSlices(); slice++) {
//        ParticleExist(slice,ptcl)=1.0;
//        dVec pos;
//        pos = 0.0;
//        for (int dim=0; dim<species.NumDim; dim++)
//          pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
//        Path(slice,ptcl) = pos;
//      }      
//       }
//       int startEmpty;
//       assert(in.ReadVar("NumRealParticles",startEmpty));

//       //      int startEmpty=1;
//       for (int ptcl=startEmpty;ptcl<NumParticles();ptcl++)
//      for (int slice=0;slice<NumTimeSlices();slice++){
//        ParticleExist(slice,ptcl)=0.0;
//      }
//       ParticleExist(NumTimeSlices()-1,startEmpty)=0.0;
//       ParticleExist(0,startEmpty)=0.0;
//       //      ParticleExist.AcceptCopy();
//       for (int ptcl=startEmpty-1;ptcl<NumParticles()-1;ptcl++){
//      Permutation(ptcl)=ptcl+1;
//       }
//       Permutation(NumParticles()-1)=startEmpty-1;
//       Permutation.AcceptCopy();
//       Path.AcceptCopy();
//       ParticleExist.AcceptCopy();
//       NowOpen.AcceptCopy();
//     }
    else if (InitPaths == "ADDVACANCIES") {
      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      //      assert (Positions.rows() >= species.NumParticles);
      //      assert (Positions.cols() == species.NumDim);
      ///      cerr<<"My extent 1 is "<<Positions.extent(1);
      ///      cerr<<"MY scale box is "<<ScaleBox<<endl;
      Positions.resizeAndPreserve(species.NumParticles,Positions.extent(1));
      ///      cerr<<"My ptcl are "<<Path.extent(1)<<" and agaisnt "<<Positions.extent(0)<<endl;
      for (int ptcl=species.FirstPtcl; 
           ptcl<=species.LastPtcl; ptcl++){
        for (int slice=0; slice<NumTimeSlices(); slice++) {
          dVec pos;
          pos = 0.0;
          for (int dim=0; dim<species.NumDim; dim++)
            pos(dim) = Positions(ptcl-species.FirstPtcl,dim)*ScaleBox;
          Path(slice,ptcl) = pos;
        }      
      }
      //      cerr<<"My species and numparticles are "<<species.NumParticles<<endl;
      //      Positions.resizeAndPreserve(Positions.extent(0),species.NumParticles);
    }    

    else if (InitPaths == "ALLPATHS") {
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
    else if (InitPaths=="RESTART"){
      string pathFile;
      assert(in.ReadVar("File",pathFile));
      Restart(in,pathFile,false,species);
    }
    else if (InitPaths == "FILE"){
      //cerr<<"I'm going to read the file now"<<endl;
      if (replicate){
        perr << "Replicate 'ON'; Using time slice 0 for all time slices." << endl;
      }
      string pathFile;
      assert(in.ReadVar("File",pathFile));
      ReadSqueeze(in,pathFile,replicate);
      //      ReadOld(pathFile,replicate);
    }
    else if (InitPaths == "TILE"){
      string pathFile;
      assert(in.ReadVar("File",pathFile));
      Tile(in,pathFile,replicate);
    }


    else if (InitPaths == "SQUEEZE"){
      string pathFile;
      assert(in.ReadVar("File",pathFile));
    }
    else if (InitPaths == "LEVIFLIGHT") {
      Array<double,2> Positions;
      double sigmaFactor = 1.0;
      assert (in.ReadVar ("Positions", Positions));
      if (in.ReadVar ("SigmaFactor", sigmaFactor, 1.0))
        perr << "SigmaFactor = " << sigmaFactor << endl;
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      Array<dVec,1> R0(species.NumParticles);
      for (int ptcl=0; ptcl<species.NumParticles; ptcl++) 
        for (int dim=0; dim<NDIM; dim++)
          R0(ptcl)[dim] = Positions(ptcl,dim);
      // NodeAvoidingLeviFlight (speciesIndex,R0);
      PhaseAvoidingLeviFlight(speciesIndex, R0, sigmaFactor);
    }
    else if (InitPaths == "RANDOMFIXED") {
      cout <<CloneStr<< " Random fixed initialization." << endl;
      InitRandomFixed (in, species);
    }
//     else if (InitPaths == "LEVIFLIGHT") {
//       int myStart, myEnd;
//       SliceRange (Communicator.MyProc(), myStart, myEnd);
//       Array<double,2> Positions;
//       assert (in.ReadVar ("Positions", Positions));
//       assert (Positions.rows() == species.NumParticles);
//       assert (Positions.cols() == species.NumDim);
//       Array<dVec,1> flight(TotalNumSlices+1);
//       for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
//      for (int dim=0; dim<NDIM; dim++) {
//        flight(0)[dim] = Positions(ptcl-species.FirstPtcl,dim);
//        flight(TotalNumSlices)[dim] = Positions(ptcl-species.FirstPtcl,dim);
//      }
//      LeviFlight (flight, species.lambda);
// //   for (int i=0; i<MyNumSlices; i++)
// //     Path(i, ptcl) = flight(i-RefSlice);
//      for (int slice=myStart; slice<=myEnd; slice++)
//        Path(slice-myStart, ptcl) = flight(slice);
//       }
//     }

    else {
      perr << "Unrecognize initialization strategy " << InitPaths << endl;
      abort();
    }
    //cerr<<Communicator.MyProc()<<" "<<species.Name<<" Init real slices"<<endl;
    InitRealSlices();
    in.CloseSection(); // Species
  }
  in.CloseSection(); // "Particles"

  string openSpeciesName;
  if (OpenPaths) {
    assert(in.ReadVar("OpenSpecies",openSpeciesName));
    OpenSpeciesNum=SpeciesNum(openSpeciesName);
    InitOpenPaths();
  }
  if (OrderN){
    cerr<<"inside order n things"<<endl;
    double cutoff;
    assert(in.ReadVar("CutoffDistance",cutoff));
    Cell.CutoffDistance=cutoff;
    //    Array<int,1>  numGrid(NDIM);
    Array<int,1>  numGrid;
    assert(in.ReadVar("CellGridSize",numGrid));
    assert(numGrid.size()==NDIM);
    //  numGrid(0)=6*2;
    //  numGrid(1)=6*2;
    //  numGrid(2)=5*2;
    ///    numGrid(0)=2;
    ///    numGrid(1)=3;
    ///    numGrid(2)=20;
    //    //    numGrid(0)=10;
    //    //    numGrid(1)=10;
    //    numGrid(2)=10;
    numGrid(0)=25;
    numGrid(1)=25;
#if NDIM==3
    numGrid(2)=30;
#endif

    //    numGrid(0)=30;
    //    numGrid(1)=30;
    //#if NDIM==3
    //    numGrid(2)=30;
    //#endif


    //  Cell=new GridClass(*this);
    Cell.Init(Box,numGrid);
    //  Cell.BuildNeighborGrids();
    //  cerr<<"I am now printing neighbor grids"<<endl;
    //  Cell.PrintNeighborGrids();
    for (int slice=0;slice<NumTimeSlices();slice++)
      Cell.BinParticles(slice);
    //  cerr<<"I have binned them"<<endl;
    //Cell.PrintParticles(0);
  }
  //Everything needs to be accepted
  Path.AcceptCopy();
  Permutation.AcceptCopy();
  Rho_k.AcceptCopy();
  BroadcastRefPath();
  RefPath.AcceptCopy();
  SignWeight.AcceptCopy();
  ExistsCoupling.AcceptCopy();
  NodeDist.AcceptCopy();
  NodeDet.AcceptCopy();

}



void
PathClass::InitRandomFixed(IOSectionClass &in,
                           SpeciesClass &species)
{
  double radius;
  assert (in.ReadVar("Radius", radius));
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
      cerr << "Exceed 1000 tries in InitRandomFixed.  "
           << "Decrease excluded radius.\n";
      abort();
    }
  }
  for (int slice=0; slice<NumTimeSlices(); slice++)
    for (int ptcl=0; ptcl<N; ptcl++) {
      Path[0](slice, ptcl + species.FirstPtcl) = R(ptcl);
      Path[1](slice, ptcl + species.FirstPtcl) = R(ptcl);
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

#include <unistd.h>

void 
PathClass::NodeAvoidingLeviFlight (int speciesNum, Array<dVec,1> &R0)
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
      if (! Actions.NodalActions(speciesNum)->IsPositive(0)) {
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
  
  int N = TotalNumSlices+1;
  for (int slice=1; slice<N; slice++) {
    int sliceOwner = SliceOwner(slice);
    int relSlice = slice-myFirstSlice;
     
    double delta = (double)(N-slice-1);
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
        if (sliceOwner == myProc ) {
          for (int ptcl=0; ptcl<numPtcls; ptcl++)
            (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
          positive = 
            Actions.NodalActions(speciesNum)->IsPositive(relSlice);
        }
        // Now broadcast whether or not I'm positive to everyone
        Communicator.Broadcast(sliceOwner, positive);
      }
      if (!positive) {
        numRejects++;
        if ((numRejects%50)==0)
          cerr << "numRejects = " << numRejects << endl;
      }
    } while (!positive);
    // Copy slice into Path if I'm the slice owner.
    if ((slice>=myFirstSlice) && (slice<=myLastSlice)) {
      for (int ptcl=0; ptcl<numPtcls; ptcl++) 
        (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
      // Check to make sure we're positive now.
      if (haveNodeAction) {
//      if ((slice==myFirstSlice)||(slice==myLastSlice)) {
//        fprintf (stderr, "myProc=%d slice=%d det=%1.8e\n", myProc, slice, 
//                 Actions.NodalActions(speciesNum)->Det(relSlice));
//        Array<double,2> matrix(numPtcls,numPtcls);
//        matrix = Actions.NodalActions(speciesNum)->GetMatrix(relSlice);
//        char fname[100];
//        snprintf (fname, 100, "matrix%d-%d.dat", slice, myProc);
//        FILE *fout=fopen (fname, "w");
//        for (int i=0; i<numPtcls; i++) {
//          for (int j=0; j<numPtcls; j++) 
//            fprintf (fout, "%1.16e ", matrix(i,j));
//          fprintf (fout, "\n");
//        }
//        fclose(fout);
//        for (int ptcl=0; ptcl<numPtcls; ptcl++) 
//          fprintf (stderr, "myProc=%d newSlice(%d)=[%10.7e %10.7e %10.7e]\n",
//                   myProc, ptcl, 
//                   newSlice(ptcl)[0], newSlice(ptcl)[1], newSlice(ptcl)[2]);
//      }
        if (!Actions.NodalActions(speciesNum)->IsPositive(relSlice)) {
          perr << "Still not postive at slice " << slice 
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
    Communicator.PrintSync();
//     perr << "myProc = " << myProc << " localAction = " 
//       << localAction << " NumTimeSlices = " << NumTimeSlices() << endl;
    double globalAction = Communicator.AllSum(localAction);
    cerr << "After AllSum.\n";
    if (Communicator.MyProc()==0)
      perr << "Nodal Action after Levi flight = " << globalAction << endl;
  }
  Actions.NodalActions(speciesNum)->AcceptCopy(0, NumTimeSlices()-1);


//   Communicator.PrintSync();
//   char fname[100];
//   snprintf (fname, 100, "%s.dat", species.Name.c_str());
//   FILE *fout;
//   if (myProc == 0)
//     fout = fopen (fname, "w");
//   else
//     fout = fopen (fname, "a");
//   for (int slice=0; slice<(NumTimeSlices()-1); slice++) {
//     for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) 
//       for (int i=0; i<NDIM; i++)
//      fprintf (fout, "%1.12e ", (*this)(slice, ptcl)[i]);
//     fprintf (fout, "\n");
//   }
//   fclose(fout);
}


void 
PathClass::PhaseAvoidingLeviFlight (int speciesNum, Array<dVec,1> &R0,
                                    double sigmaFactor)
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


void
PathClass::InitLangevin(IOSectionClass &in, 
                        SpeciesClass &species)
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
