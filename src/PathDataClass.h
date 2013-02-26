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

#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Communication/Communication.h"
#include "Common.h"
#include "SpeciesClass.h"
#include "PathClass.h"
#include "Actions/ActionsClass.h"

#include "Random/Random.h"
#include "MoleculeHelper.h"

#ifdef USE_QMC
  #include <QMCApp/QMCInterface.h>
  #include "Message/Communicate.h"
  #include "Utilities/OhmmsInfo.h"
#endif



/// This is the class that holds all of the information about the paths 
/// including all of the particles, the momized data, and the action.
/// Has routines for accepting and rejecting and for shifting data
/// between the processors.
class PathDataClass
{
private:
  int MyCloneNum, NumClones;
  /// These store the time of start and the maximum wall time in seconds.
  int StartWallTime, MaxWallTime;
  int GetWallTime();

public:
  MoleculeManagerClass Mol;

  bool RUN_QMC;
  bool IAmQMCManager;
  int Seed;
  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  /// This is for commmunication between nodes within a clone group.
  CommunicatorClass IntraComm;
  /// This is for communication between the rank 0 nodes of each clone
  /// group.  Hence, its between the clone groups.
  CommunicatorClass InterComm;
  /// This is the global MPIWORLD communicator.
  CommunicatorClass WorldComm;
  RandomClass Random;

#ifdef USE_QMC
  qmcplusplus::QMCInterface* qmc;
  CommunicatorClass MetaWorldComm;
  CommunicatorClass QMCComm;
  void AssignPtclSetStrings();
  bool useDefaultStrings;
  // user-defined parameters to control QMC run
  // read in from input file
public:
  Array<string, 1> ptclSet0, ptclSet1;
  double dt;
  int walkers, chains, steps, blocks;
  bool correlated;
  string QMCMethod;
#endif

  // global clock for number of moves made
  // incremented (by ALL moves) in MoveBase::DoEvent()
  int moveClock;

  ///////////////////////////////////////////////////////////////////
  ///                        Wall Time Data                        //
  ///////////////////////////////////////////////////////////////////
  /// This function returns true if we have exceeded the maximum wall
  /// time.
  void SetMaxWallTime(int maxWallTime);
  bool ExceededWallTime();
  /// This object computes all actions.
  //  ActionClass Action; //(MemoizedDataClass,SpeciesArrayClass);
  ActionsClass Actions;
  /// The constructor that initializes the action and such
  int &Join;
  PathClass Path;
  inline void ShiftData(int numTimeSlicesToShift){
    Path.ShiftData(numTimeSlicesToShift);
    Actions.ShiftData(numTimeSlicesToShift);
  }

  /// We are probaby going to have to move permutation
  /// information up here if we want it to notice
  /// the change in join.
  inline void MoveJoin(int newJoin)
  {
    Path.MoveJoin(Join,newJoin);
    Actions.MoveJoin(Join,newJoin);
    Join=newJoin;
  }

  // Calculates the centroid postions of all the paths
  void GetCentroids(Array<TinyVector<double,NDIM>,1> &CentPos, Array<int,1> &activeParticles);

  ////Worm Moves///
  bool SliceFullandNextSliceEmpty(int slice,int ptcl);
  bool SliceFullandPreviousSliceEmpty(int slice,int ptcl);
  void FindHead(int &headSlice,int &headPtcl);
  void FindTail(int &tailSlice,int &tailPtcl);
  void MoveTailToSliceZero();
  void Next(int &slice, int &ptcl);
  void WormInfo(int &headSlice, int &headPtcl, int &tailSlice, int &tailPtcl, int &numEmpty, int &wormSize);

  //////////
  void MoveOpenLinkToEnd();
  void MoveLinkToEnd(int linkToMove);
  /// This function shifts the path to make the reference slice
  /// equal to the absolute slice position absSlice.  The join is left
  /// on the last slice, so the permutation is out of the way.
  void MoveRefSlice (int absSlice);

  /// Returns the number of time slices.
  inline int NumTimeSlices() {  return Path.NumTimeSlices();  }

  /// Do all copies necessary to accept a move.
  inline void AcceptMove(int startTimeSlice,int endTimeSlice, const Array <int,1> &activeParticles);

  /// Do all copies necessary to accept a move.
  inline void RejectMove(int startTimeSlice,int endTimeSlice, const Array <int,1> &activeParticles);

  /// Returns the number of particle species in the path
  inline int NumSpecies(){ return Path.NumSpecies();  }

  /// Returns the total number of particles
  inline int NumParticles(){ return Path.NumParticles(); }

  /// Returns a reference to the SpeciesClass object of number species
  inline SpeciesClass& Species(int species){ return Path.Species(species); }

  inline int SpeciesNum(string speciesName);

  /// Returns the position of the particle of type species, particle
  /// number particle, and time slice timeSlice.
  inline dVec operator()(int timeSlice,int particle){
    return Path(timeSlice,particle);
  }

  /// Sets the position of the particle labeled by
  /// (species, particle, timeSlice) to r and updates the time stamp
  /// of that piece of information.
  inline void SetPos(int timeSlice, int particle, const dVec& r){
    //    Path.SetPos(timeSlice,particle,r);
    Path(timeSlice,particle) = r;
  }

  inline int GetCloneNum()  { return MyCloneNum; }
  inline int GetNumClones() { return NumClones; }

  void Read (IOSectionClass &in);

  /// Constructor
  PathDataClass() :
    /* Action(*this), */Actions(*this), Random(WorldComm), Path(IntraComm,Random, Actions), MaxWallTime(-1), Join(Path.Join)
  {
    Join = 1;
    StartWallTime = GetWallTime();
  }


};

inline int PathDataClass::SpeciesNum(string name)
{
  for (int spec=0;spec<NumSpecies();spec++){
    if (Species(spec).Name==name){
      return spec;
    }
  }
  return -1;
}

inline void PathDataClass::AcceptMove(int startTimeSlice, int endTimeSlice, const Array <int,1> &activeParticles)
{
  Actions.AcceptCopy(startTimeSlice, endTimeSlice, activeParticles);
  Path.AcceptCopy(startTimeSlice, endTimeSlice, activeParticles);
}

inline void PathDataClass::RejectMove(int startTimeSlice, int endTimeSlice, const Array <int,1> &activeParticles)
{
  Actions.RejectCopy(startTimeSlice, endTimeSlice, activeParticles);
  Path.RejectCopy(startTimeSlice, endTimeSlice, activeParticles);
}


#endif
