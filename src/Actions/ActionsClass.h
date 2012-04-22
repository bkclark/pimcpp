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

#ifndef ACTIONS_CLASS_H
#define ACTIONS_CLASS_H


#include "DiagonalActionClass.h"
#include "ShortRangeClass.h"
#include "ShortRangeOnClass.h"
#include "ShortRangeOn_diagonal_Class.h"
#include "ShortRangeApproximateClass.h"
#include "ShortRangePrimitive.h"
#include "LongRangeClass.h"
//#include "LongRangeCoulombClass.h"
#include "LongRangeRPAClass.h"
#include "ShortRangePotClass.h"
#include "LongRangePotClass.h"
#include "KineticClass.h"
#include "WaterClass.h"
//#include "MoleculeInteractionsClass.h"
//#include "ST2WaterClass.h"
//#include "TIP5PWaterClass.h"
//#include "EAMClass.h"
#include "NodalActionClass.h"
#include "FreeNodalActionClass.h"

// #include "GroundStateNodalActionClass.h"

#include "DavidLongRangeClass.h"
#include "DavidLongRangeClassYk.h"
//#include "QMCSamplingClass.h"
//#include "QBoxAction.h"
#include "OpenLoopImportance.h"
#include "StructureReject.h"
//#include "KineticRotorClass.h"
// #include "KineticSphereClass.h"
#include "TruncatedInverse.h"
#include "Mu.h"
#include "VariationalPI.h"
#include "Tether.h"
//#include "NonlocalClass.h"


/// ActionsClass is a shell of a class holding all of the necessary
/// ActionBaseClass derivatives representing the different actions.
/// It includes kinetic, short range, long range, long range RPA
/// version, and nodal actions.  It may later contain importance
/// sampling "actions" as well.
class ActionsClass
{
private:

  // new functionality
  // store pointers to ONLY action classes that are used (multiple instantiation also is allowed)
  void ReadPairActions(IOSectionClass &in);
  public:
  ActionBaseClass* GetAction(string name);
  Array<PairActionFitClass*,1> PairArray;
  private:
  // This stores pointers to the pair actions.


  PathDataClass &PathData;
  int MaxLevels; //is this the right place for this?
  void ReadNodalActions (IOSectionClass &in);
  //  FixedPhaseClass *FixedPhaseA, *FixedPhaseB;
public:
  std::list<ActionBaseClass*> ActionList;
  std::list<string> ActionLabels;
  std::list<ActionBaseClass*> EnergyList;

  Array<double,1> TauValues;
  Array<PairActionFitClass*,1> SpecificHeatPairArray;
  // This stores pointers to pair action fits for each pair of species.
  Array<PairActionFitClass*,2> PairMatrix;
  VariationalPIClass VariationalPI;
  TruncatedInverseClass TruncatedInverse;
  /// Used to keep track of the total action
  double TotalA, TotalB;

  /// Specifies whether to use long range breakups
  bool UseLongRange;
  /// Specifies whether we have a nonlocal potential
  bool UseNonlocal;


  // Actions
  OpenLoopImportanceClass OpenLoopImportance;
  StructureRejectClass StructureReject;
  /// The Kinetic action
  KineticClass Kinetic;
  //  KineticSphereClass KineticSphere;
  //KineticRotorClass KineticRotor;
  //FixedAxisRotorClass FixedAxisRotor;
  MuClass Mu;
  /// The short range part of the pair action.  This is the complete
  /// pair action in the case of short-range potententials.  The
  /// short-range action is summed in real space. 
  ShortRangeClass ShortRange;
  ShortRangeOnClass ShortRangeOn;
  ShortRangeOn_diagonal_class ShortRangeOnDiagonal;
  ShortRangeApproximateClass ShortRangeApproximate;
  ShortRangePrimitiveClass ShortRangePrimitive;
  DiagonalActionClass DiagonalAction;

  /// The long range part of the action, which is summed in k-space.  
  LongRangeClass LongRange;

  //LongRangeCoulombClass LongRangeCoulomb;
  
  /// The Random Phase Approximation-corrected form of the above.
  LongRangeRPAClass LongRangeRPA;

  ///David's Long Range Class
  DavidLongRangeClassYk DavidLongRange;

//  // Water-related stuff
//  MoleculeInteractionsClass MoleculeInteractions;
//  QBoxActionClass QBoxAction;
//  ST2WaterClass ST2Water;
//  TIP5PWaterClass TIP5PWater;

  // Nonlocal action
  //  NonlocalClass Nonlocal;

//  EAMPotentialClass EAM;
//
//#ifdef USE_QMC
//  CEIMCActionClass CEIMCAction;
//#endif
//
//  QMCSamplingClass QMCSampling;
//
//	IonIonActionClass IonInteraction;

  /// This array of actions are used for Restricted PIMC for
  /// fermions.  These effective actions ensure that the paths do not
  /// cross the nodes of some trial density matrix with respective to
  /// the reference slice.
  Array<NodalActionClass *,1> NodalActions;
  //DiagonalClass Diagonal;
  //ImportanceSampleClass ImportanceSample;

  // Potentials
  ShortRangePotClass ShortRangePot;
  LongRangePotClass  LongRangePot;
  TetherClass Tether;
  /// Stores whether we use Random Phase Approximation corrections to
  /// the long range action.
  bool UseRPA;

  /// Stores number of images to sum over for kinetic action and energy.
  int NumImages;

  Potential &GetPotential (int species1, int species2);
  /// Return the all the actions for this processor's segment of
  /// the path.  Must do global sum to get total action.
  void GetActions (double& kinetic, double &duShort, double &duLong, 
		   double &node);

  /// This function adds to F the current forces calculated from the
  /// gradient of the action.  Note that you must do an AllSum over
  /// the clone processors to get the total.
  void GetForces(const Array<int,1> &ptcls, 
		 Array<dVec,1> &Fshort, Array<dVec,1> &Flong);
  /// Finite difference version for testing.
  void GetForcesFD(const Array<int,1> &ptcls, Array<dVec,1> &F);

  /// Return the all the energies for this processor's segment of
  /// the path.  Must do global sum to get total energy.
	//void Energy(map<double>& Energies);
  void Energy (double& kinetic, double &duShort, double &duLong, 
	       double &node, double &vShort, double &vLong,
	       double &duNonlocal);
  void Energy (double& kinetic, double &duShort, double &duLong, 
	       double &node, double &vShort, double &vLong,
	       double &duNonlocal,double &residual);

  /// Read the action parameters from the input file and do the
  /// necessary initialization.  This reads the pair actions, and does
  /// long range breakups and RPA corrections if necessary.
  void Read(IOSectionClass &in);

  /// This routine does any necessary shifting of stored data for the
  /// paths.  It just calls the ShiftData
  void ShiftData (int slicesToShift);
  void MoveJoin (int oldJoinPos, int newJoinPos);

  void AcceptCopy (int startSlice, int endSlice,
		   const Array<int,1> &activeParticles);
  void RejectCopy (int startSlice, int endSlice,
		   const Array<int,1> &activeParticles);
  /// This should be called after the paths have been constructed to
  /// initialize any cached data;
  void Init();

  void Setk(Vec3 k);

  bool HaveLongRange();

  /// This function writes any pertinent information related to the
  /// actions to the output file.
  void WriteInfo (IOSectionClass &out);

  void UpdateNodalActions();

  inline int GetMaxLevels() { return MaxLevels; }

  ActionsClass(PathDataClass &pathData) : 
    Tether(pathData),
    ShortRange(pathData,PairMatrix),
    ShortRangeOn(pathData,PairMatrix),
    ShortRangeOnDiagonal(pathData,PairMatrix),
    ShortRangeApproximate(pathData,PairMatrix),
    ShortRangePrimitive(pathData,PairMatrix),
    ShortRangePot(pathData, PairMatrix),
    DiagonalAction(pathData,PairMatrix),
    LongRange(pathData,PairMatrix,PairArray), 
    //LongRangeCoulomb(pathData,PairMatrix,PairArray), 
    DavidLongRange(pathData,PairMatrix,PairArray),
    LongRangeRPA(pathData, PairMatrix, PairArray),
    LongRangePot(pathData, PairMatrix),
    OpenLoopImportance(pathData),
    Kinetic(pathData),
    //KineticRotor(pathData),
    //FixedAxisRotor(pathData),
      //    KineticSphere(pathData),
    PathData(pathData),

//    MoleculeInteractions(pathData),
//    QBoxAction(pathData),
//    ST2Water(pathData),
//    TIP5PWater(pathData),
//    EAM(pathData),
//#ifdef USE_QMC
//    CEIMCAction(pathData),
//#endif
//    QMCSampling(pathData),
//    IonInteraction(pathData),
    Mu(pathData),
    VariationalPI(pathData),
    StructureReject(pathData),
    TruncatedInverse(pathData),
    NumImages(1),
    UseLongRange(true)
      //    Nonlocal (pathData),
      //    FixedPhaseA(NULL), FixedPhaseB(NULL)
  {
    ///Do nothing for now
  }





};

#endif
