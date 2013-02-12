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

#ifndef MOVE_BASE_H
#define MOVE_BASE_H

#include "../PathDataClass.h"
#include "../Observables/ObservableVar.h"
#include "../EventClass.h"

/// This is the generic parent class for all moves, including "real moves"
/// which actually move particles and "pseudo moves", which just shift around
/// data, but don't move anything physical.
class MoveClass : public EventClass
{
protected:
  double SecondsInMove;
  /// This variable stores the acceptance ratio for the move
  ObservableDouble RatioVar;
  int DumpFreq;
public:
  /// Call this in order to make a move.
  virtual void MakeMove()=0;
  ///All moves ought to be able to read
  virtual void Read(IOSectionClass &input)=0;
  virtual double AcceptanceRatio() {return sqrt((double)-1.0);}
  virtual void WriteRatio();
  void DoEvent();
  
  /// MoveClass constructor. Sets reference to the PathData object
  MoveClass(PathDataClass &pathData, IOSectionClass &out) : 
    EventClass(pathData, out), DumpFreq(10000), 
    RatioVar("AcceptRatio", IOSection, pathData.Path.Communicator)
    {
      ///		cerr << "MoveClass construct" << endl;
      // do nothing 
    }
};
  

/// This is a specialization of MoveClass which actually physically moves
/// particles.
class ParticleMoveClass : public MoveClass
{
protected:
  /// Stores which species of particles this moves will work on
  Array<int,1> ActiveSpecies;
  /// Total number of particles in the active species
  int TotalParticles;

  /// Scratch Array holding a random subset of particles
  Array<int,1> MyParticleIndices; 
  /// Our move class takes the number of particles to move at a
  /// time. This is stored here. 
  int NumParticlesToMove;

 public:
  /// Stores the number of moves made and the number accepted
  int NumAccepted;
  /// Desired acceptance ratio
  double DesiredAcceptRatio;
  /// An accumulator used to publish the diffusion value. -jg
  double total_r_squared;
  /// This returns the Acceptance Ratio.
  inline double AcceptanceRatio() 
  {
    return (double)(NumAccepted)/(double)TimesCalled;
  }
  /// Call this to make a move
  //  virtual void MakeMove()=0;

  /// This array contains the int's of particles that you are 
  /// currently moving (i.e. NumParticlesToMove of them
  Array<int,1> ActiveParticles;
  /// When we choose particles we select the  particles (randomly)
  /// from the set of species enumerated in ActiveSpecies 
  void SetActiveSpecies(Array<int,1> ActSpecies);
  /// Function that sets the number of particles to move
  inline void SetNumParticlesToMove(int i)
  {
    NumParticlesToMove = i;
    MyParticleIndices.resize(i);
    ActiveParticles.resize(i);
  }
  /// Function that chooses the particles that you should move and
  /// places them in ActiveParticles; 
  void ChooseParticles();
  void ChooseParticlesOpen();
  ParticleMoveClass(PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass (myPathData, outSection)
  { 
    NumAccepted=0;
    /* Do nothing for now.*/  
  }
};


#endif
