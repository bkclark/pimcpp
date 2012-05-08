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

#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../Common.h"
#include "../IO/IO.h"


using namespace IO;

class PathDataClass;
  class PathClass;

/// The ActionBaseClass is an abstract base class from which all the
/// physical action classes are derived.  Its primary methods are
/// Action and its beta-derivative.
class ActionBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  double TimeSpent;
  /// This function takes a range of time slices from slice1 to
  /// slice2, inclusive, and an array of particles which are changing
  /// positions.  The level is used in bisection moves in which, at
  /// higher levels, we skip 2^level intervening time slices, building
  /// up a new path in a recursive style.
  virtual double SingleAction(int slice1, int slice2,
			      const Array<int,1> &activeParticles,
			      int level) = 0;

  /// This is the main function that gets called.  It simply returns
  /// SingleAction if we are not using correlated sampling.  If we are
  /// using correlated sampling, it returns the average of the actions
  /// for the A and B configurations.  At level 0, we must do
  /// something special, so Action should not be called with
  /// Correlated sampling on and level = 0;  Note that this is
  /// virtual, so that the nodal classes can override.
  virtual double Action(int slice1, int slice2,
			const Array<int,1> &activeParticles,
			int level);

  /// This function returns the \f$\beta\f$-derivative of the above
  /// function.  Since we are interested in total energy, it does not
  /// take a list of particles we are moving. 
  virtual double d_dBeta (int slice1, int slice2,
			  int level) = 0;

  /// This function returns the gradient of the action on a specified
  /// set of particles, summed over the timeslices from slice1 to
  /// slice2, inclusive.  The function adds its contribution to
  /// whatever is already in gradVec.  The default implementation does
  /// nothing.  
  virtual void GradAction (int slice1, int slice2, 
			   const Array<int,1> &ptcls, int level,
			   Array<dVec,1> &gradVec);
			     
  // This returns the sum over all time slices, using MPI
  // to sum over processors if necessary.
  virtual void Read (IOSectionClass &in);
  virtual void ShiftData (int slices2Shift);
  virtual void MoveJoin (int oldJoinPos, int newJoinPos);
  virtual string GetName()=0;
  virtual void AcceptCopy (int slice1, int slice2);
  virtual void RejectCopy (int slice1, int slice2);
  ActionBaseClass(PathDataClass &pathData);				   
};

/// The PotentialBaseClass is an abstract base class for storing
/// potentials that go into the energy.  Currently, the ShortRange and
/// LongRange potentials are derived from it.
class PotentialBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  virtual double V (int slice) = 0;
  PotentialBaseClass (PathDataClass &pathData);
};

#endif
