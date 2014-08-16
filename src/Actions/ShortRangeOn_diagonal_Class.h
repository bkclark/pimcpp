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

#ifndef SHORT_RANGE_ON_DIAGONAL_CLASS_H
#define SHORT_RANGE_ON_DIAGONAL_CLASS_H

#include "ActionBase.h"
#include <vector>
/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class ShortRangeOn_diagonal_class : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<bool,1> DoPtcl;
  Array<bool,1> todoIt;
  Array<int,1> ptcls;
  Array<dVec,1> gradVec;
  Array<double,1> gradVecSquared;
  vector<double> distancesNew;
  vector<double> distancesOld;
  vector<double> factorArrayNew;
  vector<double> factorArrayOld;
  vector<double> *distances;
  vector<double> *factorArray;
public:
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double SingleAction_slow (int slice1, int slice2,
			    const Array<int,1> &changedParticles,
			    int level);
  double SingleAction_fast (int slice1, int slice2,
			    const Array<int,1> &changedParticles,
			    int level);


  double d_dBeta (int slice1, int slice2, int level);
  double d_dBeta_check (int slice1, int slice2, int level);
  double d_dBeta_slow (int slice1, int slice2, int level);
  string GetName();
  double residual_energy();
  void GradAction_help(int slice1, int slice2,
		       const Array<int,1> &ptcls, 
		       int level,
		       Array<dVec,1> &gradVec,
		       Array<double,1> &gradVecSquared);
  ShortRangeOn_diagonal_class (PathDataClass &pathData,
			       Array<PairActionFitClass*, 2> &pairMatrix);
};

#endif
