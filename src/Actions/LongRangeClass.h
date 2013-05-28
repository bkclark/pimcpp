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

#ifndef LONG_RANGE_CLASS_H
#define LONG_RANGE_CLASS_H

#include "ActionBase.h"
#include "../PairAction/PAFit.h"
#include "../Ewald/OptimizedBreakup.h"
#include "../Integration/GKIntegration.h"

typedef enum {JOB_U, JOB_DU, JOB_V} JobType; 

/// The LongRangeClass is an action class responsible for the long
/// wavelength components of the action that are summed in k-space.
/// This class has member functions which perform the optimized
/// short range/long range breakup of the action and its
/// beta-derivative.  This is based on modified version of the method
/// by Natoli and Ceperley (J. Comp. Phys. 117, 171-178 [1995]).
class LongRangeClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<PairActionFitClass*,1> &PairArray;
  LinearGrid LongGrid;
  void OptimizedBreakup_U(IOSectionClass &out);
  void OptimizedBreakup_dU(IOSectionClass &out);
  void OptimizedBreakup_V(IOSectionClass &out);


  /// This calculates the quantity 
  /// \f$ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
  double CalcXk (int paIndex, int level, double k, double rc, JobType type);
  // This must be called after all of the OptimizedBreakup_x's

  int Level, ki;
public:
  int numKnots;
  bool UseBackground;
  void WriteInfo(IOSectionClass &out);
  void Init(IOSectionClass &out);
  void Read(IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2,  int level);
  void GradAction (int slice1, int slice2, const Array<int,1> &ptcls,
		   int level, Array<dVec,1> &gradVec);
  string GetName();
  LongRangeClass (PathDataClass &pathData,
		  Array<PairActionFitClass*, 2> &pairMatrix,
		  Array<PairActionFitClass*, 1> &pairArray);

};

#endif
