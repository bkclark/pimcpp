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

#ifndef LONG_RANGE_COULOMB_CLASS_H
#define LONG_RANGE_COULOMB_CLASS_H

#include "ActionBase.h"
#include "../PairAction/PAFit.h"

//typedef enum {JOB_U, JOB_DU, JOB_V} JobType; 

/// The LongRangeClass is an action class responsible for the long
/// wavelength components of the action that are summed in k-space.
/// This class has member functions which perform the optimized
/// short range/long range breakup of the action and its
/// beta-derivative.  This is based on modified version of the method
/// by Natoli and Ceperley (J. Comp. Phys. 117, 171-178 [1995]).
class LongRangeCoulombClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<PairActionFitClass*,1> &PairArray;
  //LinearGrid LongGrid;
  //void OptimizedBreakup_U(int numKnots,  IOSectionClass &out);
  //void OptimizedBreakup_dU(int numKnots, IOSectionClass &out);
  //void OptimizedBreakup_V(int numKnots,  IOSectionClass &out);

  // jg i need these!
	Array<bool,1> activeSpecies;
  double alpha, sq_alpha;
  Array<double,1> k2;
  Array<double,1> phi;
  Array<double,1> PtclCharge;
  double self;
  double volume;
  double prefactor;
  double elementary_charge, N_Avogadro, kcal_to_joule,
		epsilon_not, angstrom_to_m, SI, k_B, erg_to_eV, joule_to_eV, bohr_per_angstrom;

  bool initPhi, LRfirstTime;
  bool doRealSpace;

  /// This calculates the quantity 
  /// \f$ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
  //double CalcXk (int paIndex, int level, double k, double rc, JobType type);
  // This must be called after all of the OptimizedBreakup_x's

  int Level, ki;
public:
  //bool UseBackground;
  //void Init(IOSectionClass &in, IOSectionClass &out);
  void Read (IOSectionClass &in);
  double ComputeEnergy (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2,  int level);
  void GradAction (int slice1, int slice2, const Array<int,1> &ptcls,
		   int level, Array<dVec,1> &gradVec);
  string GetName();
  LongRangeCoulombClass (PathDataClass &pathData,
		  Array<PairActionFitClass*, 2> &pairMatrix,
		  Array<PairActionFitClass*, 1> &pairArray);

};

#endif
