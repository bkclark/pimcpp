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

#ifndef CUMMINGS_WATER_POTENTIAL_CLASS_H
#define CUMMINGS_WATER_POTENTIAL_CLASS_H

#include "ActionBase.h"
#include "../Moves/MoveUtils.h"

class CummingsWaterPotentialClass : public ActionBaseClass
{


	bool ReadComplete;
	bool firstTime, with_truncations;
	enum Interactions{LJ,CHARGE};
	Array<bool, 2> Interacting;
	Array<string, 1> LJSpecies;
	Array<string, 1> ChargeSpecies;
  Array<double,1> k2;
  Array<double,3> phi;
	double Dyn2kcal;

  double volfactor;
  // model parameters
  double sigma_COM, alpha, tolerance, gamma, mu_0;
  // reaction field parameters
  double Rcut_RF, dielectric_RF, prefactor_RF;

  // a 3d array of 2d (3x3) matrices!
  Array<Array<double,2>, 3> T;

  // containers for vectors
  Array<dVec,2> Eq;
  Array<dVec,2> Ep;
  Array<dVec,2> p, p_permanent;
  Array<dVec,2> RF_q, RF_p;

public:
 	// Parameters that can be specified on input
	// Defaults in constructor
	double prefactor, conversion, prefactor_Efield;
	double CUTOFF;

	// hardwired units and unit conversions
	// set in constructor
  double elementary_charge, N_Avogadro, kcal_to_joule,
		epsilon_not, angstrom_to_m, SI, k_B, erg_to_eV, joule_to_eV, bohr_per_angstrom;

	double sigma_over_cutoff, offset;

  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
	double ComputeEnergy(int slice1, int slice2,
						const Array<int,1> &activeParticles, int level, bool isAction);
	void SetNumImages(int num);
  string GetName();
 
  void Update_T(int slice);
  void Update_Ep(int slice);
  void Update_Eq(int slice);
  void Update_Eq(int slice, int i);
  void Update_RF_p(int slice);
  void Update_RF_q(int slice);
  void SolveDipole(int slice);
	CummingsWaterPotentialClass (PathDataClass &pathData); };

#endif
