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

#ifndef MOLECULE_INTERACTIONS_CLASS_H
#define MOLECULE_INTERACTIONS_CLASS_H

#include "ActionBase.h"
#include "../Moves/MoveUtils.h"

class MoleculeInteractionsClass : public ActionBaseClass
{
	bool ReadComplete;
	bool withS;
	bool IntraMolecular;
	bool TruncateAction, TruncateEnergy;
  bool useFirstImage;
  bool hacksOn;
	enum Interactions{LJ,CHARGE,SPRING};
	Array<bool, 2> Interacting;
	Array<string, 1> LJSpecies;
	Array<string, 1> ChargeSpecies;
	Array<string, 1> SpringSpecies;
	Array<string, 1> KineticSpecies;
	Array<string, 1> QuadSpecies;
	Array<string, 1> PairSpecies;
	Array<string, 1> CoreSpecies;
	Array<bool,2> Updated;
	Array<double,2> COMTable;
	Array<dVec,2> COMVecs;
	Array<double, 1> lambdas;
	Array<double,1> PtclCharge;
	Array<double,1> PtclEpsilon;
	Array<double,1> PtclSigma;
	int NumImages;
	// parameters for ST2 modulation function
	double RL, RU;
	// parameters for SPC quadratic intramolecular potential
	double b, c, d, D, rho, R_OH_0, R_HH_0, alpha;
	double Dyn2kcal;
  double omega;
	// this is kind of a hack
	ofstream outfile;
	bool special;

  Array<CubicSpline*, 2> PairVTable;
  Grid* grid;
  double pairCutoff;

  // parameters for hard core potential
  double StartCore;
  double A, LO, HI;
  string boundary;

public:
  // hack
  int TIPPIMC;
 	// Parameters that can be specified on input
	// Defaults in constructor
	double prefactor, conversion;
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
						const Array<int,1> &activeParticles, int level,
						bool with_truncations, bool isAction);
  dVec Force(int slice, int ptcl);
  dVec Force(int slice, int ptcl, int ptcl2);
  dVec Force(int slice, int ptcl, Array<int,1>& activePtcls);
  double CalcCutoff(int ptcl1, int ptcl2, int slice, double Rcmag);
	double COMSeparation(int slice,int ptcl1,int ptcl2);
	double S(double r);
	dVec COMVelocity(int sliceA, int sliceB, int ptcl);
	void SetNumImages(int num);
  string GetName();
 
	MoleculeInteractionsClass (PathDataClass &pathData); };

#endif
