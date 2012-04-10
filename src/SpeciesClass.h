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

#ifndef SPECIES_CLASS_H
#define SPECIES_CLASS_H

#include "Common.h"
#include <Common/IO/IO.h>

using namespace IO;

typedef enum {FERMION, BOSON, BOLTZMANNON, ANYON} ParticleType;


/// This is an base class that holds all the information about
/// identical particles.  It may be specialized to hold specialized
/// information about particular types of particles.
class SpeciesClass
{
public:
  string Name, Type, NodeType;
  /// FirstPtcl and LastPtcl are inclusive
  int LastPtcl;
  int FirstPtcl;
  int NumDim;
  int NumParticles;
  Array<int,1> Ptcls;
  TinyVector <bool,NDIM> DimensionActive;
  virtual bool Read(IOSectionClass &inSection);
  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  dVec assymetric_lambda;
  double Charge;
  double chargeSpread;
  // I added this parameter for use with empirical potentials
  // in a pre-rejection context, where having the "real" charge
  // is also necessary - john
  double pseudoCharge;
  /// sigma and epsilon are parameters for the Lennard-Jones potential used in the TIP5P water model. -jg
  double Sigma;
  double Epsilon; 

  // OBSOLETE
	////  If specified in the input file,
	////  molecule holds the name of a "molecule"
  ////  to which particles of this species belong
	//// 	formula holds an integer to designate
	////	the corresponding	chemical formula
	//string molecule; 
	//int formula;
	//bool AssignMoleculeIndex;

  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual ParticleType GetParticleType() = 0;
};

SpeciesClass* ReadSpecies(IOSectionClass &inSection);



///This is the inherited class that holds the information about
///the electrons. Eventually this will probably be turned into
///a fermion class.  Currently, no actual information about the electrons
///are actually included.
class FermionClass : public SpeciesClass
{
public:
  bool Read(IOSectionClass &inSection);
  FermionClass()  { 
    NumDim=NDIM;  
    DimensionActive=true;
  };
  ~FermionClass() { };
  ParticleType GetParticleType(){ return FERMION; }
};


///This is the inherited class that holds the information about
///the protons.  Currently no useful information about protons is contained
///here
class BosonClass : public SpeciesClass
{
public:
  ParticleType GetParticleType(){ return BOSON; }
  ///Just returns 0 until we do something more intelligble.
  BosonClass() 
  {
    NumDim=NDIM;
    DimensionActive=true;
  }
  ~BosonClass()
  {
  }
  
};

class BoltzmannonClass : public SpeciesClass
{
public:
  ParticleType GetParticleType(){ return BOLTZMANNON; }

  BoltzmannonClass()
  {
    NumDim=NDIM;  
    DimensionActive=true;
  }
  
  ~BoltzmannonClass()
  {
  }
  
};

#endif

