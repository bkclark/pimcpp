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

#include "SpeciesClass.h"

bool SpeciesClass::Read(IOSectionClass &inSection)
{
  assert(inSection.ReadVar("Name",Name));
  assert(inSection.ReadVar("lambda",lambda));
  Charge = 0.0;
  inSection.ReadVar("Charge",Charge);
  pseudoCharge = Charge;
  inSection.ReadVar("PseudoCharge",pseudoCharge);
  inSection.ReadVar("ChargeSpread",chargeSpread);
  inSection.ReadVar("Epsilon",Epsilon);
  inSection.ReadVar("Sigma",Sigma);
  if(!inSection.ReadVar("isIon",isIon))
    isIon = false;
  assert(inSection.ReadVar("NumParticles",NumParticles));
  assert(inSection.ReadVar("NumDim",NumDim));
  assert(inSection.ReadVar("Type",Type));
  return true;
}

bool FermionClass::Read(IOSectionClass &inSection)
{
  bool success = SpeciesClass::Read(inSection);
  assert (inSection.ReadVar ("NodeType", NodeType));
  return success;
}

SpeciesClass* ReadSpecies(IOSectionClass &inSection)
{
  string statisticsString;
  SpeciesClass *mySpecies;
  assert(inSection.ReadVar("Statistics",statisticsString));
  if (statisticsString == "FERMION")
    mySpecies = new FermionClass();
  else if (statisticsString == "BOSON")
    mySpecies = new BosonClass();
  else if (statisticsString == "BOLTZMANNON")
    mySpecies = new BoltzmannonClass();
  else {
    cerr<<"Species Statistics Unknown "<<statisticsString<<endl;
    exit(1);
  }
  assert(mySpecies->Read(inSection));
  return mySpecies;
}





