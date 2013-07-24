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

#include "HBond.h"
#include <iostream>
#include <fstream>
#include <string>
#include "../Moves/MoveUtils.h"

//////////////////////////////////////
///  Hydrogen Bond Analysis Class  ///
//////////////////////////////////////

bool HbondClass::IsHBond(int slice, int OA, int OB){
  bool bond = false;
	dVec OO;
	double OOmag;
  PathData.Path.DistDisp(slice,OA,OB,OOmag,OO);
	if(OOmag < OOlimit){
 		if (CheckPair(slice,OA,OB,Protons(OB,0)))
  	  bond = true;
  	else if (CheckPair(slice,OA,OB,Protons(OB,1)))
  	  bond = true;
  	else if (CheckPair(slice,OB,OA,Protons(OA,0)))
  	  bond = true;
  	else if (CheckPair(slice,OB,OA,Protons(OA,1)))
  	  bond = true;
	}
  return bond;
}

bool HbondClass::CheckPair(int slice, int obond, int ohome, int p){
  bool bond = false;
  double OHmag,OHbondmag;
  dVec OH,OHbond;
  PathData.Path.DistDisp(slice,p,obond,OHbondmag,OHbond);
  PathData.Path.DistDisp(slice,p,ohome,OHmag,OH);
  double theta = GetAngle(OHbond,OH);
  if (theta > HOHangle)
    bond = true;
  return bond; 
}

void HbondClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("Frequency",Frequency));
	totalSlices = PathData.NumTimeSlices()-1;
	totalMol = PathData.Mol.NumMol();
  AgeTable.resize(totalSlices, totalMol, totalMol);
  AgeTable=0;
	cerr << "resized AgeTable: " << AgeTable.size() << endl;

  BondCount.resize(PathData.Mol.NumMol(),PathData.Mol.NumMol());
  BondCount=0;
	cerr << "resized BondCount: " << BondCount.size() << endl;

	string ProtonSpecies = "H";
	in.ReadVar("ProtonSpecies",ProtonSpecies);
	Protons.resize(PathData.Mol.NumMol(), 2);
	int totalFound = 0;
	for(int m=0; m<PathData.Mol.NumMol(); m++){
		int foundIndex = 0;
  	for (int a = 0; a < PathData.Mol.MembersOf(m).size(); a++){
  	  int ptcl = PathData.Mol.MembersOf(m)(a);
			if(PathData.Path.ParticleSpeciesNum(ptcl) == PathData.Path.SpeciesNum(ProtonSpecies)){
				Protons(m, foundIndex) = ptcl;
				foundIndex++;
				totalFound++;
			}
		}
	}
	assert(totalFound != 0);

  bool readStartGrid=in.ReadVar("start",gridStart);
  if (!readStartGrid)
    gridStart=0.0;
  assert(in.ReadVar("BinSize",numGridPoints));
	gridEnd = gridStart + numGridPoints*Frequency;

  Histogram.resize(numGridPoints);
  Histogram = 0;
  LifetimeHist.resize(numGridPoints);
  LifetimeHist = 0;
  AgeHist.resize(numGridPoints);
  AgeHist = 0;

  //in.CloseSection();
}

void HbondClass::WriteBlock()
{
	cerr << "HBOND WRITEBLOCK totalslices is " << totalSlices << " and totalSamples is " << TotalSamples << endl;
  double AgeNorm = 1.0/(totalSlices*TotalSamples);
  double LifetimeNorm = 1.0/(totalSlices);
	double totalHBS = 0.0;
	double hbs = 0.0;

  Array<double,1> LifetimeArray(LifetimeHist.size());
  for (int i=0; i<LifetimeArray.size(); i++){
    LifetimeArray(i) = (double)LifetimeHist(i)*LifetimeNorm;
		totalHBS += LifetimeArray(i);
  }
  LifetimeVar.Write(LifetimeArray);
  LifetimeVar.Flush();
	//cerr << "HBOND: Histogram contains " << totalHBS << " HBonds; hist is " << LifetimeHist << ", norm is " << LifetimeNorm << endl;

  Array<double,1> AgeArray(AgeHist.size());
  for (int i=0; i<AgeArray.size(); i++){
    AgeArray(i) = (double)AgeHist(i)*AgeNorm;
		hbs += AgeArray(i);
  }
	cerr << "HBOND: Snapshot has " << hbs << " bonds" << endl;
  AgeVar.Write(AgeArray);
  AgeVar.Flush();
	AgeHist = 0;
	TotalSamples = 0;
}

void HbondClass::Accumulate()
{
	TimesCalled++;
	TotalSamples++;
	// Measure and tabulate hbonds at the curent time step
	int bondCt = 0;
	int totalPair = 0;
	int died=0;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
  	/// loop over molecules 
    for (int mol1=0;mol1<PathData.Mol.NumMol()-1;mol1++){
      for (int mol2=mol1+1;mol2<PathData.Mol.NumMol();mol2++){
				totalPair++;
        if(IsHBond(slice,mol1,mol2)){
          AgeTable(slice, mol1, mol2)++;
					int age = AgeTable(slice,mol1,mol2);
					if(age>oldest)
						oldest = age;
					if(age>=numGridPoints)
						age = numGridPoints - 1;
					AgeHist(age)++;
					bondCt++;
          //BondCount(mol1,mol2)++;
          //BondCount(mol2,mol1)++;
        }
				else{
					int index = AgeTable(slice, mol1, mol2);
					if(index >= numGridPoints){
						index = numGridPoints-1;
						overflow++;
					}
					if(index>0){
						LifetimeHist(index)++;
						died++;
					}
					AgeTable(slice, mol1, mol2) = 0;
				}
      }
    }
	}

	if(TimesCalled%100 == 0){
		cerr << "HBond: " << overflow << " are longer than the histogram size" << endl;
		cerr << TimesCalled << ": HBOND counted " << bondCt << " of " << totalPair << " possible HBonds" << " and " << died << " died with lifetime>0" << endl;
		cerr << "Oldest HB is " << oldest << endl;
  //	// Write files
  //	WriteBlock();
	}

  //// Tabulate the number of hbonds and accumulate in Histogram
  //int count = 0;
  //for (int mol = 0; mol < BondCount.extent(0); mol++){
  //  // Sum entries over rows
  //  for (int m = 0; m < BondCount.extent(0); m++)
  //    count += BondCount(m,mol);
  //  // Sum entries over columns
  //  for (int n = 0; n < BondCount.extent(1); n++)
  //    count += BondCount(mol,n);
  //  Histogram(count)++;
  //  count = 0;
  //}
  //BondCount = 0;
  //Histogram = 0;
}

void HbondClass::Initialize()
{
  TimesCalled = 0;
  TotalSamples = 0;
	overflow = 0;
	oldest = 0;

	// Here we use HBond criteria r_OO < 3.5 and HOH bond angle > 145 deg, after Artacho (2004)
  OOlimit = 3.5;
  HOHangle = 2*M_PI*145/360;
}
