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
//Changed Something again
#include "ObservableEnergy.h"

// These are included for the new runtime 
// specification of energy observables 
// to compute; see below -John
#include "../Actions/MoleculeInteractionsClass.h"
#include "../Actions/EAMClass.h"
#include "../Actions/ST2WaterClass.h"
#include "../Actions/QMCSamplingClass.h"
#include "../Actions/DavidLongRangeClassYk.h"
#include "../Actions/ShortRangeOn_diagonal_displace_Class.h"


// Fix to include final link between link M and 0
void EnergyClass::Accumulate()
{
  TimesCalled++;
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices() - 1);

  // Get the Full Weight from sign and importance sampling
  double FullWeight = CalcFullWeight();
  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, dUNonlocal, residual;
  PathData.Actions.Energy(kinetic, dUShort, dULong, node, vShort, vLong,
                          dUNonlocal, residual);
  double localSum = (kinetic + dUShort + dULong + node + dUNonlocal) * FullWeight;
  TotalSum += localSum;
  KineticSum += kinetic * FullWeight; /* * PathData.Path.Weight */ ;
  dUShortSum += dUShort * FullWeight; /* * PathData.Path.Weight */ ;
  dULongSum += dULong * FullWeight; /* * PathData.Path.Weight */ ;
  NodeSum += node * FullWeight; /* * PathData.Path.Weight */ ;
  VShortSum += vShort * FullWeight; /* * PathData.Path.Weight */ ;
  VLongSum += vLong * FullWeight; /* * PathData.Path.Weight */ ;
  dUNonlocalSum += dUNonlocal * FullWeight;
  Residual += residual;

  // Permutation Counting
  if(CountPerms) {
    int PermSector, PermNumber;
    vector<int> ThisPerm;
    GetPermInfo(ThisPerm,PermSector,PermNumber);
    if (Path.Communicator.MyProc() == 0) {
      PermEnergy(PermSector) += localSum;
      PermHist(PermSector) += 1;
      EnergyVals(PermNumber) += localSum;
    }
  }

  // Other Energies
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices() - 1;
  for (int n = 0; n < numEnergies; n++) {
    double otherE = OtherActions[n]->d_dBeta(slice1, slice2, 0);
    OtherSums[n] += otherE;
    localSum += otherE;
    TotalSum += otherE;
  }

  // Energy Histogram
  double completeSum = PathData.Path.Communicator.Sum(localSum) /
                       (double) PathData.Path.TotalNumSlices;
  EnergyHistogram.add(PathData.Path.Communicator.Sum(localSum) /
                      (double) PathData.Path.TotalNumSlices, 1.0);
}


void EnergyClass::ShiftData(int NumTimeSlices)
{
  // Do nothing
}


void EnergyClass::WriteBlock()
{
  if (FirstTime) {
    FirstTime = false;
    Array<double,1> vtail;
    vtail.resize(PathData.Actions.PairArray.size());
    double longrange_vtail = 0.0;
    for (int i = 0; i < PathData.Actions.PairArray.size(); i++)
      vtail(i) = ((DavidPAClass *) (PathData.Actions.PairArray(i)))->Vimage;
    if (PathData.Path.DavidLongRange) {
      DavidLongRangeClassYk *lr = (DavidLongRangeClassYk *) (&(PathData.Actions.DavidLongRange));
      //longrange_vtail = 0.5 * lr->yk_zero(0) * PathData.Path.NumParticles() / Path.GetVol();
      longrange_vtail = 0.5 * lr->yk_zero(0) * PathData.Path.NumParticles();
    }
    VTailSRVar.Write(vtail);
    VTailLRVar.Write(longrange_vtail);
    HistStart.Write(EnergyHistogram.startVal);
    HistEnd.Write(EnergyHistogram.endVal);
    NumPoints.Write(EnergyHistogram.NumPoints);
  }

  int nslices = PathData.Path.TotalNumSlices;
  double norm = 1.0 / ((double) NumSamples * (double) nslices);

  TotalVar.Write(Prefactor * PathData.Path.Communicator.Sum(TotalSum) * norm);
  KineticVar.Write(Prefactor * PathData.Path.Communicator.Sum(KineticSum) * norm);
  dUShortVar.Write(Prefactor * PathData.Path.Communicator.Sum(dUShortSum) * norm);
  dULongVar.Write(Prefactor * PathData.Path.Communicator.Sum(dULongSum) * norm);
  NodeVar.Write(Prefactor * PathData.Path.Communicator.Sum(NodeSum) * norm);
  VShortVar.Write(Prefactor * PathData.Path.Communicator.Sum(VShortSum) * norm);
  VLongVar.Write(Prefactor * PathData.Path.Communicator.Sum(VLongSum) * norm);
  dUNonlocalVar.Write(Prefactor * PathData.Path.Communicator.Sum(dUNonlocalSum) * norm);
  ResidualVar.Write(Prefactor * PathData.Path.Communicator.Sum(Residual) * norm);
  for (int n = 0; n < numEnergies; n++) {
    OtherVars[n]->Write(Prefactor * PathData.Path.Communicator.Sum(OtherSums[n]) * norm);
    OtherSums[n] = 0.0;
  }

  // Permutation Counting
  if (CountPerms) {
    EnergyVals = Prefactor * EnergyVals * norm;
    EnergyValsVar.Write(EnergyVals);
    for (int i = 0; i < PermEnergy.size(); i++) {
      PermEnergy(i) = Prefactor * PathData.Path.Communicator.Sum(PermEnergy(i) / ((double)nslices*(double)PermHist(i)));
      if (isnan(PermEnergy(i)))
        PermEnergy(i) = 0.0;
      PermHist(i) = Prefactor * PathData.Path.Communicator.Sum(PermHist(i) / (double)NumSamples);
    }
    PermEnergyVar.Write(PermEnergy);
    //PermHistVar.Write(PermHist); May want to do this later
    PermEnergy = 0.0;
    PermHist = 0;
  }

  // Energy Histogram
  for (int i = 0; i < EnergyHistogram.histogram.size(); i++)
    EnergyHistogram.histogram[i] = Prefactor * EnergyHistogram.histogram[i] * norm * nslices;
  Array <double,1> EnergyHistogramTemp(&(EnergyHistogram.histogram[0]), shape(EnergyHistogram.histogram.size()), neverDeleteData);
  EnergyHistogramVar.Write(EnergyHistogramTemp);

  if (PathData.Path.Communicator.MyProc() == 0)
    IOSection.FlushFile();

  TotalSum = 0.0;
  KineticSum = 0.0;
  dUShortSum = 0.0;
  dULongSum = 0.0;
  NodeSum = 0.0;
  VShortSum = 0.0;
  VLongSum = 0.0;
  dUNonlocalSum = 0.0;
  Residual = 0.0;
  EnergyVals = 0.0;
  NumSamples = 0;
  EnergyHistogram.Clear();
}


void EnergyClass::Read(IOSectionClass & in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("Frequency", Freq));

  if (PathData.Path.Communicator.MyProc() == 0) {
    WriteInfo();
    IOSection.WriteVar("Type", "Scalar");
  }

  // Sign Tracking
  if(!in.ReadVar("TrackSign", TrackSign))
    TrackSign = 0;

  // Perm Sector Counting
  if(!in.ReadVar("CountPerms",CountPerms))
    CountPerms = 0;
  if(CountPerms) {
    int N = PathData.Path.NumParticles();
    // Maximum number of permutation sectors tracked
    int MaxNSectors;
    if(!in.ReadVar("MaxNSectors", MaxNSectors))
      MaxNSectors = 0; // 0 -> Track all sectors
    SetupPermSectors(N,MaxNSectors);
    PermEnergy.resize(PossPerms.size());
    PermEnergy = 0.0;
    PermHist.resize(PossPerms.size());
    PermHist = 0;
    EnergyVals.resize(N*2);
    EnergyVals = 0.0;
    int PermSector, PermNumber;
    vector<int> ThisPerm;
    GetPermInfo(ThisPerm,PermSector,PermNumber);
    cout << PathData.Path.CloneStr << " Starting in Perm Sector " << PermSector << " of " << PossPerms.size()-1 << endl;
  }

  // Other Energies
  Array <string,1> EnergyStrings(0);
  in.ReadVar("ComputeEnergies", EnergyStrings);
  numEnergies = EnergyStrings.size();
  OtherActions.resize(numEnergies);
  OtherVars.resize(numEnergies);
  OtherSums.resize(numEnergies);
  for (int n = 0; n < numEnergies; n++) {
    OtherActions[n] = PathData.Actions.GetAction(EnergyStrings(n));
    cerr << "Energy observable added action with label " << EnergyStrings(n) << endl;
    OtherVars[n] = new ObservableDouble(EnergyStrings(n), IOSection, PathData.Path.Communicator);
    OtherSums[n] = 0.0;
  }

  // Energy Histogram
  double histStart = 0.0;
  double histEnd = 1.0;
  int histPoints = 1;
  in.ReadVar("HistStart", histStart);
  in.ReadVar("HistEnd", histEnd);
  in.ReadVar("HistPoints", histPoints);
  EnergyHistogram.Init(histPoints, histStart, histEnd);
  EnergyHistogramSum.resize(EnergyHistogram.histogram.size());
}
