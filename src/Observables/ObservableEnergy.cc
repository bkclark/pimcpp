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
#include <algorithm>
#include <numeric>

// These are included for the new runtime 
// specification of energy observables 
// to compute; see below -John
#include "../Actions/MoleculeInteractionsClass.h"
#include "../Actions/EAMClass.h"
#include "../Actions/ST2WaterClass.h"
#include "../Actions/QMCSamplingClass.h"
#include "../Actions/DavidLongRangeClassYk.h"
#include "../Actions/ShortRangeOn_diagonal_displace_Class.h"


void EnergyClass::SetupPermSectors(int n)
{
  vector<int> a;
  a.resize(n);
  for (int i=0; i<n; i++) {
    a[i] = 0;
  }
  int k = 1;
  int y = n-1;
  while (k != 0) {
    int x = a[k-1] + 1;
    k -= 1;
    while (2*x <= y) {
      a[k] = x;
      y -= x;
      k += 1;
    }
    int l = k+1;
    while (x <= y) {
      a[k] = x;
      a[l] = y;
      vector<int> b;
      for (vector<int>::size_type j=0; j!=k+2; j++)
        b.push_back(a[j]);
      PossPerms.push_back(b);
      x += 1;
      y -= 1;
    }
    a[k] = x+y;
    y = x+y-1;
    vector<int> c;
    for (vector<int>::size_type j=0; j!=k+1; j++)
      c.push_back(a[j]);
    PossPerms.push_back(c);
  }

  for (vector<int>::size_type j=0; j != PossPerms.size(); j++) {
    sort(PossPerms[j].begin(),PossPerms[j].end());
  }

}


void EnergyClass::GetPermInfo(int &PermSector, int &PermNumber)
{
  vector<int> ThisPerm;
  PermNumber = 0;
  PathClass & Path = PathData.Path;
  int N = PathData.Path.NumParticles();
  if (CountedAlready.size() != N) {
    CountedAlready.resize(N);
    TotalPerm.resize(N);
  }
  PathData.Path.TotalPermutation(TotalPerm);
  CountedAlready = false;
  int ptcl = 0;
  /// Only proc 0 gets TotalPerm
  if (Path.Communicator.MyProc() == 0) {
    while (ptcl < N) {
      if (!CountedAlready(ptcl)) {
        int startPtcl = ptcl;
        int roamingPtcl = ptcl;
        int cycleLength = 0;
        roamingPtcl = TotalPerm(roamingPtcl);
        while (roamingPtcl != startPtcl) {
          CountedAlready(roamingPtcl) = true;
          cycleLength++;
          roamingPtcl = TotalPerm(roamingPtcl);
        }
        ThisPerm.push_back(cycleLength+1);
        PermNumber += cycleLength;
      }
      ptcl++;
    }
  } else
    return;

  sort(ThisPerm.begin(),ThisPerm.end());
  for (vector<int>::size_type i=0; i != PossPerms.size(); i++)
    if (ThisPerm == PossPerms[i]) {
      PermSector = i;
      return;
    }

  // Broken Permutation!
  cerr << "Broken Permutation: " << endl;
  for (vector<int>::size_type i=0; i != ThisPerm.size(); i++)
    cerr << ThisPerm[i] << " ";
  cerr << endl;

  exit(1);
  return;
}


// Fix to include final link between link M and 0
void EnergyClass::Accumulate()
{
  TimesCalled++;
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices() - 1);
  double FullWeight;
  if (TrackSign) {
    double currWeight = PathData.Path.Weight;
    PathData.Path.Communicator.GatherProd(currWeight, FullWeight, 0);
  } else
    FullWeight = 1;
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
  int PermSector, PermNumber;
  GetPermInfo(PermSector,PermNumber);
  if (Path.Communicator.MyProc() == 0) {
    PermTotalSum(PermSector) += localSum;
    EnergyVals(PermNumber) += localSum;
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
      //      DavidLongRangeClassYk *lr = (DavidLongRangeClassYk*)(PathData.Actions.GetAction("DavidLongRange"));
      DavidLongRangeClassYk2 *lr = (DavidLongRangeClassYk2 *) (&(PathData.Actions.DavidLongRange));
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
  EnergyVals = Prefactor * EnergyVals * norm;
  EnergyValsVar.Write(EnergyVals);
  for (int i = 0; i < PermTotalSum.size(); i++)
    PermTotalSum(i) = Prefactor * PathData.Path.Communicator.Sum(PermTotalSum(i) * norm);
  PermTotalSumVar.Write(PermTotalSum);
  for (int i = 0; i < EnergyHistogram.histogram.size(); i++)
    EnergyHistogram.histogram[i] = Prefactor * EnergyHistogram.histogram[i] * norm * nslices;
  Array <double,1> EnergyHistogramTemp(&(EnergyHistogram.histogram[0]), shape(EnergyHistogram.histogram.size()), neverDeleteData);
  EnergyHistogramVar.Write(EnergyHistogramTemp);
  for (int n = 0; n < numEnergies; n++) {
    OtherVars[n]->Write(Prefactor * PathData.Path.Communicator.Sum(OtherSums[n]) * norm);
    OtherSums[n] = 0.0;
  }

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
  PermTotalSum = 0.0;
  Residual = 0.0;
  EnergyVals = 0.0;
  NumSamples = 0;
  EnergyHistogram.Clear();
}


void EnergyClass::Read(IOSectionClass & in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("Frequency", Freq));
  if(!in.ReadVar("TrackSign", TrackSign))
    TrackSign = 0;
  if (PathData.Path.Communicator.MyProc() == 0) {
    WriteInfo();
    IOSection.WriteVar("Type", "Scalar");
  }

  // Perm Sector Counting
  SetupPermSectors(PathData.Path.NumParticles());
  PermTotalSum.resize(PossPerms.size());
  PermTotalSum = 0.0;
  int PermSector, PermNumber;
  GetPermInfo(PermSector,PermNumber);
  cout << "Starting in Perm Sector " << PermSector << " of " << PossPerms.size()-1 << endl;

  // Cycle Length Counting
  EnergyVals.resize(PathData.Path.NumParticles() * 2);
  EnergyVals = 0.0;

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
