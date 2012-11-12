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

#include "ObservableBase.h"
#include "time.h"


// void
// ObservableVar::Flush()
// {
//   if (Comm.MyProc() == 0)
//     Out.FlushFile();
// }


void ObservableClass::WriteInfo()
{
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.WriteVar("Description",Description);
}


void ObservableClass::Read(IOSectionClass &in)
{
  in.ReadVar("Prefactor", Prefactor);
  if (in.ReadVar("TemporalFrequency",TemporalFrequency)) {
    gettimeofday(&starttime, &tzone);
    Frequency = -1; // Probably a better way
  } else {
    assert(in.ReadVar("Frequency",Frequency));
    TemporalFrequency = -1;
  }
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
}


void ObservableClass::DoEvent()
{
  TimesCalled++;
  struct timeval endtime;
  gettimeofday(&endtime, &tzone);
  double TimeDiff = (double)(endtime.tv_sec-starttime.tv_sec) + 1.0e-6*(double)(endtime.tv_usec-starttime.tv_usec);
  if ((Frequency > 0 && (TimesCalled % Frequency) == 0) || (TemporalFrequency > 0 && TimeDiff > TemporalFrequency)) {
    Accumulate();
    gettimeofday(&starttime,&tzone);
  }
}


// Calculates weight from sign and importance sampling
double ObservableClass::CalcFullWeight()
{
  double TempSign;
  double currSign = PathData.Path.SignWeight;
  PathData.Path.Communicator.GatherProd(currSign, TempSign, 0);
  PathData.Path.Communicator.Broadcast(0, TempSign);
  double FullSign = TempSign;

  double FullNodeWeight = 0.0;
  if (PathData.Path.UseNodeImportance) {
    double TempNodeWeight;
    double NodeWeight = 0.0;
    PathData.Actions.GetNodalActions(NodeWeight);
    PathData.Path.Communicator.GatherProd(NodeWeight, TempNodeWeight, 0);
    PathData.Path.Communicator.Broadcast(0, TempNodeWeight);
    FullNodeWeight = TempNodeWeight;
  }

  double FullWeight = exp(-FullNodeWeight)*FullSign;
  //cout << PathData.Path.CloneStr << " FW : " << FullWeight << " " << NodeWeight << " " << FullSign << endl;
  return 1.0/FullWeight;
}


// Permutation Counting Things

void ObservableClass::SetupPermSectors(int n, int MaxNSectors)
{
  if (PathData.Path.Communicator.MyProc() == 0) {
    cout << PathData.Path.CloneStr << " Setting up Permutation Sectors" << endl;
    vector<int> a;
    a.resize(n);
    for (int i=0; i<n; i++) {
      a[i] = 0;
    }
    int k = 1;
    int y = n-1;
    vector< vector<int> > tmpPossPerms;
    while (k != 0 && (MaxNSectors > PossPerms.size() || !MaxNSectors)) {
      int x = a[k-1] + 1;
      k -= 1;
      while (2*x <= y) {
        a[k] = x;
        y -= x;
        k += 1;
      }
      int l = k+1;
      while (x <= y && (MaxNSectors > PossPerms.size() || !MaxNSectors)) {
        a[k] = x;
        a[l] = y;
        vector<int> b;
        for (vector<int>::size_type j=0; j!=k+2; j++)
          b.push_back(a[j]);
        tmpPossPerms.push_back(b);
        x += 1;
        y -= 1;
      }
      a[k] = x+y;
      y = x+y-1;
      vector<int> c;
      for (vector<int>::size_type j=0; j!=k+1; j++)
        c.push_back(a[j]);
      tmpPossPerms.push_back(c);
    }


    cout << PathData.Path.CloneStr << " Putting sectors in map" << endl;
    for (vector<int>::size_type j=0; j != tmpPossPerms.size(); j++) {
      sort(tmpPossPerms[j].begin(),tmpPossPerms[j].end());
      PossPerms[tmpPossPerms[j]] = j;
    }
  }

  //sort(PossPerms.begin(),PossPerms.end(),CompareVectors());

}


void ObservableClass::GetPermInfo(vector<int> &ThisPerm, int &PermSector, int &PermNumber)
{
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
  PossPermsIterator = PossPerms.find(ThisPerm);
  if (PossPermsIterator == PossPerms.end()) {
    cerr << "Broken Permutation: " << endl;
    for (vector<int>::size_type i=0; i != ThisPerm.size(); i++)
      cerr << ThisPerm[i] << " ";
    cerr << endl;
    exit(1);
  } else {
    PermSector = PossPermsIterator->second;
    return;
  }
}
