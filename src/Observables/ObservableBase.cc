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
  PathData.Path.Communicator.Broadcast(0, TimeDiff);
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


void ObservableClass::GetStats(vector<double>& xs, double& mean, double& err, double& N)
{
  vector<double>::iterator it;
  double totx = 0.0;
  double totx2 = 0.0;
  for(it = xs.begin(); it != xs.end(); it++) {
    double x = *it;
    totx += x;
    totx2 += x*x;
  }
  N = xs.size();
  mean = totx/N;
  double mean2 = totx2/N;
  double sampVar = mean2 - mean*mean;
  double var = 0.0;
  if (N > 1)
    var = ((N+0.0)/(N-1.0)) * sampVar;
  double sigma = sqrt(var);
  err = sigma/sqrt(N);
}
