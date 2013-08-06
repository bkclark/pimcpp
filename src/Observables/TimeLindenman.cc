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

#include "TimeLindenman.h"


void 
TimeLindenmanClass::ProduceTimeMatrix(int mcStep)
{
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
    uj_minus_ujp(mcStep,ptcl).clear();
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    ///now loop over nearby ptcls
    for (int nearPtcl=ptcl+1;nearPtcl<PathData.Path.NumParticles();nearPtcl++){
      double dist; dVec diff;
      dVec r1=CentroidPos(ptcl);
      dVec r2=CentroidPos(nearPtcl);
      diff =r2-r1;
      PathData.Path.PutInBox(diff);
      dist=sqrt(dot(diff,diff));
      if (dist<DistCutoff){ //the particles are considered nearby
	pair<int, dVec> nearPair(nearPtcl,diff);
	dVec mDiff=diff;
	mDiff=-1.0*mDiff;
	pair<int, dVec> nearPair2(ptcl,mDiff);
	uj_minus_ujp(mcStep,ptcl).insert(nearPair);
	uj_minus_ujp(mcStep,nearPtcl).insert(nearPair2);
      }
    }
  }
}

///centroidPos must be resized to the number of particles
///when this function is called

///not sure this gets the right centroid in parallel
void 
TimeLindenmanClass::CalculateCentroid()
{
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    dVec centroid=0.0;
    for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
      dVec disp=PathData.Path.Velocity(0,slice,ptcl);
      centroid += disp;
    }
    centroid = centroid * (1.0/ (PathData.Path.NumTimeSlices()-1));
    centroid = centroid + PathData.Path(0,ptcl);
    CentroidPos(ptcl)=centroid;
  }  
}


///centroidPos must be resized to the number of particles
///when this function is called

///not sure this gets the right centroid in parallel
void 
TimeLindenmanClass::CalculateCentroid_parallel()
{


  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    dVec centroid;
    dVec zeroVec=PathData.Path(0,ptcl);
    PathData.Path.Communicator.Broadcast(0,zeroVec);
    dVec localCentroid=0.0;
    for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++){
      dVec disp=PathData.Path.MinImageDisp(zeroVec,PathData.Path(slice,ptcl));
      localCentroid += disp;
    }
    for (int dim=0;dim<NDIM;dim++)
      centroid(dim)=PathData.Path.Communicator.Sum(localCentroid(dim));
    centroid = centroid * (1.0/ (PathData.Path.TotalNumSlices));
    centroid = centroid + zeroVec;
    CentroidPos(ptcl)=centroid;
  }  
}



double 
TimeLindenmanClass::CalcValue(int time1, int time2)
{
  map<int,dVec>::iterator ptclIter;
  double term1=0.0; double term2=0.0; double term3=0.0;
  for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++)
    for (ptclIter=uj_minus_ujp(time2,ptcl1).begin();ptclIter!=uj_minus_ujp(time2,ptcl1).end();
	 ptclIter++){
      int ptcl2=(*ptclIter).first;
      dVec vec_12_t2=(*ptclIter).second;
      if (uj_minus_ujp(time1,ptcl1).count(ptcl2)>0){
	term1+=dot(vec_12_t2,vec_12_t2);
	dVec vec_12_t1=uj_minus_ujp(time1,ptcl1)[ptcl2];
	term2+=dot(vec_12_t1,vec_12_t1);
	term3+=0.0-2.0*dot(vec_12_t1,vec_12_t2);
      }
    }
  return term1+term2+term3;
}


void
TimeLindenmanClass::Accumulate()
{


  if (!FullyWrapped){
    TotalCurrentData++;
    if (TotalCurrentData>=NumStoredMCSteps)
      FullyWrapped=true;
  }
  
//   CalculateCentroid();
//   Array<dVec,1> Centroid_backup(CentroidPos.size());
//   for (int i=0;i<CentroidPos.size();i++)
//     Centroid_backup(i)=CentroidPos(i);
  CalculateCentroid_parallel();
//   for (int i=0;i<CentroidPos.size();i++)
//     for (int dim=0;dim<NDIM;dim++)
//       assert(Centroid_backup(i)[dim]-CentroidPos(i)[dim]<1e-10);
  CurrTime=(CurrTime+1) % NumStoredMCSteps;
  ProduceTimeMatrix(CurrTime);
  int maxToGo=min(NumStoredMCSteps,TotalCurrentData);
  for (int i=0;i<maxToGo;i++){
    TimeDisp(i)+=CalcValue(CurrTime,(CurrTime-i +NumStoredMCSteps) % NumStoredMCSteps);
    NumStepArray(i)=NumStepArray(i)+1.0;
  }
  
}


void
TimeLindenmanClass::WriteBlock()
{
//   for (int i=0;i<TimeDisp.size();i++)
//     if (NumStepArray(i)!=0)
//       TimeDisp(i)=TimeDisp(i)/NumStepArray(i);
  TimeDispVar.Write(TimeDisp);
  NumStepVar.Write(NumStepArray);
  TimeDisp=0.0;
  NumStepArray=0.0;
  ////  CurrTime=0;
}

void 
TimeLindenmanClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  IOSection.NewSection("grid");
  grid.Write(IOSection);
  IOSection.CloseSection();

  int numBins = grid.NumPoints-1;
  Array<double,1> r(numBins);
  for (int i=0; i<numBins; i++) {
    double ra = grid(i);
    double rb = grid(i+1);
    r(i) = 0.6666 * (rb*rb*rb-ra*ra*ra)/(rb*rb-ra*ra);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "g(r)");
  IOSection.WriteVar("Type","CorrelationFunction");
  IOSection.WriteVar("Cumulative", false);
}


void 
TimeLindenmanClass::ReadGrid(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  assert(gridType=="Linear");
  bool readStartGrid=in.ReadVar("start",gridStart);
  bool readEndGrid=in.ReadVar("end",gridEnd);
  if (!readStartGrid)
    gridStart=0.0;
  if (!readEndGrid){
    if (PathData.Path.GetPeriodic()[0]){
      gridEnd=PathData.Path.GetBox()[0];
    }
    else {
      cerr<<"I don't know where you want me to end this grid"<<endl;
      assert(1==2);
    }
  }
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  NumSamples=0;
  Histogram.resize(numGridPoints-1);
  HistSum.resize(Histogram.size());
  HistDouble.resize(Histogram.size());
  Histogram=0;
  in.CloseSection();
}

void
TimeLindenmanClass::Read (IOSectionClass &in)
{
  NumSamples=0;
  ParticleOrder.resize(PathData.Path.NumParticles());
  ///It's probably important that the grid is the same grid that is in
  ///the pair correlation function. Not sure how to authenticate this.
  ReadGrid(in);

  if (!in.ReadVar("Centroid",Centroid))
      Centroid=false;
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();
  assert(in.ReadVar("TimeWindow",TimeWindow));
  assert(in.ReadVar("TimeResolution",TimeResolution));
  NumStoredMCSteps=TimeWindow/TimeResolution;
  uj_minus_ujp.resize(NumStoredMCSteps,PathData.Path.NumParticles());
  TimeDisp.resize(NumStoredMCSteps);
  NumStepArray.resize(NumStoredMCSteps);
  CentroidPos.resize(PathData.Path.NumParticles());
}
