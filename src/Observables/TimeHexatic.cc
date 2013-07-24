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

#include "TimeHexatic.h"
#include <complex>

double TimeHexaticClass::unit2angle(double x,double y)
{
  double angle=atan(y/x);
  if (x<0)
    angle=angle+M_PI;
  return angle;
}

complex<double>
TimeHexaticClass::OrderParamater(Array<dVec,1> &centroidPos,int ptcl)
{
  complex<double> op=0.0;
  for (int nearPtcl=0;nearPtcl<PathData.Path.NumParticles();
       nearPtcl++){
    double r12dist;
    dVec r12disp;
    if (nearPtcl!=ptcl){
      r12disp=centroidPos(nearPtcl)-centroidPos(ptcl);
      r12dist=sqrt(dot(r12disp,r12disp));
      if (r12dist<DistCutoff){
	r12disp=r12disp * (1.0/r12dist);
	if (abs(dot(r12disp,r12disp)-1.0)>=0.001)
	  cerr<<dot(r12disp,r12disp);
	if (!((dot(r12disp,r12disp)-1.0)<=0.001))
	  cerr<<"GRRRL: "<<dot(r12disp,r12disp)<<" "<<r12dist<<" "<<ptcl<<" "<<nearPtcl<<endl;
	assert(abs(dot(r12disp,r12disp)-1.0)<0.001);
	double theta_12=unit2angle(r12disp(0),r12disp(1));
	op=op+complex<double>(cos(theta_12*q),sin(theta_12*q));
      }
    }
  }
  return op;
}




void 
TimeHexaticClass::ProduceTimeMatrix(int mcStep)
{
  assert(mcStep<PosArray.extent(0));
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    OPArray(mcStep,ptcl)=OrderParamater(CentroidPos,ptcl);
    PosArray(mcStep,ptcl)=CentroidPos(ptcl);
  }
}

///centroidPos must be resized to the number of particles
///when this function is called

///not sure this gets the right centroid in parallel
void 
TimeHexaticClass::CalculateCentroid()
{
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    dVec centroid=0.0;
    for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
      dVec disp=PathData.Path.Velocity(0,slice,ptcl);
      centroid += disp;
    }
    centroid = centroid  * (1.0/ (PathData.Path.NumTimeSlices()-1));
    centroid = centroid + PathData.Path(0,ptcl);
    CentroidPos(ptcl)=centroid;
  }  
}


///centroidPos must be resized to the number of particles
///when this function is called

///not sure this gets the right centroid in parallel
void 
TimeHexaticClass::CalculateCentroid_parallel()
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




void 
TimeHexaticClass::BuildDistanceTable(int time1, int time2)
{
  for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
      dVec pos1=PosArray(time1,ptcl1);
      dVec pos2=PosArray(time2,ptcl2);
      dVec disp = pos1-pos2;
      PathData.Path.PutInBox(disp);
      double dist = sqrt(dot(disp,disp));
      DistanceTable(ptcl1,ptcl2)=dist;
    }
  }
}


pair<complex<double> , double> 
TimeHexaticClass::CalcValue(int time1, int time2,int ptcl1,int ptcl2)
{
  complex<double> total=0.0;
  complex<double> val1=OPArray(time1,ptcl1);
  complex<double> val2=OPArray(time2,ptcl2);
  total+=val1*conj(val2);
  double dist = DistanceTable(ptcl1,ptcl2);
  pair<complex<double> , double> temp(total,dist);
  return temp;
}


void
TimeHexaticClass::Accumulate()
{
  if (!FullyWrapped){
    TotalCurrentData++;
    if (TotalCurrentData>=NumStoredMCSteps)
      FullyWrapped=true;
  }
  CalculateCentroid_parallel();
  for (int i=0;i<CentroidPos.size();i++)
    for (int dim=0;dim<NDIM;dim++)
      CentroidPos_write(i,dim)=CentroidPos(i)[dim];
  CentroidPosVar.Write(CentroidPos_write);
  return;
  CurrTime=(CurrTime+1) % NumStoredMCSteps;
  ProduceTimeMatrix(CurrTime);
  int maxToGo=min(NumStoredMCSteps,TotalCurrentData);
  int numParticles=PathData.Path.NumParticles();
  pair<complex<double> ,double> histVals;
  for (int i=0;i<maxToGo;i+=TimeResolution){
    BuildDistanceTable(CurrTime,(CurrTime-i +NumStoredMCSteps) % NumStoredMCSteps);
    for (int ptcl1=0;ptcl1<numParticles;ptcl1++)
      for (int ptcl2=0;ptcl2<numParticles;ptcl2++){
	histVals=CalcValue(CurrTime,(CurrTime-i +NumStoredMCSteps) % NumStoredMCSteps,ptcl1,ptcl2);
	double dist = histVals.second;
	complex<double> val = histVals.first;
	if (dist<grid.End && dist>1.0){
	  int index=grid.ReverseMap(dist);      
	  if (i/TimeResolution>=TimeDisp.extent(1) || 
	      index>=TimeDisp.extent(0))
	    cerr<<"ERROR! "<<endl;
	  assert(i/TimeResolution<TimeDisp.extent(1));
	  assert(index<TimeDisp.extent(0));
	  TimeDisp(index,i/TimeResolution)+=val;
	  gofrDisp(index,i/TimeResolution)+=1;
	}
      }
    assert(i/TimeResolution<NumStepArray.size());
    NumStepArray(i/TimeResolution)=NumStepArray(i/TimeResolution)+1.0;
  }
}


void
TimeHexaticClass::WriteBlock()
{
  return; 
//   for (int i=0;i<TimeDisp.size();i++)
//     if (NumStepArray(i)!=0)
//       TimeDisp(i)=TimeDisp(i)/NumStepArray(i);
  for (int i=0;i<TimeDisp.extent(0);i++)
    for (int j=0;j<TimeDisp.extent(1);j++)
      TimeDisp_double(i,j)=(TimeDisp(i,j)).real();
  TimeDispVar.Write(TimeDisp_double);
  gofrDispVar.Write(gofrDisp);
  NumStepVar.Write(NumStepArray);
  TimeDisp=0.0;
  gofrDisp=0.0;
  NumStepArray=0.0;
  ////  CurrTime=0;
}

void 
TimeHexaticClass::WriteInfo()
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
TimeHexaticClass::ReadGrid(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  //  int numGridPoints;
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

  in.CloseSection();
}

void
TimeHexaticClass::Read (IOSectionClass &in)
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
  NumStoredMCSteps=TimeWindow; ///TimeResolution;
  OPArray.resize(NumStoredMCSteps,PathData.Path.NumParticles());
  PosArray.resize(NumStoredMCSteps,PathData.Path.NumParticles());
  TimeDisp.resize(numGridPoints,NumStoredMCSteps/TimeResolution);
  gofrDisp.resize(numGridPoints,NumStoredMCSteps/TimeResolution);
  TimeDisp_double.resize(numGridPoints,NumStoredMCSteps/TimeResolution);
  NumStepArray.resize(NumStoredMCSteps/TimeResolution);
  CentroidPos.resize(PathData.Path.NumParticles());
  CentroidPos_write.resize(PathData.Path.NumParticles(),NDIM);
  ParticleOrder.resize(PathData.Path.NumParticles());
  DistanceTable.resize(PathData.Path.NumParticles(),PathData.Path.NumParticles());
}
