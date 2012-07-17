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

#include "Hexatic.h"
#include <complex>

double unit2angle(double x,double y)
{
  double angle=atan(y/x);
  if (x<0)
    angle=angle+M_PI;
  return angle;
}

complex<double>
HexaticClass::OrderParamater(Array<dVec,1> centroidPos,int ptcl)
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
	if (!((dot(r12disp,r12disp)-1.0)<0.001))
	  cerr<<"GRRRL: "<<dot(r12disp,r12disp)<<" "<<r12dist<<" "<<ptcl<<" "<<nearPtcl<<endl;
	///	cerr<<"hmm: "<<dot(r12disp,r12disp)<<endl;
	assert(abs(dot(r12disp,r12disp)-1.0)<0.001);
	double theta_12=unit2angle(r12disp(0),r12disp(1));
	op=op+complex<double>(cos(theta_12*q),sin(theta_12*q));
      }
    }
  }
  return op;
}



complex<double>
HexaticClass::OrderParamater(int slice,int ptcl)
{
  complex<double> op=0.0;
  for (int nearPtcl=0;nearPtcl<PathData.Path.NumParticles();
       nearPtcl++){
    double r12dist;
    dVec r12disp;
    if (nearPtcl!=ptcl){
      PathData.Path.DistDisp(slice,ptcl,nearPtcl,r12dist,r12disp);
      if (r12dist<DistCutoff){
	r12disp=r12disp *(1.0/r12dist);
	if (abs(dot(r12disp,r12disp)-1.0)>=0.001)
	  cerr<<dot(r12disp,r12disp);
	if (!((dot(r12disp,r12disp)-1.0)<0.001))
	  cerr<<"GRRR2: "<<dot(r12disp,r12disp)<<endl;
	//	cerr<<"hmm2: "<<dot(r12disp,r12disp)<<endl;
	assert(abs(dot(r12disp,r12disp)-1.0)<0.001);
	double theta_12=unit2angle(r12disp(0),r12disp(1));
	op=op+complex<double>(cos(theta_12*q),sin(theta_12*q));
      }
    }
  }
  return op;
}

complex<double>
Conj21(complex<double> a)
{
  complex<double> b(a.real(),-a.imag());
  return b;
}


void 
HexaticClass::CalculateCentroid_parallel()
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
    centroid = centroid * (1.0/(PathData.Path.TotalNumSlices));
    centroid = centroid + zeroVec;
    CentroidPos(ptcl)=centroid;
  }  
}


void
HexaticClass::Accumulate_old()
{
  
  //  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++){
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      ParticleOrder(ptcl)=OrderParamater(slice,ptcl);

    }
    for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
      for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
	double dist;
	dVec disp;
	PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
	if (dist<grid.End){
	  int index=grid.ReverseMap(dist);
//           if (index==0){
//             cerr<<"GRRR: "<<slice<<" "<<ptcl1<<" "<<ptcl2<<" "<<index<<" "
//                 <<Histogram(index)<<" "<<ParticleOrder(ptcl1)<<" "
//                 <<ParticleOrder(ptcl1)<<" "<<ParticleOrder(ptcl2)<<" "
//                 <<Conj21(ParticleOrder(ptcl2))
//                 <<" "<<ParticleOrder(ptcl1)*conj(ParticleOrder(ptcl2))<<endl;
//           }
	  Histogram(index)=Histogram(index)+
	    ParticleOrder(ptcl1)*conj(ParticleOrder(ptcl2));
// 	  if (index==0){
//             cerr<<"GRRR: "<<slice<<" "<<ptcl1<<" "<<ptcl2<<" "
//                 <<Histogram(index)<<" "<<ParticleOrder(ptcl1)<<" "
//                 <<ParticleOrder(ptcl1)<<" "<<ParticleOrder(ptcl2)<<" "
// 		<<Conj21(ParticleOrder(ptcl2))
//                 <<" "<<ParticleOrder(ptcl1)*conj(ParticleOrder(ptcl2))<<endl;
// 	  }

	}
      }
    }
  }
//   cerr<<"proc"<<PathData.Path.Communicator.MyProc()<<"HERE"<<endl;
//   for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
//     cerr<<"proc"<<PathData.Path.Communicator.MyProc()<<"Ptcl "<<ParticleOrder(ptcl)<<" "<<endl;
//   }
//   for (int i=0;i<Histogram.size();i++){
//     cerr<<"proc"<<PathData.Path.Communicator.MyProc()<<" "<<Histogram(i)<<" "<<endl;
//   }
  NumSamples++;
}

void
HexaticClass::Accumulate()
{
  if (!Centroid){
    Accumulate_old();
    return;
  }
  CalculateCentroid_parallel();
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    ParticleOrder(ptcl)=OrderParamater(CentroidPos,ptcl);
  }
  for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
      dVec disp=CentroidPos(ptcl2)-CentroidPos(ptcl1);
      PathData.Path.PutInBox(disp);
      double dist=sqrt(dot(disp,disp));
      if (dist<grid.End){
	int index=grid.ReverseMap(dist);
	Histogram(index)=Histogram(index)+
	  ParticleOrder(ptcl1)*conj(ParticleOrder(ptcl2));
      }
    }
  }
  
  NumSamples++;
}


void
HexaticClass::WriteBlock()
{
  cerr<<"Hexatic class writing"<<endl;
  PathClass &Path= PathData.Path;
  double norm=1.0/((double)NumSamples*Path.TotalNumSlices);
  for (int counter=0;counter<Histogram.size();counter++){
  //    cerr<<"Proc"<<PathData.Path.Communicator.MyProc()<<" "<<counter<<Histogram(counter).real()<<endl;
    HistDouble(counter)=Histogram(counter).real();
  }
  Path.Communicator.Sum(HistDouble,HistSum);
  HistSum=HistSum*norm;
  HexaticRealVar.Write(HistSum);

  for (int counter=0;counter<Histogram.size();counter++)
    HistDouble(counter)=Histogram(counter).imag();
  Path.Communicator.Sum(HistDouble,HistSum);
  HistSum=HistSum*norm;
  HexaticImagVar.Write(HistSum);

  HexaticRealVar.Flush();
  Histogram=0.0;
  HistSum=0.0;
  HistDouble=0.0;
  NumSamples=0;
  cerr<<"Hexatic class writing done"<<endl;
}

void 
HexaticClass::WriteInfo()
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
HexaticClass::ReadGrid(IOSectionClass &in)
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
HexaticClass::Read (IOSectionClass &in)
{
  complex<double> test(1.0,2.0);
  cerr<<"CONJ: "<<Conj21(test)<<endl;
  NumSamples=0;
  ParticleOrder.resize(PathData.Path.NumParticles());
  CentroidPos.resize(PathData.Path.NumParticles());
  ///It's probably important that the grid is the same grid that is in
  ///the pair correlation function. Not sure how to authenticate this.
  ReadGrid(in);
  if (!in.ReadVar("Centroid",Centroid))
      Centroid=false;
  ObservableClass::Read(in);
  cerr<<"Centroid is "<<Centroid<<endl;
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();

}
