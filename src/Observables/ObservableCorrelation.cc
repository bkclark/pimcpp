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

#include "ObservableCorrelation.h"

////////////////////////////////////////
///Pair Correlation Class           ///
///////////////////////////////////////

void PairCorrelationClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  string species1Name;
  string species2Name;
  Species1=-1;
  Species2=-1;
  assert(in.ReadVar("Species1",species1Name));
  assert(in.ReadVar("Species2",species2Name));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==species2Name){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  assert(gridType == "Linear");
  bool readStartGrid = in.ReadVar("start",gridStart);
  bool readEndGrid = in.ReadVar("end",gridEnd);
  if (!readStartGrid)
    gridStart = 0.0;
  if (!readEndGrid) {
    if (PathData.Path.GetPeriodic()[0])
      gridEnd = PathData.Path.GetBox()[0]/2.;
    else {
      cerr<<"I don't know where you want me to end this grid"<<endl;
      assert(1==2);
    }
  }
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);

  // Measure Ref Slice by itself
  if(!in.ReadVar("RefSliceOnly", RefSliceOnly))
    RefSliceOnly = 0;

  HaveRefSlice = ((PathData.Path.Species(Species1).GetParticleType() == FERMION &&
                   PathData.Actions.NodalActions(Species1) != NULL &&
                   !PathData.Actions.NodalActions(Species1)->IsGroundState()) ||
                  (PathData.Path.Species(Species2).GetParticleType() == FERMION &&
                   PathData.Actions.NodalActions(Species2) != NULL &&
                   !PathData.Actions.NodalActions(Species2)->IsGroundState()));

  // Record g(r) for each time slice
  if(!in.ReadVar("RefSliceOnly", RefSliceOnly))
    RefSliceOnly = 0;

  // Record g(r) for each time slice
  if(!in.ReadVar("PerSlice", PerSlice))
    PerSlice = 0;

  // Sign Tracking
  if(!in.ReadVar("TrackSign", TrackSign))
    TrackSign = 0;

  Histogram.resize(numGridPoints-1);
  Histogram = 0;
  if (RefSliceOnly) {
    HistogramRefSlice.resize(numGridPoints-1);
    HistogramRefSlice = 0;
  }
  if (PerSlice) {
    HistogramPerSlice.resize(PathData.Path.TotalNumSlices,numGridPoints-1);
    HistogramPerSlice = 0;
  }

  in.CloseSection();

  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();

  TotalCounts=0;
}


void PairCorrelationClass::WriteInfo()
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
    r(i) = 0.75 * (rb*rb*rb*rb-ra*ra*ra*ra)/(rb*rb*rb-ra*ra*ra);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "g(r)");
  IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
  IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
  IOSection.WriteVar("Type","CorrelationFunction");
  IOSection.WriteVar("Cumulative", false);
}


void PairCorrelationClass::WriteBlock()
{
  PathClass &Path= PathData.Path;
  int M = PathData.Path.TotalNumSlices;
  double V;
  dVec IsPeriodic = PathData.Path.GetPeriodic();
  if (IsPeriodic(0)) // HACK: only checks 1 dimension
    V = PathData.Path.GetVol();
  else
    V = 1.;
  int N1 = PathData.Species(Species1).NumParticles;
  int N2 = PathData.Species(Species2).NumParticles;

  // Sum up histograms
  Array<int,1> HistSum(Histogram.size());
  Array<int,1> HistSumRefSlice(HistogramRefSlice.size());
  Array<int,2> HistSumPerSlice(M,Histogram.size());
  Path.Communicator.Sum(Histogram, HistSum);
  if (RefSliceOnly)
    Path.Communicator.Sum(HistogramRefSlice, HistSumRefSlice);
  if (PerSlice)
    Path.Communicator.Sum(HistogramPerSlice, HistSumPerSlice);

  // Normalize histograms
  Array<double,1> gofrArray(HistSum.size());
  Array<double,1> gofrArrayRefSlice(HistSumRefSlice.size());
  Array<double,2> gofrArrayPerSlice(M,HistSum.size());
  double norm = 0.0;
  if (Species1==Species2) //Normalizes things when species are same
    norm = 0.5*(double)TotalCounts * (double)(N1*(N1-1.0))/V;
  else
    norm = (double)TotalCounts * (double)(N1*N2)/V;
  for (int i=0; i<grid.NumPoints-1; i++){
    double r1 = grid(i);
    double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
    double r = 0.5*(r1+r2);
#if NDIM==3
    double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
#endif
#if NDIM==2
    double binVol = M_PI * (r2*r2-r1*r1);
#endif
    gofrArray(i) = (double) HistSum(i) / (binVol*norm*M);
    if (RefSliceOnly)
      gofrArrayRefSlice(i) = (double) HistSumRefSlice(i) / (binVol*norm);
    if (PerSlice) {
      for (int j = 0; j < M; j++)
        gofrArrayPerSlice(j,i) = (double) HistSumPerSlice(j,i) / (binVol*norm*M);
    }
  }

  // Write out arrays
  gofrVar.Write(gofrArray);
  gofrVar.Flush();
  Histogram = 0;
  if (RefSliceOnly) {
    gofrVarRefSlice.Write(gofrArrayRefSlice);
    gofrVarRefSlice.Flush();
    HistogramRefSlice = 0;
  }
  if (PerSlice) {
    gofrVarPerSlice.Write(gofrArrayPerSlice);
    gofrVarPerSlice.Flush();
    HistogramPerSlice = 0;
  }
  TotalCounts = 0;
}


void PairCorrelationClass::Print()
{
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      double gofr = (double) Histogram(i) / (vol*TotalCounts);
      fprintf (stderr, "%1.12e %1.12e\n", r, gofr);
    }
}


/// Fix me to accumulate data only between the two species I'm
/// interested in.
void PairCorrelationClass::Accumulate()
{

  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);

  TotalCounts++;

  double FullWeight = CalcFullWeight();
  if (!TrackSign)
    FullWeight = 1;

  int localRef = -1;
  if (HaveRefSlice && RefSliceOnly) {
    int myProc = PathData.Path.Communicator.MyProc();
    int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
    if (procWithRefSlice == myProc) {
      /// Note:  Pair Correlation only defined on reference slice
      int firstSlice, lastSlice;
      Path.SliceRange (myProc, firstSlice, lastSlice);
      localRef = Path.GetRefSlice() - firstSlice;
    }
  }

  /// Note:  we make sure we don't count that last times slice
  /// we have.  This prevents double counting "shared" slices.
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++) {
      int ptcl2a = species2.FirstPtcl;
      int ptcl2b = species2.LastPtcl;
      if (Species1==Species2) {
        ptcl2a = ptcl1+1;
        ptcl2b = species1.LastPtcl;
      }
      for (int ptcl2=ptcl2a;ptcl2<=ptcl2b;ptcl2++) {
        dVec disp;
        double dist;
        PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
        if (dist<grid.End) {
          int index=grid.ReverseMap(dist);
          Histogram(index) += FullWeight*1;
          if (RefSliceOnly && slice == localRef)
            HistogramRefSlice(index) += FullWeight*1;
        }
      }
    }
  }

  /// Note:  we make sure we don't count that last times slice
  /// we have.  This prevents double counting "shared" slices.
  if (PerSlice) {
    for (int sliceDiff=0; sliceDiff<PathData.Path.NumTimeSlices()-1; sliceDiff++) {
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++) {
        int ptcl2a = species2.FirstPtcl;
        int ptcl2b = species2.LastPtcl;
        if (Species1==Species2) {
          ptcl2a = ptcl1+1;
          ptcl2b = species1.LastPtcl;
        }
        for (int ptcl2=ptcl2a;ptcl2<=ptcl2b;ptcl2++) {
          dVec disp;
          double dist;
          for (int slice1=0; slice1<PathData.Path.NumTimeSlices()-1; slice1++) {
            int slice2 = (slice1 + sliceDiff) % PathData.Path.TotalNumSlices;
            PathData.Path.DistDisp(slice1,slice2,ptcl1,ptcl2,dist,disp);
            if (dist<grid.End) {
              int index=grid.ReverseMap(dist);
              HistogramPerSlice(sliceDiff,index) += FullWeight*1;
            }
          }
        }
      }
    }
  }
}


void PairCorrelationClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled = 0;
}


////////////////////////////////////////
///N(r) Correlation Class            ///
///////////////////////////////////////

///Initializes the class. The grid should always be overwritten with
///some other size. 
void nofrClass::Initialize()
{
  TotalCounts=0;
  TimesCalled=0;
}


///Reads in this classes input.
void nofrClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);


  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  // Hack (Forced Linear)
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
  TotalCounts=0;
  Histogram.resize(numGridPoints-1);
  Histogram=0;
  in.CloseSection();

}

///Writes the data relevant for this classes output including its
///name, axis to plot, etc.
void nofrClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  IOSection.NewSection("grid");
  grid.Write(IOSection);
  IOSection.CloseSection();

  int numBins = grid.NumPoints-1;
  Array<double,1> r(numBins);
  Array<double,1> r2(numBins);
  for (int i=0; i<numBins; i++) {
    double ra = grid(i);
    double rb = grid(i+1);
    //r(i) = 0.75 * (rb*rb*rb*rb-ra*ra*ra*ra)/(rb*rb*rb-ra*ra*ra);
    r(i) = 0.5*(ra+rb);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "n(r)");
  IOSection.WriteVar("Type","CorrelationFunction");
  IOSection.WriteVar("Cumulative", false);
}

///Writes a block of histogram data.  Does not compensate for volume
///effects or for importance sampling. This has to be done after the
///fact currently
void nofrClass::WriteBlock()
{

  PathClass &Path= PathData.Path;
  Array<double,1> HistSum(Histogram.size());
  double norm=0.0;
  norm = TotalCounts/PathData.Path.GetVol();
  Path.Communicator.Sum(Histogram, HistSum);
  ///This will only work in serial because I'm not summing!!!!
  //\\  Histogram3d=Histogram3d/norm;
  //\\  HistSum3d.Write(Histogram3d);
  if (Path.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
    Array<double,1> nofrArray(HistSum.size());
    for (int i=0; i<grid.NumPoints-1; i++){
      double r1 = grid(i);
      double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
      double r = 0.5*(r1+r2);
#if NDIM==3
      double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
#endif
#if NDIM==2
      double binVol = M_PI * (r2*r2-r1*r1);
#endif
      ////////////////)/////////
      // This line does not normalize by volume for a dimer, e.g. -jg
      // nofrArray(i) = (double) HistSum(i) / (norm);
      //////////////////////////
      nofrArray(i) = (double) HistSum(i) / (binVol*norm);
      if (isnan(nofrArray(i)) || isinf(nofrArray(i)))
        nofrArray(i) = 0.0;
    }
    nofrVar.Write(nofrArray);
    nofrVar.Flush();
  }
  TotalCounts=0;
  Histogram=0;
}


void nofrClass::Print()
{
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
#if NDIM==3
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
#endif
#if NDIM==2
      double vol = M_PI * (r2*r2-r1*r1);
#endif
      double nofr = (double) Histogram(i) / (vol*TotalCounts);
      fprintf (stderr, "%1.12e %1.12e\n", r, nofr);
    }
}


///Bins the distance between the head and the tail of the open
///loop. This does not compensaite for volume effects or importance sampling
void nofrClass::Accumulate()
{
  
  int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);  
  if (procWithRefSlice == PathData.Path.Communicator.MyProc()) {

    dVec disp=0.0;
    double dist2;
    int openLink=(int)(PathData.Path.OpenLink);
    int openPtcl=(int)(PathData.Path.OpenPtcl);
    PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
                           dist2,disp); //This is distance between head and tail!
    //cerr<<"DIST IS "<<dist2<<endl;
    int numLinks=PathData.Path.NumTimeSlices()-1;
    disp=0.0;
    int currSlice=openLink;
    int currPtcl=openPtcl;
    int nextSlice=-1;
    int nextPtcl=-1;

    while (nextSlice!=openLink || nextPtcl!=openPtcl){
      nextSlice = (currSlice + 1) % PathData.Path.NumTimeSlices();
      //    if (nextSlice==0)
      //      nextSlice=numLinks+1;
      if (currSlice==PathData.Join)
        nextPtcl=PathData.Path.Permutation(currPtcl);
      else 
        nextPtcl=currPtcl;
      dVec linkDisp=PathData.Path.VelocityBetweenPtcl(currSlice,currPtcl,nextSlice,nextPtcl);
      disp=disp+linkDisp;
      currSlice=nextSlice;
      currPtcl=nextPtcl;
    }
    double dist=sqrt(dot(disp,disp));
    //  if (((dist-dist2)/PathData.Path.GetBox()[0])-floor((dist-dist2)/PathData.Path.GetBox()[0]+0.1)<1e-5)
    //    cerr<<"dist, dist2, diff: "<<dist<<" "<<dist2<<" "<<dist-dist2<<endl;

    //  for (int slice=0;slice<numLinks;slice++) {
    //    int realSlice=(openLink+slice) % numLinks;
    //    int realSlicep1=(openLink+slice+1) % numLinks;
    //    dVec linkDisp;
    //    linkDisp=PathData.Path.Velocity(realSlice,realSlicep1,openPtcl);
    //    disp =disp+linkDisp;
    //  }
    //  double  dist=sqrt(dot(disp,disp));
    //  if (dist-dist2>1e-5)
    //cerr<<"dist, dist2, diff: "<<dist<<" "<<dist2<<" "<<dist-dist2<<endl;
    if (dist<grid.End){
      int index=grid.ReverseMap(dist);
      //      if (PathData.Actions.OpenLoopImportance.ImpChoice==DISPXIMP){
      //        int myProc=PathData.GetCloneNum();
      //        double shift=(myProc % 16)+0.5;
      //        Histogram(index)=Histogram(index)+exp(-(disp(0)-shift)*(disp(0)-shift));
      //      }
      //      else {
        //    Histogram(index)=Histogram(index)+(0.5)/(dist*dist)+(0.9*exp(-dist*dist)+0.1);
        Histogram(index)=Histogram(index)+1.0;
        //   }
    //      if (disp(0)<grid.End && disp(1)<grid.End && disp(2)<grid.End){
    //  int index0=grid.ReverseMap(disp(0));
    //  int index1=grid.ReverseMap(disp(1));
    //  int index2=grid.ReverseMap(disp(2));
        //\\      Histogram3d(index0,index1,index2)=Histogram3d(index0,index1,index2)+1.0;
        //      }
    }
  }

  TotalCounts++;
  //  cerr<<"done accumulating"<<endl;
  return; 
}
