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

#include "ObservableRefCorrelation.h"

////////////////////////////////////////
///Pair Correlation Class           ///
///////////////////////////////////////

void RefPairCorrelationClass::Read(IOSectionClass& in)
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

  HaveRefSlice = ((PathData.Path.Species(Species1).GetParticleType() == FERMION &&
                   PathData.Actions.NodalActions(Species1) != NULL &&
                   !PathData.Actions.NodalActions(Species1)->IsGroundState()) ||
                  (PathData.Path.Species(Species2).GetParticleType() == FERMION &&
                   PathData.Actions.NodalActions(Species2) != NULL &&
                   !PathData.Actions.NodalActions(Species2)->IsGroundState()));

  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();
}



void RefPairCorrelationClass::WriteInfo()
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
  IOSection.WriteVar("HaveRefSlice", HaveRefSlice);
}



void RefPairCorrelationClass::WriteBlock()
{
  PathClass &Path= PathData.Path;
  Array<int,1> HistSum(Histogram.size());
  double norm=0.0;
  int N1 = PathData.Species(Species1).NumParticles;
  int N2 = PathData.Species(Species2).NumParticles;

  if (HaveRefSlice) {
    if (Species1==Species2) //Normalizes things when species are same
      norm = 0.5*(double)TotalCounts * (double)(N1*(N1-1.0))/PathData.Path.GetVol();
    else
      norm = (double)TotalCounts * (double)(N1*N2)/PathData.Path.GetVol();
  }

  Path.Communicator.Sum(Histogram, HistSum);
  Array<double,1> gofrArray(HistSum.size());
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
                //////////////////////////
                // This line does not normalize by volume for a dimer, e.g. -jg
    //gofrArray(i) = (double) HistSum(i) / (norm);
                //////////////////////////
    gofrArray(i) = (double) HistSum(i) / (binVol*norm);
  }
  gofrVar.Write(gofrArray);
  gofrVar.Flush();
  Histogram = 0;
  TotalCounts = 0;
}


void RefPairCorrelationClass::Print()
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
void RefPairCorrelationClass::Accumulate()
{

  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);

  TotalCounts++;
  if (HaveRefSlice) {
    int myProc = PathData.Path.Communicator.MyProc();
    int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
    if (procWithRefSlice == myProc) {
      /// Note:  Pair Correlation only defined on reference slice
      int firstSlice, lastSlice;
      Path.SliceRange (myProc, firstSlice, lastSlice);
      int localRef = Path.GetRefSlice() - firstSlice;
      if (Species1==Species2) {
        for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++)
          for (int ptcl2=ptcl1+1;ptcl2<=species1.LastPtcl;ptcl2++) {
            dVec disp;
            double dist;
            PathData.Path.DistDisp(localRef,ptcl1,ptcl2,dist,disp);
            if (dist<grid.End) {
              int index=grid.ReverseMap(dist);
              Histogram(index)++;
            }
          }
      } else {
        for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++) {
          for (int ptcl2=species2.FirstPtcl;ptcl2<=species2.LastPtcl;ptcl2++) {
            dVec disp;
            double dist;
            PathData.Path.DistDisp(localRef,ptcl1,ptcl2,dist,disp);
            if (dist<grid.End) {
              int index=grid.ReverseMap(dist);
              Histogram(index)++;
            }
          }
        }
      }
    }
  }
}


void RefPairCorrelationClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
}

