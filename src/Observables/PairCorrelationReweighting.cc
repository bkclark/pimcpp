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

#include "PairCorrelationReweighting.h"

////////////////////////////////////////
///Pair Correlation Class           ///
///////////////////////////////////////

void PairCorrelationReweightingClass::Read(IOSectionClass& in)
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
  Correction.resize(numGridPoints-1);
  Correction=0;
  NormWeight=0;
  CorrectWeight=0;
  delta_S = 0.0;
  in.CloseSection();
  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc()==0)
    WriteInfo();
  SampledActions.resize(0);
  ReweightActions.resize(0);
  Array<string,1> sampledList, reweightList;
  assert(in.ReadVar("SampledActions",sampledList));
  for(int a=0; a<sampledList.size(); a++) {
    ActionBaseClass* newAction = PathData.Actions.GetAction(sampledList(a));
    SampledActions.push_back(newAction);
  }
  assert(in.ReadVar("ReweightingActions",reweightList));
  for(int a=0; a<reweightList.size(); a++) {
    ActionBaseClass* newAction = PathData.Actions.GetAction(reweightList(a));
    ReweightActions.push_back(newAction);
  }
  activeParticles.resize(PathData.Path.NumParticles());
  for(int p=0; p<activeParticles.size(); p++)
    activeParticles(p) = p;
}



void PairCorrelationReweightingClass::WriteInfo()
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



void PairCorrelationReweightingClass::WriteBlock()
{
  PathClass &Path= PathData.Path;
  Array<double,1> HistSum(Histogram.size());
  Array<double,1> CorrectSum(Correction.size());
  Array<double,1> FinalHist(Correction.size());
  double norm=0.0;
  int N1 = PathData.Species(Species1).NumParticles;
  int N2 = PathData.Species(Species2).NumParticles;
  if (Species1==Species2) //Normalizes things when species are same
    norm = 0.5*(double)TotalCounts * PathData.Path.TotalNumSlices*
      (double)(N1*(N1-1.0))/PathData.Path.GetVol();
  else
    norm = (double)TotalCounts * PathData.Path.TotalNumSlices*
      (double)(N1*N2)/PathData.Path.GetVol();
    //norm = NormWeight * PathData.Path.TotalNumSlices*
    //norm = NormWeight/PathData.Path.GetVol();

  Path.Communicator.Sum(Histogram, HistSum);
  Path.Communicator.Sum(Correction, CorrectSum);
  Array<double,1> gofrArray(HistSum.size());
  Array<double,1> correctArray(CorrectSum.size());
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
    double mean_delta_S = delta_S/(double)TotalCounts;
    gofrArray(i) = HistSum(i) / (binVol*norm);
    correctArray(i) = CorrectSum(i) / (binVol*norm);
    FinalHist(i) = gofrArray(i) * (1 + mean_delta_S) - correctArray(i);
  }
  //gofrVar.Write(gofrArray);
  gofrVar.Write(FinalHist);
  gofrVar.Flush();
  Histogram = 0;
  Correction = 0;
  NormWeight = 0;
  TotalCounts = 0;
  delta_S = 0;
}


void PairCorrelationReweightingClass::Print()
{
  cerr << "UHHHH I called PairCorr Print !!?" << endl;
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
void PairCorrelationReweightingClass::Accumulate()
{

  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);
  // obtain weights
  int BisectionLevel = 0;
  double S_sampled = 0.0;
  cerr << "PairCorrReweight slices " << 0 << " " << PathData.NumTimeSlices()-2 << " activeP " <<  activeParticles.size() << endl;
  list<ActionBaseClass*>::iterator actionIter=SampledActions.begin();
  while (actionIter!=SampledActions.end()){
    S_sampled += ((*actionIter)->Action(0, PathData.NumTimeSlices()-2, activeParticles, BisectionLevel));
    actionIter++;
  }
  double S_reweight = 0.0;
  actionIter=ReweightActions.begin();
  while (actionIter!=ReweightActions.end()){
    S_reweight += ((*actionIter)->Action(0, PathData.NumTimeSlices()-2, activeParticles, BisectionLevel));
    actionIter++;
  }
  double tmp_delta_S = S_reweight - S_sampled;
  cerr << "  S_sampled " << S_sampled << " S_reweight " << S_reweight << " delta_S " << tmp_delta_S << endl;

  //NormWeight += reweight;
  TotalCounts++;
  delta_S += tmp_delta_S;
  if (Species1==Species2) {
    /// Note:  we make sure we don't count that last times slice
    /// we have.  This prevents double counting "shared" slices.
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++)  {
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++) {
        for (int ptcl2=ptcl1+1;ptcl2<=species1.LastPtcl;ptcl2++) {
          dVec disp;
          double dist;
          PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
          if (dist<grid.End) {
            int index=grid.ReverseMap(dist);
            Histogram(index) += 1;
            Correction(index) += tmp_delta_S;
          } 
        }
      }
    }
  }
  else {
    /// Note:  we make sure we don't count that last times slice
    /// we have.  This prevents double counting "shared" slices.
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
      for (int ptcl1=species1.FirstPtcl;ptcl1<=species1.LastPtcl;ptcl1++) {
        for (int ptcl2=species2.FirstPtcl;ptcl2<=species2.LastPtcl;ptcl2++){
          dVec disp;
          double dist;
          PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
          if (dist<grid.End) {
            int index=grid.ReverseMap(dist);
            Histogram(index) += 1;
            Correction(index) += tmp_delta_S;
          }
        }
      }
    }
  }
}


void PairCorrelationReweightingClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
}
