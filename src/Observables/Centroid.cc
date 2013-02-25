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

#include "Centroid.h"


/// Writes the data relevant for this classes output including
/// its name, axis to plot, etc.
void CentroidClass::WriteInfo()
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
    r(i) = 0.5*(ra+rb);
  }
  IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "r");
  IOSection.WriteVar("ylabel", "n(r)");
  IOSection.WriteVar("Type","CorrelationFunction");
  IOSection.WriteVar("Cumulative", false);
}


/// Writes a block of histogram data.
void CentroidClass::WriteBlock()
{
  PathClass &Path = PathData.Path;
  Array<int,1> centHistSum(centHistogram.size());
  Path.Communicator.Sum(centHistogram, centHistSum);
  double norm = NumSamples*Path.NumParticles()*Path.TotalNumSlices/Path.GetVol();
  if (Path.Communicator.MyProc() == 0) {
    Array<double,1> centArray(centHistSum.size());
    for (int i=0; i<grid.NumPoints-1; i++) {
      double r1 = grid(i);
      double r2 = (i<(grid.NumPoints-1)) ? grid(i+1):(2.0*grid(i)-grid(i-1));
      double r = 0.5*(r1+r2);
#if NDIM==3
      double binVol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
#endif
#if NDIM==2
      double binVol = M_PI * (r2*r2-r1*r1);
#endif
      centArray(i) = (double) centHistSum(i) / (binVol*norm);
      if (isnan(centArray(i)) || isinf(centArray(i)))
        centArray(i) = 0.0;
    }
    CentroidSpreadVar.Write(centArray);
    CentroidSpreadVar.Flush();
  }
  NumSamples = 0;
  centHistogram = 0;
}


void CentroidClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);

  // Setup the grid and histograms
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  assert(gridType == "Linear"); // HACK: Forced linear grid
  bool readStartGrid = in.ReadVar("start",gridStart);
  bool readEndGrid = in.ReadVar("end",gridEnd);
  if (!readStartGrid)
    gridStart=0.0;
  if (!readEndGrid) {
    if (PathData.Path.GetPeriodic()[0]) {
      gridEnd=PathData.Path.GetBox()[0];
    } else {
      cerr << PathData.Path.CloneStr << " I don't know where you want me to end this centroid grid." << endl;
      abort();
    }
  }
  assert(in.ReadVar("NumPoints",numGridPoints));
  grid.Init(gridStart,gridEnd,numGridPoints);
  centHistogram.resize(numGridPoints-1);
  centHistogram = 0;
  in.CloseSection();

  /// Now write the one-time output variables
  if (PathData.Path.Communicator.MyProc() == 0) {
    WriteInfo();
  }

}

// Check out MetaMoves.cc! Shift was changed to 0 b/c for some reason it's messing this procedure up!
void CentroidClass::Accumulate()
{
  NumSamples++;
  const int N = PathData.Path.NumParticles();

  // Get the centroid positions
  Array<TinyVector<double,NDIM>,1> CentPos(N);
  PathData.GetCentroids(CentPos);

  // Calculate the variance in each direction for each particle
  const int M = PathData.Path.NumTimeSlices()-1;
  const int totM = PathData.Path.TotalNumSlices;
  for (int ptcl = 0; ptcl < N; ptcl++) {

    // Build the covariance matrix
    const int D = NDIM;
    Array<double,2> C(D,D); // Will be the covariance matrix.
    C = 0;                  // Eigenvectors will be major & minor axes. Eigenvalues will be lengths of axes, squared.
    for (int k = 0; k < M; k++) {
      dVec diff = PathData.Path(k,ptcl) - CentPos(ptcl);
      PathData.Path.PutInBox(diff);
      for (int d1 = 0; d1 < D; d1++)
        for (int d2 = 0; d2 < D; d2++)
          C(d1,d2) += diff(d1)*diff(d2);
    }

    // Gather all covariance matrix totals and normalize
    Array<double,2> totC(D,D);
    PathData.Path.Communicator.AllSum(C,totC);
    totC = totC/(totM-1);

    // Diagonalize the covariance matrix
    Array<double,1> Vals;
    Array<double,2> Vecs;
    SymmEigenPairs(totC,D,Vals,Vecs);
    //cout << CentPos(ptcl) << totC << D << Vals << Vecs << endl;
  }

  // Calculate spread from centroid (PROBABLY GOING TO END UP REMOVING THIS!!!)
  bool usedPtcl[N];
  for (int ptcl = 0; ptcl < N; ptcl++)
    usedPtcl[ptcl] = false;
  for (int ptcl = 0; ptcl < N; ptcl++) {
    if (!usedPtcl[ptcl]) { // Only visit each particle once
      int currSlice = 0;
      int currPtcl = ptcl;
      int nextSlice = -1;
      int nextPtcl = -1;
      while (nextSlice != 0 || nextPtcl != ptcl) {
        nextSlice = (currSlice + 1) % M;
        if (currSlice == PathData.Join)
          nextPtcl = PathData.Path.Permutation(currPtcl);
        else
          nextPtcl = currPtcl;
        usedPtcl[nextPtcl] = true;
        dVec centDisp;
        double centDist;
        PathData.Path.DistDispPos(currSlice,currPtcl,CentPos(ptcl),centDist,centDisp);

        // Add distance from centroid to histogram
        if (centDist < grid.End) {
          int index = grid.ReverseMap(centDist);
          centHistogram(index)++;
        }
        currSlice = nextSlice;
        currPtcl = nextPtcl;
      }
    }
  }

}
