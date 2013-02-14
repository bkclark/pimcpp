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
  int N = PathData.Path.NumParticles();

  // Loop through particles, making sure to visit each only once
  int initSlice = 0;
  bool usedPtcl[N];
  Array<TinyVector<double,3>,1> tmpCentPos(N);
  for (int ptcl = 0; ptcl < N; ptcl++) {
    usedPtcl[ptcl] = false;
    tmpCentPos(ptcl) = 0.0;
  }
  for (int ptcl = 0; ptcl < N; ptcl++) {
    if (!usedPtcl[ptcl]) { // Only visit each particle once

      // Calculate centroid position
      int currSlice = initSlice;
      int currPtcl = ptcl;
      int nextSlice = -1;
      int nextPtcl = -1;
      dVec centPos = 0.0;
      int numSlices = 0;
      int numPtcls = 0;
      while (nextSlice != initSlice || nextPtcl != ptcl) {
        nextSlice = (currSlice + 1) % (PathData.Path.NumTimeSlices()-1);
        if (currSlice == PathData.Join) {
          nextPtcl = PathData.Path.Permutation(currPtcl);
          numPtcls++;
        } else
          nextPtcl = currPtcl;
        usedPtcl[nextPtcl] = true;
        dVec slicePos = PathData.Path(currSlice,currPtcl);
        //PathData.Path.PutInBox(slicePos); // HACK: DO I NEED THIS???
        centPos = centPos + slicePos;
        currSlice = nextSlice;
        currPtcl = nextPtcl;
        numSlices++;
      }
      tmpCentPos(ptcl) = centPos;
    }
  }

  // Gather up all the position data to get the centroids
  Array<TinyVector<double,3>,1> totCentPos(N);
  PathData.Path.Communicator.AllSum(tmpCentPos,totCentPos);

  // Calculate spread from centroid
  for (int ptcl = 0; ptcl < N; ptcl++)
    usedPtcl[ptcl] = false;
  for (int ptcl = 0; ptcl < N; ptcl++) {
    if (!usedPtcl[ptcl]) { // Only visit each particle once
      dVec centPos = totCentPos(ptcl)/PathData.Path.TotalNumSlices;
      int currSlice = initSlice;
      int currPtcl = ptcl;
      int nextSlice = -1;
      int nextPtcl = -1;
      while (nextSlice != initSlice || nextPtcl != ptcl) {
        nextSlice = (currSlice + 1) % (PathData.Path.NumTimeSlices()-1);
        if (currSlice == PathData.Join)
          nextPtcl = PathData.Path.Permutation(currPtcl);
        else
          nextPtcl = currPtcl;
        usedPtcl[nextPtcl] = true;
        dVec centDisp;
        double centDist;
        PathData.Path.DistDispPos(currSlice,currPtcl,centPos,centDist,centDisp);

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


