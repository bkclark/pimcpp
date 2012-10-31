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

// #include "../MPI/Communication.h"
#include "FreeNodalActionClass.h"
#include "../PathDataClass.h"
#include "../MatrixOps/MatrixOps.h"
#include "ctime"
#include "sys/time.h"

double FreeNodalActionClass::ActionImageSum (double L, double lambdaBeta,  double disp)
{
  int numImages = 10;
  double sum = 0.0;
  double fourLambdaBetaInv = (lambdaBeta!=0.0) ?  1.0/(4.0*lambdaBeta) : 0.0;
  // If the images won't contributed anything, let's not worry
  // about image sums.
  if ((disp*disp*fourLambdaBetaInv) > 50.0)
    return (disp*disp*fourLambdaBetaInv);
  for (int image=-numImages; image<=numImages; image++) {
    double x = disp + (double)image*L;
    sum += exp (-(x*x)*fourLambdaBetaInv);
  }
  return (-log(sum));
}

double FreeNodalActionClass::ActionkSum (double L, double lambdaBeta, double disp)
{
  double kmax = sqrt (50.0/lambdaBeta);
  double kInc = 2.0*M_PI/L;

  if (lambdaBeta == 0.0)
    return (0.0);

  double sum = 0.5;
  for (double k=kInc; k<=kmax; k+=kInc)
    sum += cos(k*disp)*exp(-lambdaBeta*k*k);
  return (-log(sum));
}


void FreeNodalActionClass::Init()
{
  FirstDistTime = 1;
  FirstDetTime = 1;
  // Initialize NodeDist and/or NodeDet
  if (PathData.Path.UseNodeDist||PathData.Path.UseNodeDet) {
    SetMode(NEWMODE);
    SpeciesClass &species = PathData.Path.Species(SpeciesNum);
    double lambda = species.lambda;
    double levelTau = PathData.Path.tau;
    int myStart, myEnd;
    Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
    int refSlice = Path.GetRefSlice() - myStart;
    int startSlice = 0;
    int endSlice = PathData.Path.NumTimeSlices()-1;
    for (int slice=startSlice; slice<endSlice; slice++) {
      int sliceDiff = abs(slice-refSlice);
      sliceDiff = min (sliceDiff, PathData.Path.TotalNumSlices-sliceDiff);
      if (sliceDiff!=0) {
        if (PathData.Path.UseNodeDist)
          PathData.Path.NodeDist(slice,SpeciesNum) = HybridDist(slice, lambda*levelTau);
        if (PathData.Path.UseNodeDet)
          PathData.Path.NodeDet(slice,SpeciesNum) = Det(slice);
      }
    }
    if (PathData.Path.UseNodeDist)
      PathData.Path.NodeDist[OLDMODE](Range(startSlice,endSlice), SpeciesNum) =
        PathData.Path.NodeDist[NEWMODE](Range(startSlice,endSlice), SpeciesNum);
    if (PathData.Path.UseNodeDet)
      PathData.Path.NodeDet[OLDMODE](Range(startSlice,endSlice), SpeciesNum) =
        PathData.Path.NodeDet[NEWMODE](Range(startSlice,endSlice), SpeciesNum);
  }
}

void 
FreeNodalActionClass::SetupFreeActions()
{
  //cerr<<"SETTING UP FREE ACTION"<<endl;
  const int nPoints = 1000;
  // Setup grids
  for (int i=0; i<NDIM; i++)
    ActionGrids[i].Init (-0.5*Path.GetBox()[i], 0.5*Path.GetBox()[i], nPoints);

  Array<double,1> actionData(nPoints);
  int nSplines = Path.TotalNumSlices/2 + (Path.TotalNumSlices%2)+1;
  double lambdaTau = Path.tau * Path.Species(SpeciesNum).lambda;


  /// DEBUG
  //FILE *fout;
  //fout = fopen ("FreeImageActions.dat", "w");
  // Now, setup up actions
  ActionSplines.resize(nSplines);
  for (int spline=0; spline<nSplines; spline++) {
    double lambdaBeta = lambdaTau * (double)spline;
    for (int dim=0; dim<NDIM; dim++) {
      double L = Path.GetBox()[dim];
      for (int i=0; i<nPoints; i++) {
        double disp = ActionGrids[dim](i);
        actionData(i) = ActionImageSum (L, lambdaBeta, disp) - ActionImageSum(L, lambdaBeta, 0.0);
        //fprintf (fout, "%1.12e ", actionData(i));
      }
      //fprintf (fout, "\n");
      // Since the action is periodic, the slope should be zero
      // at the boundaries
      ActionSplines(spline)[dim].Init (&ActionGrids[dim], actionData, 0.0, 0.0);
    }
  }
  //fclose (fout);

}


FreeNodalActionClass::FreeNodalActionClass (PathDataClass &pathData, int speciesNum) :
  NodalActionClass (pathData), Path (pathData.Path), SpeciesNum (speciesNum)
{
  int N = Path.Species(speciesNum).LastPtcl - Path.Species(speciesNum).FirstPtcl+1;
  DetMatrix.resize(N,N);
  Cofactors.resize(N,N);
  GradVec.resize(N);
  SavePath.resize(N);
  NumLineDists = NumGradDists = 0;
  nSingular = 0;
}


double FreeNodalActionClass::Det (int slice, Array<dVec,1> &tempPath)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);

  int N = last - first + 1;
  Array<double,2> detMatrix(N,N);
  // Fill up determinant matrix
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl-first, dist, diff, tempPath);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        action += ActionSplines(sliceDiff)[dim](diff[dim]);
      detMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  return Determinant (detMatrix);
}


double FreeNodalActionClass::Det (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);

  int N = last - first + 1;
  Array<double,2> DetMatrix2(N,N);
  // Fill up determinant matrix
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        action += ActionSplines(sliceDiff)[dim](diff[dim]);
      DetMatrix2(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  return Determinant (DetMatrix2);
}


bool FreeNodalActionClass::IsPositive (int slice)
{
  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = slice - refSlice;
  sliceDiff = min(sliceDiff, Path.TotalNumSlices-sliceDiff);
  //cerr << "slice=" << slice << " sliceDiff = " << sliceDiff << endl;
  if (sliceDiff == 0)
    return (true);
  else
    return (Det(slice) > 0.0);
}


void FreeNodalActionClass::GradientDet (int slice, double &det, Array<dVec,1> &gradient, Array<dVec,1> &tempPath)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);

  int N = last - first + 1;
  Array<double,2> detMatrix(N,N), cofactors(N,N);

  // Fill up determinant matrix
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl-first, dist, diff, tempPath);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        action += ActionSplines(sliceDiff)[dim](diff[dim]);
      detMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  // Compute determinant
  det = Determinant (detMatrix);
  if (det < 0.0 && !PathData.Path.UseNodeImportance) {
    return;
  }

  // Check if singular
  if (det == 0.0 || isnan(det)) {
    if (((nSingular)%100000) == 99999) {
      cerr << "Warning: Num Singular Matrices = " << nSingular << endl;
    }
    cerr << "Warning: Singular Matrix at slice: " << slice << endl;
    nSingular++;
    det = -1.0;
    gradient(0)[0] = sqrt(-1.0);
    return;
  }

  cofactors = detMatrix;
  GJInverse (cofactors);
  Transpose (cofactors);
  cofactors = det * cofactors;

  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl-first) = 0.0;
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl-first, dist, diff, tempPath);
      dVec gradPhi;
      for (int dim=0; dim<NDIM; dim++)
        gradPhi[dim] = -ActionSplines(sliceDiff)[dim].Deriv(diff[dim]) * detMatrix(refPtcl-first, ptcl-first);
      //dVec gradPhi = -2.0*C*diff*DetMatrix(refPtcl-first,ptcl-first);

      gradient(ptcl-first) = gradient(ptcl-first)+ gradPhi*cofactors(refPtcl-first, ptcl-first);
    }
  }
}



void FreeNodalActionClass::GradientDet (int slice, double &det, Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);

  int N = last - first + 1;
  Array<double,2> DetMatrix2(N,N), Cofactors2(N,N);

  // Fill up determinant matrix
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        action += ActionSplines(sliceDiff)[dim](diff[dim]);
      DetMatrix2(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  // Compute determinant
  det = Determinant (DetMatrix2);
  if (det < 0.0 && !PathData.Path.UseNodeImportance) {
    return;
  }

  // Check if singular
  if (det == 0.0 || isnan(det)) {
    if (((nSingular)%100000) == 99999) {
      cerr << "Num Singular Matrices = " << nSingular << endl;
    }
    // cerr << "Singular Matrix at slice: " << slice << endl;
    nSingular++;
    det = -1.0;
    gradient(0)[0] = sqrt(-1.0);
    return;
  }

  Cofactors2 = DetMatrix2;
  GJInverse (Cofactors2);
  Transpose (Cofactors2);
  Cofactors2 = det * Cofactors2;

  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl-first) = 0.0;
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      dVec gradPhi;
      for (int dim=0; dim<NDIM; dim++)
        gradPhi[dim] = -ActionSplines(sliceDiff)[dim].Deriv(diff[dim]) * DetMatrix2(refPtcl-first, ptcl-first);
      //dVec gradPhi = -2.0*C*diff*DetMatrix(refPtcl-first,ptcl-first);

      gradient(ptcl-first) = gradient(ptcl-first)+ gradPhi*Cofactors2(refPtcl-first, ptcl-first);
    }
  }
}



void 
FreeNodalActionClass::GradientDetFD (int slice, double &det, Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  double t = abs(refSlice-slice) * PathData.Path.tau;
  double beta = PathData.Path.TotalNumSlices * PathData.Path.tau;
  t = min (t, fabs(beta-t));
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);

  // HACK HACK HACK for now;  should work for serial mode.
//   if (Path.GetRefSlice() < Path.NumTimeSlices())
//     for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
//       Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
        action += ActionSplines(sliceDiff)[dim](diff[dim]);
      DetMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }


  // Compute determinant
  det = Determinant (DetMatrix);

  dVec disp;
  double dist;

  double eps = 1.0e-5;
  for (int ptcl=first; ptcl <= last; ptcl++) {
    dVec delta = 0.0;
    dVec &r = Path(slice, ptcl);
    double dplus, dminus;
    for (int dim=0; dim<NDIM; dim++) {
      delta = 0.0;
      delta[dim] = 0.5*eps;
      for (int ref=first; ref <= last; ref++) {
        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
        dVec diff = disp + delta;
        double action = 0.0;
        for (int dim=0; dim<NDIM; dim++)
          action += ActionSplines(sliceDiff)[dim](diff[dim]);
        DetMatrix(ref-first, ptcl-first) = exp(-action);
      }
      dplus = Determinant (DetMatrix);
      for (int ref=first; ref <= last; ref++) {
        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
        dVec diff = disp - delta;
        double action = 0.0;
        for (int dim=0; dim<NDIM; dim++)
          action += ActionSplines(sliceDiff)[dim](diff[dim]);
        DetMatrix(ref-first, ptcl-first) = exp(-action);
      }
      dminus = Determinant(DetMatrix);
      for (int ref=first; ref <= last; ref++) {
        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
        dVec diff = disp;
        double action = 0.0;
        for (int dim=0; dim<NDIM; dim++)
          action += ActionSplines(sliceDiff)[dim](diff[dim]);
        DetMatrix(ref-first, ptcl-first) = exp(-action);
      }
      gradient(ptcl-first)[dim] = (dplus-dminus)/eps;
    }
  }
}

/// This simply returns the 1st Newton-Raphson estimate of the nodal distance
double FreeNodalActionClass::NodalDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  double det;
  int N = last-first+1;
  GradientDet (slice, det, GradVec);
  double grad2 = 0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec(i), GradVec(i));
  double dist = det/sqrt(grad2);
  return (dist);
}


/// HybridDist first computes a distance with the 1/grad(ln(det))
/// method.  If this method says we're far from the nodes, we go with
/// that.  If it says we're close, call LineSearchDist to get a more
/// accurate value.
double FreeNodalActionClass::HybridDist (int slice, double lambdaTau)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  double det;
  int N = last-first+1;

  Array<dVec,1> GradVec2(N);
  GradientDet (slice, det, GradVec2);

  if (det < 0.0 && !PathData.Path.UseNodeImportance)
    return -1.0;

  double grad2 = 0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec2(i), GradVec2(i));

  double gradDist = det/sqrt(grad2);

  //if (((NumGradDists+NumLineDists)%1000000) == 999999) {
  //  cerr << "Percent line searches = " << (double)NumLineDists/(NumGradDists+NumLineDists) << endl;
  //}

  // gradDist will almost always be a lower bound to the real
  // distance.  Therefore, if says we are far from the nodes, we
  // probably are and we can just use its value.
  if (gradDist > sqrt(4.0*lambdaTau)) {
    NumGradDists++;
    return (gradDist);
  }
  // However, if gradDist says we are close, we should check to see if
  // it is correct with a more accurate bisection search.
  else {
    NumLineDists++;
    return LineSearchDist(slice);
  }
  //double lineSearchDist = NewtonRaphsonDist(slice);// LineSearchDist(slice);
  //cout << slice << " " << sqrt(4.0*lambdaTau) << " " << gradDist << " " << lineSearchDist << " " << NewtonRaphsonDist(slice) << endl;
  //return lineSearchDist;
}


double FreeNodalActionClass::MaxDist(int slice)
{
  if (Det(slice) < 0.0)
    return (-1.0);

  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  dVec disp;
  double dist;
  double minDist = 1.0e300;
  for (int ptcl1=first; ptcl1<last; ptcl1++) {
    for (int ptcl2=ptcl1+1; ptcl2<=last; ptcl2++) {
      Path.DistDisp (slice, ptcl1, ptcl2, dist, disp);
      minDist = (dist < minDist) ? dist : minDist;
    }
  }
  return minDist*M_SQRT1_2;
}




double FreeNodalActionClass::LineSearchDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  double epsilon = 1.0e-4 * sqrt (4.0*species.lambda*Path.tau);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  //  double maxDist = MaxDist (slice);
  double det0, det;
  int N = last-first+1;
  double retVal;

  Array<dVec,1> gradVec(N), tempPath(N), savePath(N);

  // Save the current path
  for (int i=0; i < N; i++) {
    savePath(i) = Path(slice,i+first);
    tempPath(i) = Path(slice,i+first);
  }

  GradientDet (slice, det0, gradVec, tempPath);
  if (det0 < 0.0 && !PathData.Path.UseNodeImportance)
    return -1.0;

  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (gradVec(i), gradVec(i));
  double gradMag = sqrt (grad2);

  for (int i=0; i<N; i++)
    gradVec(i) = (1.0/gradMag)*gradVec(i);

  double dist = abs(det0/gradMag);

  double minFactor, maxFactor, tryFactor, newDet;
  minFactor = 0.0;
  maxFactor = 0.5;

  bool done = false;
  // First, find first sign change
  det = det0;
  double maxDist=sqrt(dot(PathData.Path.GetBox(),PathData.Path.GetBox()));
  while ((det*det0 > 0.0) && ((maxFactor*dist)<maxDist)) {
    maxFactor *= 2.0;
    for (int i=0; i<N; i++)
      tempPath(i) = savePath(i) - maxFactor*dist*gradVec(i);
    det = Det(slice,tempPath);
  }

  if (det*det0 >= 0.0)
    retVal = maxDist;
  else {
    // Now, do a bisection search for the sign change.
    while (((maxFactor-minFactor)*dist > epsilon) && (minFactor*dist < maxDist)) {
      tryFactor = 0.5*(maxFactor+minFactor);
      for (int i=0; i<N; i++)
        tempPath(i) = savePath(i) - tryFactor*dist*gradVec(i);
      det = Det (slice,tempPath);
      if (det*det0 > 0.0)
        minFactor = tryFactor;
      else
        maxFactor = tryFactor;
    }
    if (minFactor*dist >= maxDist)
      retVal = maxDist;
    else
      retVal = dist * tryFactor;
  }

  // Restore original Path position
  //for (int i=0; i<N; i++)
  //  Path (slice, i+first) = SavePath2(i);
  return retVal; //min (maxDist, retVal);
}



double FreeNodalActionClass::NewtonRaphsonDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  const int maxIter = 15;

  double det;
  int N = last-first+1;
  double retVal = 0.0;
  double maxBox = sqrt (dot(Path.GetBox(), Path.GetBox()));

  if (Det(slice) < 0.0)
    return (-1.0);

  // Save the current path
  Array<dVec,1> savePath(N), tempPath(N), lastPath(N);
  for (int i=0; i<N; i++) {
    savePath(i) = Path(slice,i+first);
    tempPath(i) = Path(slice,i+first);
    lastPath(i) = Path(slice,i+first);
  }

  bool done = false;
  int numIter = 0;
  // Do Newton-Raphson iterations
  GradientDet (slice, det, GradVec, tempPath);
  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec(i), GradVec(i));
  double gradMag = sqrt (grad2);
  double firstDist = det/gradMag;
  double dist;

  while (!done && (numIter < maxIter)) {
    dist = 2.0*det/gradMag;
    do {
      dist *= 0.5;
      for (int i=0; i<N; i++)
        tempPath(i) = lastPath(i) -(dist/gradMag)*GradVec(i);
    } while (Det (slice, tempPath) < 0.0);

    for (int i=0; i<N; i++)
      lastPath(i) = tempPath(i);

    if (dist < 0.001*firstDist)
      done = true;
    GradientDet (slice, det, GradVec, tempPath);
    // If we come to a singular determinant matrix...
    if (isnan(GradVec(0)[0])) {
      //for (int i=0; i<N; i++)
      //  Path(slice,first+i) = savePath(i);
      //cout << "hi" << endl;
      return (LineSearchDist(slice));
    }
    grad2 = 0.0;
    for (int i=0; i<N; i++)
      grad2 += dot (GradVec(i), GradVec(i));
    gradMag = sqrt (grad2);
    numIter++;
  }

  double totalDist = 0.0;
  for (int i=0; i<N; i++) {
    dVec diff = tempPath(i) - savePath(i);
    totalDist += dot (diff,diff);
  }
  totalDist = sqrt (totalDist);

  // Restore original Path
  //for (int i=0; i<N; i++)
  //  Path(slice,i+first)  = savePath(i);

  return totalDist;
}



void FreeNodalActionClass::Read (IOSectionClass &in)
{
  TimeSpent = 0.0;
  SetupFreeActions();
}


double FreeNodalActionClass::SingleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  if (PathData.Path.Equilibrate)
    return SimpleAction(startSlice,endSlice,changePtcls,level);
  else if (PathData.Path.UseNodeImportance)
    return NodeImportanceAction(startSlice,endSlice,changePtcls,level);
  else
    return PreciseAction(startSlice,endSlice,changePtcls,level);
}


/// Return essentially 0 or infinity
double FreeNodalActionClass::SimpleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  int skip = 1<<level;
  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;
  //std::cout << "FreeNodalAction --------------" << endl;

  int totalSlices = Path.TotalNumSlices;
  int numSlices = (endSlice - startSlice)/skip + 1;
  double deter[numSlices];
  for (int i=0; i<numSlices; i++)
    deter[i] = 0.0;
  double uNode = 0.0;
  bool abort = 0;
  #pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    #pragma omp flush (abort)
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    if (!sliceIsRef&&!abort) {
      int i = (slice - startSlice)/skip;
      if (((GetMode()==NEWMODE)||FirstDetTime)||!PathData.Path.UseNodeDet)
        deter[i] = Det(slice);
      else
        deter[i] = PathData.Path.NodeDet(slice,SpeciesNum);
      if (deter[i] <= 0.0) {
        #pragma omp critical
        {
            abort = 1;
        }
      }
      //std::cout<<slice<<" "<<deter<<endl;
    }
  }
  #pragma omp barrier

  if(abort)
    uNode = 1.0e100;
  else {
    if (((level==0 && GetMode()==NEWMODE) || FirstDetTime) && PathData.Path.UseNodeDet) {
      for (int slice=startSlice; slice <= endSlice; slice+=skip) {
        int i = (slice - startSlice)/skip;
        bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
        if (!sliceIsRef)
          PathData.Path.NodeDet(slice,SpeciesNum) = deter[i];
      }
      FirstDetTime = 0;
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  //std::cout << "-----------------" << level << " " << uNode << " " << deter << " " << brokenSlice << endl;
  return uNode;
}


double FreeNodalActionClass::NodeImportanceAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double uNode = 0.0;

  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Path.tau * (double)skip;

  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;

  int totalSlices = Path.TotalNumSlices;
  int numSlices = (endSlice - startSlice)/skip + 1;
  if (numSlices < 2)
    cerr << "ERROR: numSlices < 2 in FreeNodalAction" << endl;
  double dist[numSlices];

  for (int i=0; i<numSlices; i++)
    dist[i] = 0.0;
  #pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    if (!sliceIsRef) {
      int i = (slice - startSlice)/skip;
      if (((GetMode()==NEWMODE)||FirstDistTime)||!PathData.Path.UseNodeDist)
        dist[i] = HybridDist (slice, lambda*levelTau);
      else
        dist[i] = PathData.Path.NodeDist(slice,SpeciesNum);
      //cout << PathData.Path.CloneStr << " " << SpeciesNum << " " << refSlice << " " << i << " " << dist[i] << " " << slice << endl;
    }
  }
  #pragma omp barrier

  for (int slice = startSlice; slice < endSlice; slice+=skip) {
    int i = (slice - startSlice)/skip;

    bool slice1IsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    bool slice2IsRef = (slice+skip == refSlice) || (slice+skip == refSlice+totalSlices);
    double dist1 = abs(dist[i]);
    double dist2 = abs(dist[i+1]);

    if (slice1IsRef || (dist1==0.0))
      uNode -= log1p(-exp(-dist2*dist2/(lambda*levelTau)));
    else if (slice2IsRef || (dist2==0.0))
      uNode -= log1p(-exp(-dist1*dist1/(lambda*levelTau)));
    else
      uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
    if (((level==0 && GetMode()==NEWMODE) || FirstDistTime) && PathData.Path.UseNodeDist) {
      PathData.Path.NodeDist(slice,SpeciesNum) = dist1;
      PathData.Path.NodeDist(slice+skip,SpeciesNum) = dist2;
      FirstDistTime = 0;
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);

  double eps = PathData.Path.NodeImpEps;
  if (eps != 0.0)
    uNode = log(eps) + log1p(exp(uNode)/eps);

  return uNode;
}


double FreeNodalActionClass::PreciseAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double uNode = 0.0;

  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Path.tau * (double)skip;

  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;

  int totalSlices = Path.TotalNumSlices;
  int numSlices = (endSlice - startSlice)/skip + 1;
  if (numSlices < 2)
    cerr << "ERROR: numSlices < 2 in FreeNodalAction" << endl;
  double dist[numSlices];

  for (int i=0; i<numSlices; i++)
    dist[i] = 0.0;
  bool abort = 0;
  //if (GetMode()==OLDMODE)
  //  cout << PathData.Path.CloneStr << " start, end " << startSlice << ", " << endSlice << "!!!!!!!!!!!!!!!!!!!ACTION(OLD)" << endl;
  //else
  //  cout << PathData.Path.CloneStr << " start, end " << startSlice << ", " << endSlice << "!!!!!!!!!!!!!!!!!!!ACTION(NEW)" << endl;
  #pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    #pragma omp flush (abort)
    if (!abort) {
      bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
      if (!sliceIsRef&&!abort) {
        int i = (slice - startSlice)/skip;
        if (((GetMode()==NEWMODE)||FirstDistTime)||!PathData.Path.UseNodeDist)
          dist[i] = HybridDist (slice, lambda*levelTau);
        else
          dist[i] = PathData.Path.NodeDist(slice,SpeciesNum);
        //cout << PathData.Path.CloneStr << " " << SpeciesNum << " " << refSlice << " " << i << " " << dist[i] << " " << slice << endl;
        if (dist[i] < 0.0) {
          #pragma omp critical
          {
            abort = 1;
          }
        }
      }
    }
  }
  #pragma omp barrier

  if (abort)
    uNode = 1.0e100;
  else {
    for (int slice = startSlice; slice < endSlice; slice+=skip) {
      int i = (slice - startSlice)/skip;

      bool slice1IsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
      bool slice2IsRef = (slice+skip == refSlice) || (slice+skip == refSlice+totalSlices);
      double dist1 = dist[i];
      double dist2 = dist[i+1];

      if (!slice1IsRef && (dist1<0.0))
        abort = 1;
      else if (!slice2IsRef && (dist2<0.0))
        abort = 1;
      else if (slice1IsRef || (dist1==0.0))
        uNode -= log1p(-exp(-dist2*dist2/(lambda*levelTau)));
      else if (slice2IsRef || (dist2==0.0))
        uNode -= log1p(-exp(-dist1*dist1/(lambda*levelTau)));
      else
        uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
      if (!abort && ((level==0 && GetMode()==NEWMODE) || FirstDistTime) && PathData.Path.UseNodeDist) {
        PathData.Path.NodeDist(slice,SpeciesNum) = dist1;
        PathData.Path.NodeDist(slice+skip,SpeciesNum) = dist2;
        FirstDistTime = 0;
      }
    }
  }

  if (abort)
    uNode = 1.0e100;

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return uNode;
}


double FreeNodalActionClass::d_dBeta (int slice1, int slice2, int level)
{
  Array<int,1> changedPtcls(1);

  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Path.tau * (double)skip;

  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;
  int sliceDiff1 = abs(slice1-refSlice);
  sliceDiff1 = min (sliceDiff1, PathData.Path.TotalNumSlices-sliceDiff1);

  int totalSlices = Path.TotalNumSlices;
  int numSlices = (slice2 - slice1)/skip + 1;
  if (numSlices < 2)
    cerr << "ERROR: numSlices < 2 in FreeNodalAction" << endl;
  double dist[numSlices];
  for (int i=0; i<numSlices; i++)
    dist[i] = 0.0;
  bool abort = 0;
  //if (GetMode()==OLDMODE)
  //  cout << PathData.Path.CloneStr << " start, end " << slice1 << ", " << slice2 << "!!!!!!!!!!!!!!!!!!!ENERGY(OLD)" << endl;
  //else
  //  cout << PathData.Path.CloneStr << " start, end " << slice1 << ", " << slice2 << "!!!!!!!!!!!!!!!!!!!ENERGY(NEW)" << endl;
  #pragma omp parallel for
  for (int slice=slice1; slice <= slice2; slice+=skip) {
    #pragma omp flush (abort)
    if (!abort) {
      bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
      if (!sliceIsRef&&!abort) {
        int i = (slice - slice1)/skip;
        if (FirstDistTime||!PathData.Path.UseNodeDist)
          if (PathData.Path.Equilibrate)
            dist[i] = NodalDist (slice);
          else
            dist[i] = HybridDist (slice,lambda*levelTau);
        else
          dist[i] = PathData.Path.NodeDist(slice,SpeciesNum);
        //cout << PathData.Path.CloneStr << " " << SpeciesNum << " " << refSlice << " " << i << " " << dist[i] << " " << slice << endl;
        if (dist[i] < 0.0) {
          #pragma omp critical
          {
            cerr << PathData.Path.CloneStr << " ERROR: dist = " << dist[i] << " skip = " << skip << " slice2 = " << slice+skip << " refSlice = " << refSlice << " species = " << species.Name << endl;
            dist[i] = 0.0;
          }
        }
      }
    }
  }
  //#pragma omp barrier

  double uNode = 0.0;
  int i = 0;
  for (int slice=slice1; slice < slice2; slice+=skip) {
    int sliceDiff = slice - refSlice;
    sliceDiff = min (sliceDiff, PathData.Path.TotalNumSlices - sliceDiff);
    bool slice1IsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    bool slice2IsRef = (slice+skip == refSlice) || (slice+skip == refSlice+totalSlices);
    double dist1 = dist[i];
    double dist2 = dist[i+1];
    if (FirstDistTime && PathData.Path.UseNodeDist) {
      PathData.Path.NodeDist(slice,SpeciesNum) = dist1;
      PathData.Path.NodeDist(slice+skip,SpeciesNum) = dist2;
      FirstDistTime = 0;
    }
    double prod;
    if (slice1IsRef || (dist1==0.0))
      prod = dist2*dist2;
    else if (slice2IsRef || (dist2==0.0))
      prod = dist1*dist1;
    else
      prod = dist1*dist2;

    double prod_llt = prod/(lambda*levelTau);
    double exp_m1 = expm1 (prod_llt);
    if (fpclassify(exp_m1)==FP_NORMAL)
      if (!isnan(exp_m1) && !isinf(exp_m1))
        if (isnormal(exp_m1))
          if ((!isnan(prod_llt)) && (exp_m1 != 0.0))
            uNode += prod_llt / (levelTau*exp_m1);
    if (isnan(uNode))
      cerr << "ERROR: uNode broken again!\n";
    i += 1;
  }

  if (abort)
    uNode = 1.0e100;

  //  return 0.0;
  return uNode/(double)Path.TotalNumSlices;
}


NodeType FreeNodalActionClass::Type()
{
  return FREE_PARTICLE;
}


bool FreeNodalActionClass::IsGroundState()
{
  return (false);
}

void FreeNodalActionClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "FREE_PARTICLE");
}

string FreeNodalActionClass::GetName()
{
  return "FreeNodal";
}
