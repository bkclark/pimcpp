/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the SHO Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////


#include "SHONodalActionClass.h"
#include "../PathDataClass.h"
#include "../MatrixOps/MatrixOps.h"

void SHONodalActionClass::SetupSHOActions()
{
  double tau = Path.tau;
  c1.resize(Path.TotalNumSlices);
  c2.resize(Path.TotalNumSlices);
  c3.resize(Path.TotalNumSlices);
  for (int i = 1; i < Path.TotalNumSlices; ++i) {
    c1(i) = pow((omega/(2.0*pi*sinh(i*tau*omega))),NDIM/2.0);
    c2(i) = omega/(2.0*sinh(i*tau*omega));
    c3(i) = cosh(i*tau*omega);
  }
}


SHONodalActionClass::SHONodalActionClass (PathDataClass &pathData,int speciesNum) :
  NodalActionClass (pathData),
  Path (pathData.Path),
  SpeciesNum (speciesNum)
{
  int N = Path.Species(speciesNum).LastPtcl - Path.Species(speciesNum).FirstPtcl+1;
  DetMatrix.resize(N,N);
  Cofactors.resize(N,N);
  GradVec.resize(N);
  SavePath.resize(N);
  NumLineDists = NumGradDists = 0;
}


double SHONodalActionClass::GetAction(int refSlice, int slice, int sliceDiff, int refPtcl, int ptcl)
{
  double r0r0 = dot(Path(refSlice,refPtcl),Path(refSlice,refPtcl));
  double r1r1 = dot(Path(slice,ptcl),Path(slice,ptcl));
  double r0r1 = dot(Path(refSlice,refPtcl),Path(slice,ptcl));
  //cerr << refSlice << " " << slice << " " << sliceDiff << " " << refPtcl << " " << ptcl << " " << r0r0 << " " << r1r1 << " " << r0r1 << endl;

  return c2(sliceDiff)*(c3(sliceDiff)*(r0r0+r1r1) - 2.*r0r1) + log(c1(sliceDiff));
}


double SHONodalActionClass::Det (int slice)
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

  // Fill up determinant matrix
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      double action = GetAction(refSlice,slice,sliceDiff,refPtcl,ptcl);
      DetMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  return Determinant (DetMatrix);
}


bool SHONodalActionClass::IsPositive (int slice)
{
  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = slice - refSlice;
  sliceDiff = min(sliceDiff, Path.TotalNumSlices-sliceDiff);
  cerr << "slice=" << slice << " sliceDiff = " << sliceDiff << endl;
  if (sliceDiff == 0)
    return (true);
  else
    return (Det(slice) > 0.0);
}


void SHONodalActionClass::GradientDet (int slice, double &det, Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  // Fill up determinant matrix
  int myStartSlice, myEndSlice;
  int myProc = PathData.Path.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);

  bool singular = false;
  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      double action = GetAction(refSlice,slice,sliceDiff,refPtcl,ptcl);
      DetMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }

  if (singular)
    cerr << "Singular at slice=" << slice << endl;

  // Compute determinant
  det = Determinant (DetMatrix);
  if (singular) {
    gradient(0)[0] = sqrt(-1.0);
    return;
  }
  if (fabs(det)<1.0e-40) {
    cerr << "DetMatrix = " << DetMatrix << endl;
    for (int ptcl=first; ptcl<=last; ptcl++)
      cerr << "ptcl " << ptcl << " = " << Path(slice, ptcl) << endl;
    cerr << "Permutation = " << Path.Permutation.data() << endl;
  }
  Cofactors = DetMatrix;
  GJInverse (Cofactors);
  Transpose (Cofactors);
  Cofactors = det * Cofactors;


  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl-first) = 0.0;
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec gradPhi;
      for (int dim=0; dim<NDIM; dim++)
        gradPhi[dim] = GetAction(refSlice,slice,sliceDiff,refPtcl,ptcl) * DetMatrix(refPtcl-first, ptcl-first);
      //dVec gradPhi = -2.0*C*diff*DetMatrix(refPtcl-first,ptcl-first);
      gradient(ptcl-first) = gradient(ptcl-first) + gradPhi*Cofactors(refPtcl-first, ptcl-first);
    }
  }
}


double SHONodalActionClass::NodalDist (int slice)
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
double SHONodalActionClass::HybridDist (int slice, double lambdaTau)
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

  if (det < 0.0)
    return -1.0;

  double gradDist = det/sqrt(grad2);

  if (((NumGradDists+NumLineDists)%1000000) == 999999) {
    cerr << "Percent line searches = " << (double)NumLineDists/(NumGradDists+NumLineDists) << endl;
  }

  // gradDist will almost always be a lower bound to the real
  // distance.  Therefore, if says we are far from the nodes, we
  // probably are and we can just use its value.
  if (gradDist > sqrt(5.0*lambdaTau)) {
    NumGradDists++;
    return (gradDist);
  }
  // However, if gradDist says we are close, we should check to see if
  // it is correct with a more accurate bisection search.
  else {
    NumLineDists++;
    return LineSearchDist(slice);
  }
}


double SHONodalActionClass::MaxDist(int slice)
{
  if (Det(slice) < 0.0)
    return (-1.0);

  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  dVec disp;
  double dist;
  double minDist = 1.0e300;
  for (int ptcl1=first; ptcl1<last; ptcl1++)
    for (int ptcl2=ptcl1+1; ptcl2<=last; ptcl2++) {
      Path.DistDisp (slice, ptcl1, ptcl2, dist, disp);
      minDist = (dist < minDist) ? dist : minDist;
    }
  return minDist*M_SQRT1_2;
}


double SHONodalActionClass::LineSearchDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  double epsilon = 1.0e-4 * sqrt (4.0*species.lambda*Path.tau);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  // double maxDist = MaxDist (slice);
  double det0, det;
  int N = last-first+1;
  double retVal;

  // Save the current path
  for (int i=0; i < N; i++)
    SavePath(i) = Path(slice,i+first);

  GradientDet (slice, det0, GradVec);
  if (det0 < 0.0)
    return -1.0;

  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec(i), GradVec(i));
  double gradMag = sqrt (grad2);

  for (int i=0; i<N; i++)
    GradVec(i) = (1.0/gradMag)*GradVec(i);

  double dist = det0/gradMag;


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
      Path (slice, i+first) = SavePath(i) - maxFactor*dist*GradVec(i);
    det = Det(slice);
  }

  if (det*det0 >= 0.0)
    retVal = maxDist;
  else {
    // Now, do a bisection search for the sign change.
    while (((maxFactor-minFactor)*dist > epsilon) && (minFactor*dist < maxDist)) {
      tryFactor = 0.5*(maxFactor+minFactor);
      for (int i=0; i<N; i++)
        Path (slice, i+first) = SavePath(i) - tryFactor*dist*GradVec(i);
      det = Det (slice);
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
  for (int i=0; i<N; i++)
    Path (slice, i+first) = SavePath(i);
  return retVal; //min (maxDist, retVal);
}


double SHONodalActionClass::NewtonRaphsonDist (int slice)
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
  Array<dVec,1> savePath(N), lastPath(N);
  for (int i=0; i<N; i++) {
    savePath(i) = Path(slice,i+first);
    lastPath(i) = Path(slice,i+first);
  }

  bool done = false;
  int numIter = 0;
  // Do Newton-Raphson iterations
  GradientDet (slice, det, GradVec);
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
	Path(slice, first+i) = lastPath(i) -(dist/gradMag)*GradVec(i);
    } while (Det (slice) < 0.0);

    for (int i=0; i<N; i++)
      lastPath(i) = Path(slice, first+i);
    
    if (dist < 0.001*firstDist)
      done = true;
    GradientDet (slice, det, GradVec);
    // If we come to a singular determinant matrix...
    if (isnan(GradVec(0)[0])) {
      for (int i=0; i<N; i++)
	Path(slice,first+i) = savePath(i);
      return (LineSearchDist(slice));
    }
    grad2 = 0.0;
    for (int i=0; i<N; i++)
      grad2 += dot (GradVec(i), GradVec(i));
    gradMag = sqrt (grad2);
    numIter++;
  }
    

//     if (det < 0.0) {
//       for (int i=0; i<N; i++)
// 	Path(slice, first+i) -= (dist/gradMag)*GradVec(i);
//       done = true;
//     }
//     else {
//       for (int i=0; i<N; i++)
// 	lastPath(i) = Path(slice, first+i);
//       double grad2=0.0;
//       for (int i=0; i<N; i++)
// 	grad2 += dot (GradVec(i), GradVec(i));
//       double gradMag = sqrt (grad2);
//       double dist = det/gradMag;
//       for (int i=0; i<N; i++)
// 	Path(slice, first+i) -= (dist/gradMag)*GradVec(i);
//       if (numIter > maxIter)
// 	done = true;
      
//       if ((det > 0.0) && (det < 0.0001*det0))
// 	done = true;
//       numIter++;
//     }
//   }
//   if (retVal != -1.0) {
  double totalDist = 0.0;
  for (int i=0; i<N; i++) {
    dVec diff = Path(slice, first+i) - savePath(i);
    totalDist += dot (diff,diff);
  }
  totalDist = sqrt (totalDist);

  // Restore original Path
  for (int i=0; i<N; i++)
    Path(slice,i+first)  = savePath(i);

  return totalDist;
}


void SHONodalActionClass::Read (IOSectionClass &in)
{
  if (!in.ReadVar ("omega",omega))
    omega = 1.0;

  // Decide which nodal distance function to use
  if (!in.ReadVar ("UseNoDist",UseNoDist))
    UseNoDist = 0;
  if (!in.ReadVar ("UseHybridDist",UseHybridDist))
    UseHybridDist = 0;
  if (!in.ReadVar ("UseNewtonRaphsonDist",UseNewtonRaphsonDist))
    UseNewtonRaphsonDist = 0;
  if (!in.ReadVar ("UseLineSearchDist",UseLineSearchDist))
    UseLineSearchDist = 0;
  if (!in.ReadVar ("UseMaxDist",UseMaxDist))
    UseMaxDist = 0;

  cout << "Nodal Distance Functions: " << UseHybridDist << " " << UseNewtonRaphsonDist << " " << UseLineSearchDist << " " << UseMaxDist << endl;

  SetupSHOActions();
}


void SHONodalActionClass::Init()
{
  FirstDistTime = 1;
  FirstDetTime = 1;
  // Initialize NodeDist and/or NodeDet
  if (PathData.Path.StoreNodeDist||PathData.Path.StoreNodeDet) {
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
        if (PathData.Path.StoreNodeDist)
          PathData.Path.NodeDist(slice,SpeciesNum) = HybridDist(slice, lambda*levelTau);
        if (PathData.Path.StoreNodeDet)
          PathData.Path.NodeDet(slice,SpeciesNum) = Det(slice);
      }
    }
    if (PathData.Path.StoreNodeDist)
      PathData.Path.NodeDist[OLDMODE](Range(startSlice,endSlice), SpeciesNum) =
        PathData.Path.NodeDist[NEWMODE](Range(startSlice,endSlice), SpeciesNum);
    if (PathData.Path.StoreNodeDet)
      PathData.Path.NodeDet[OLDMODE](Range(startSlice,endSlice), SpeciesNum) =
        PathData.Path.NodeDet[NEWMODE](Range(startSlice,endSlice), SpeciesNum);
  }
}


/// Return essentially 0 or infinity
double SHONodalActionClass::SimpleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  int skip = 1<<level;
  int myStart, myEnd;
  Path.SliceRange(PathData.Path.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;

  int totalSlices = Path.TotalNumSlices;
  int numSlices = (endSlice - startSlice)/skip + 1;
  double deter[numSlices];
  for (int i=0; i<numSlices; i++)
    deter[i] = 0.0;
  double uNode = 0.0;
  bool abort = 0;
  //#pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    //#pragma omp flush (abort)
    if (!sliceIsRef&&!abort) {
      int i = (slice - startSlice)/skip;
      if (((GetMode()==NEWMODE)||FirstDetTime)||!PathData.Path.StoreNodeDet)
        deter[i] = Det(slice);
      else
        deter[i] = PathData.Path.NodeDet(slice,SpeciesNum);
      //cerr << deter[i] << endl;
      if (deter[i] <= 0.0) {
        //#pragma omp critical
        {
            abort = 1;
        }
      }
    }
  }
  //#pragma omp barrier

  if(abort)
    uNode = 1.0e100;
  else {
    if (((level==0 && GetMode()==NEWMODE) || FirstDetTime) && PathData.Path.StoreNodeDet) {
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


double SHONodalActionClass::SingleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  if (PathData.Path.Equilibrate||UseNoDist)
    return SimpleAction(startSlice,endSlice,changePtcls,level);
  else
    return PreciseAction(startSlice,endSlice,changePtcls,level);
}


double SHONodalActionClass::PreciseAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
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
    cerr << "ERROR: numSlices < 2 in SHONodalAction" << endl;
  double dist[numSlices];

  for (int i=0; i<numSlices; i++)
    dist[i] = 0.0;
  bool abort = 0;
  //#pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    //#pragma omp flush (abort)
    if (!sliceIsRef&&!abort) {
      int i = (slice - startSlice)/skip;
      dist[i] = GetNodeDist(slice,lambda,levelTau,SpeciesNum);
      //cout << PathData.Path.CloneStr << " " << SpeciesNum << " " << refSlice << " " << i << " " << dist[i] << " " << slice << endl;
      if (dist[i] < 0.0) {
        //#pragma omp critical
        //{
          abort = 1;
        //}
      }
    }
  }
  //#pragma omp barrier

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
      if (!abort && ((level==0 && GetMode()==NEWMODE) || FirstDistTime) && PathData.Path.StoreNodeDist) {
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


double SHONodalActionClass::GetNodeDist(int slice, double lambda, double levelTau, int SpeciesNum)
{
 // cout << "HD: " << HybridDist(slice,lambda*levelTau) << " NR: " << NewtonRaphsonDist(slice) << " MD: " << MaxDist(slice) << " ND: " << NodalDist(slice) << endl;
  if (PathData.Path.NumParticles() == 1)
    return 1e100;

  double dist;
  if ((GetMode()==NEWMODE||FirstDistTime)||!PathData.Path.StoreNodeDist)
    if (UseHybridDist)
      dist = HybridDist(slice,lambda*levelTau); // Single Newton-Raphson, then Line Search
    else if (UseNewtonRaphsonDist)
      dist = NewtonRaphsonDist(slice); // Iterative Newton-Raphson
    else if (UseLineSearchDist)
      dist = LineSearchDist(slice); // Bisective Line Search
    else if (UseMaxDist)
      dist = MaxDist(slice); // Maximum Distance
    else
      dist = NodalDist(slice); // Single Newton-Raphson
  else
    dist = PathData.Path.NodeDist(slice,SpeciesNum);
  return dist;
}


double SHONodalActionClass::d_dBeta (int slice1, int slice2, int level)
{
  if (PathData.Path.Equilibrate||UseNoDist)
    return 0.0;

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
  //#pragma omp parallel for
  for (int slice=slice1; slice <= slice2; slice+=skip) {
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    if (!sliceIsRef) {
      int i = (slice - slice1)/skip;
      dist[i] = GetNodeDist(slice,lambda,levelTau,SpeciesNum);
      //cout << PathData.Path.CloneStr << " " << SpeciesNum << " " << refSlice << " " << i << " " << dist[i] << " " << slice << endl;
      if (dist[i] < 0.0) {
        cerr << PathData.Path.CloneStr << " ERROR: dist = " << dist[i] << " skip = " << skip << " slice2 = " << slice+skip << " refSlice = " << refSlice << " species = " << species.Name << endl;
        dist[i] = 0.0;
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
    if (FirstDistTime && PathData.Path.StoreNodeDist) {
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

  return uNode/(double)Path.TotalNumSlices;
}


NodeType SHONodalActionClass::Type()
{
  return FREE_PARTICLE;
}


bool SHONodalActionClass::IsGroundState()
{
  return (false);
}


void SHONodalActionClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "SHO");
}


string SHONodalActionClass::GetName()
{
  return "SHONodal";
}
