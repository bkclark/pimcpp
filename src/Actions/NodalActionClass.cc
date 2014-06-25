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

#include "NodalActionClass.h"
#include "../PathDataClass.h"
#include "../MatrixOps/MatrixOps.h"
#include "ctime"
#include "sys/time.h"

void NodalActionClass::AcceptCopy (int slice1, int slice2){}
void NodalActionClass::RejectCopy (int slice1, int slice2){}
void NodalActionClass::ChangeModel(int tmpModel){}
int NodalActionClass::GetModel(){}
int NodalActionClass::GetNumModels(){}
void NodalActionClass::SetupActions(){}
void NodalActionClass::Setk (Vec3 k){}
void NodalActionClass::Update(){}
double NodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl){}
double NodalActionClass::GetRhoij(int slice, int sliceDiff, int refPtcl, int ptcl, Array<dVec,1> &tempPath){}
void NodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix){}
void NodalActionClass::GetActionDeriv(int slice, int sliceDiff, int refPtcl, int ptcl, dVec &gradPhi, Array<double,2> &detMatrix, Array<dVec,1> &tempPath){}


void NodalActionClass::Init()
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


void NodalActionClass::Read (IOSectionClass &in)
{
  TimeSpent = 0.0;

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

  SetupActions();
}


/// This simply returns the 1st Newton-Raphson estimate of the nodal distance
double NodalActionClass::NodalDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  double det;
  int N = last-first+1;
  Array<dVec,1> gradVec(N);
  GradientDet (slice, det, gradVec);

  if (det < 0 && !PathData.Path.UseNodeImportance)
    return -1;

  double grad2 = 0.0;
  for (int i=0; i<N; i++)
    grad2 += dot(gradVec(i), gradVec(i));
  double dist = abs(det)/sqrt(grad2);
  return (dist);
}


bool NodalActionClass::IsPositive (int slice)
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


/// HybridDist first computes a distance with the 1/grad(ln(det))
/// method.  If this method says we're far from the nodes, we go with
/// that.  If it says we're close, call LineSearchDist to get a more
/// accurate value.
double NodalActionClass::HybridDist (int slice, double lambdaTau)
{
  //if (((NumGradDists+NumLineDists)%1000000) == 999999) {
  //  cerr << "Percent line searches = " << (double)NumLineDists/(NumGradDists+NumLineDists) << endl;
  //}

  double gradDist = NodalDist(slice);

  if (gradDist < 0 && !PathData.Path.UseNodeImportance)
    return -1;

  // gradDist will almost always be a lower bound to the real
  // distance.  Therefore, if says we are far from the nodes, we
  // probably are and we can just use its value.
  if (gradDist > sqrt(4.0*lambdaTau)) {
    // MaxDist is the distance to the nearest particle
    double maxDist = MaxDist(slice);
    if (gradDist < maxDist) {
      NumGradDists++;
      return (gradDist);
    } else
      return (maxDist);
  }
  // However, if gradDist says we are close, we should check to see if
  // it is correct with a more accurate bisection search.
  else {
    NumLineDists++;
    return LineSearchDist(slice);
  }
}


/// Returns the shortest distance between particles
double NodalActionClass::MaxDist(int slice)
{
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


/// Does a bisection line-search in the direction of the gradient
double NodalActionClass::LineSearchDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

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
    return -1;

  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (gradVec(i), gradVec(i));
  double gradMag = sqrt (grad2);

  for (int i=0; i<N; i++)
    gradVec(i) = (1.0/gradMag)*gradVec(i);

  double dist = det0/gradMag;

  double minFactor, maxFactor, tryFactor, newDet;
  minFactor = 0.0;
  maxFactor = 0.5;

  bool done = false;
  // First, find first sign change
  det = det0;
  double maxDist = MaxDist(slice);
  while ((det*det0 > 0.0) && ((maxFactor*dist)<maxDist)) {
    maxFactor *= 2.0;
    for (int i=0; i<N; i++)
      tempPath(i) = savePath(i) - maxFactor*dist*gradVec(i);
    det = Det(slice,tempPath);
  }

  double epsilon = 1.0e-4 * sqrt(4.0*species.lambda*Path.tau);
  if (det*det0 >= 0.0)
    retVal = maxDist;
  else {
    // Now, do a bisection search for the sign change.
    while (((maxFactor-minFactor)*abs(dist) > epsilon) && (minFactor*abs(dist) < maxDist)) {
      tryFactor = 0.5*(maxFactor+minFactor);
      for (int i=0; i<N; i++)
        tempPath(i) = savePath(i) - tryFactor*dist*gradVec(i);
      det = Det (slice,tempPath);
      if (det*det0 > 0.0)
        minFactor = tryFactor;
      else
        maxFactor = tryFactor;
    }
    if (minFactor*abs(dist) >= maxDist)
      retVal = maxDist;
    else
      retVal = abs(dist) * tryFactor;
  }

  return retVal;
}


double NodalActionClass::NewtonRaphsonDist (int slice)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  const int maxIter = 5;

  double det;
  int N = last-first+1;
  double retVal = 0.0;
  double maxBox = sqrt (dot(Path.GetBox(), Path.GetBox()));

  // Save the current path
  Array<dVec,1> gradVec(N), savePath(N), tempPath(N), lastPath(N);
  for (int i=0; i<N; i++) {
    savePath(i) = Path(slice,i+first);
    tempPath(i) = Path(slice,i+first);
    lastPath(i) = Path(slice,i+first);
  }

  bool done = false;
  int numIter = 0;
  // Do Newton-Raphson iterations
  GradientDet (slice, det, gradVec, tempPath);
  if (det<0.0 && !PathData.Path.UseNodeImportance)
    return (-1.0);

  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (gradVec(i), gradVec(i));
  double gradMag = sqrt(grad2);
  double firstDist = det/gradMag;
  double dist;

  while (!done && (numIter < maxIter)) {
    dist = 2.0*det/gradMag;
    do {
      dist *= 0.5;
      for (int i=0; i<N; i++)
        tempPath(i) = lastPath(i) - (dist/gradMag)*gradVec(i);
    } while (Det(slice, tempPath) < 0.0);

    for (int i=0; i<N; i++)
      lastPath(i) = tempPath(i);

    if (dist < 0.001*firstDist) {
      done = true;
    }
    GradientDet (slice, det, gradVec, tempPath);
    // If we come to a singular determinant matrix...
    if (isnan(gradVec(0)[0])) {
      cout << "Warning, singular matix" << endl;
      return (LineSearchDist(slice));
    }
    grad2 = 0.0;
    for (int i=0; i<N; i++)
      grad2 += dot (gradVec(i), gradVec(i));
    gradMag = sqrt (grad2);
    numIter++;
  }

  double totalDist = 0.0;
  for (int i=0; i<N; i++) {
    dVec diff = tempPath(i) - savePath(i);
    totalDist += dot (diff,diff);
  }
  totalDist = sqrt (totalDist);

  return totalDist;
}


double NodalActionClass::GetNodeDist(int slice, double lambda, double levelTau, int SpeciesNum)
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


double NodalActionClass::SingleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
{
  if (PathData.Path.Equilibrate||UseNoDist)
    return SimpleAction(startSlice,endSlice,changePtcls,level);
  else if (PathData.Path.UseNodeImportance)
    return NodeImportanceAction(startSlice,endSlice,changePtcls,level);
  else
    return PreciseAction(startSlice,endSlice,changePtcls,level);
}


/// Return essentially 0 or infinity
double NodalActionClass::SimpleAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
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
      if (deter[i] <= 0.0) {
        //#pragma omp critical
        //{
            abort = 1;
        //}
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


double NodalActionClass::NodeImportanceAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
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
    cerr << "ERROR: numSlices < 2 in NodalAction" << endl;
  double dist[numSlices];

  for (int i=0; i<numSlices; i++)
    dist[i] = 0.0;
  //#pragma omp parallel for
  for (int slice=startSlice; slice <= endSlice; slice+=skip) {
    bool sliceIsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    if (!sliceIsRef) {
      int i = (slice - startSlice)/skip;
      dist[i] = GetNodeDist(slice,lambda,levelTau,SpeciesNum);
      /// Method 4, minimum distance
      if (dist[i] < PathData.Path.NodeImpEps)
        dist[i] = PathData.Path.NodeImpEps;
    }
  }
  //#pragma omp barrier

  for (int slice = startSlice; slice < endSlice; slice+=skip) {
    int i = (slice - startSlice)/skip;

    bool slice1IsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
    bool slice2IsRef = (slice+skip == refSlice) || (slice+skip == refSlice+totalSlices);
    double dist1 = dist[i];
    double dist2 = dist[i+1];

    if (slice1IsRef || (dist1==0.0))
      uNode -= log1p(-exp(-dist2*dist2/(lambda*levelTau)));
    else if (slice2IsRef || (dist2==0.0))
      uNode -= log1p(-exp(-dist1*dist1/(lambda*levelTau)));
    else
      uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
    if (((level==0 && GetMode()==NEWMODE) || FirstDistTime) && PathData.Path.StoreNodeDist) {
      PathData.Path.NodeDist(slice,SpeciesNum) = dist1;
      PathData.Path.NodeDist(slice+skip,SpeciesNum) = dist2;
      FirstDistTime = 0;
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);

  ///// Node Importance part
  //double eps = PathData.Path.NodeImpEps;
  ////// Method 1, constant shift
  //if (eps != 0.0)
  //  uNode = -log(eps) - log1p(exp(-uNode)/eps);
  ///// Method 2, max value
  //if (exp(-uNode) < eps)
  //  uNode = -log(eps);
  ///// Method 3, calculate uNode from mean interparticle spacing
  //double uNodeTmp = 0.0;
  //for (int slice = startSlice; slice < endSlice; slice+=skip) {
  //  double dist1 = PathData.Path.NodeImpEps;
  //  double dist2 = PathData.Path.NodeImpEps;
  //  uNodeTmp -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
  //}
  //if (uNode > uNodeTmp)
  //  return uNodeTmp;

  return uNode;
}


double NodalActionClass::PreciseAction (int startSlice, int endSlice, const Array<int,1> &changePtcls, int level)
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
    cerr << "ERROR: numSlices < 2 in NodalAction" << endl;
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


double NodalActionClass::d_dBeta (int slice1, int slice2, int level)
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
    cerr << "ERROR: numSlices < 2 in NodalAction" << endl;
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


double NodalActionClass::Det (int slice, Array<dVec,1> &tempPath)
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
  double det = 0.;
  double scale = 1.;
  //do {

    // Fill up determinant matrix
    for (int refPtcl=first; refPtcl<=last; refPtcl++) {
      for (int ptcl=first; ptcl<=last; ptcl++) {
        double rhoij = GetRhoij(slice,sliceDiff,refPtcl,ptcl,tempPath);
        detMatrix(refPtcl-first, ptcl-first) = scale*rhoij;
      }
    }

    // Take determinant
    det = Determinant (detMatrix);

    // Decide scaling
    if (fabs(det) < 1.e-10)
      scale *= 2;
    else if (fabs(det) > 1.e10)
      scale = pow(scale,-1./N);
    //if (scale > 1e100 || scale < 1e-100)
    //  return 0;// = 2.e-10;

  //} while (fabs(det) < 1.e-10 || fabs(det) > 1.e10);

  return det;
}


double NodalActionClass::Det (int slice)
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
  double det = 0.;
  double scale = 1.;
  //do {

    // Fill up determinant matrix
    for (int refPtcl=first; refPtcl<=last; refPtcl++) {
      for (int ptcl=first; ptcl<=last; ptcl++) {
        double rhoij = GetRhoij(slice,sliceDiff,refPtcl,ptcl);
      //  cout << refPtcl << " " << ptcl << " " << rhoij << endl;
        detMatrix(refPtcl-first, ptcl-first) = scale*rhoij;
      }
    }

    // Take determinant
    det = Determinant (detMatrix);
    //cout << detMatrix << endl;

  //cout << scale << " " << slice << " " << det << endl;
    // Decide scaling
    if (fabs(det) < 1.e-50)
      scale *= 2;
    else if (fabs(det) > 1.e50)
      scale = pow(scale,1./N);

    //if (scale > 1e100 || scale < 1e-100) {
    //  det = 2.e-10;
    //}

  //} while (fabs(det) < 1.e-50 || fabs(det) > 1.e50);
  return det;
}


void NodalActionClass::GradientDet (int slice, double &det, Array<dVec,1> &gradient, Array<dVec,1> &tempPath)
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
  for (int refPtcl=first; refPtcl<=last; refPtcl++) {
    for (int ptcl=first; ptcl<=last; ptcl++) {
      double rhoij = GetRhoij(slice,sliceDiff,refPtcl,ptcl,tempPath);
      detMatrix(refPtcl-first, ptcl-first) = rhoij;
    }
  }

  // Compute determinant
  cofactors = detMatrix;
  det = GJInverse(cofactors);
  // Check if negative
  //if (det < 0.0 && !PathData.Path.UseNodeImportance) {
  //  return;
  //}
  // Check if singular
  if (det == 0.0 || isnan(det)) {
    if (((nSingular)%100000) == 99999) {
      cerr << "Warning: Num Singular Matrices = " << nSingular << endl;
    }
    //cerr << "Warning: Singular Matrix at slice: " << slice << endl;
    nSingular++;
    det = -1.0;
    gradient(0)[0] = sqrt(-1.0);
    return;
  }
  Transpose (cofactors);
  cofactors = det * cofactors;

  // Now compute gradient of determinant
  for (int ptcl=first; ptcl<=last; ptcl++) {
    gradient(ptcl-first) = 0.0;
    dVec gradPhi;
    for (int refPtcl=first; refPtcl<=last; refPtcl++) {
      GetActionDeriv(slice,sliceDiff,refPtcl,ptcl,gradPhi,detMatrix,tempPath);
      gradient(ptcl-first) = gradient(ptcl-first)+ gradPhi*cofactors(refPtcl-first, ptcl-first);
    }
  }
}


void NodalActionClass::GradientDet (int slice, double &det, Array<dVec,1> &gradient)
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
  for (int refPtcl=first; refPtcl<=last; refPtcl++) {
    for (int ptcl=first; ptcl<=last; ptcl++) {
      double rhoij = GetRhoij(slice,sliceDiff,refPtcl,ptcl);
      detMatrix(refPtcl-first, ptcl-first) = rhoij;
    }
  }

  // Compute determinant
  cofactors = detMatrix;
  det = GJInverse(cofactors);
  // Check if negative
  if (det < 0.0 && !PathData.Path.UseNodeImportance) {
    return;
  }
  // Check if singular
  if (det == 0.0 || isnan(det)) {
    if (((nSingular)%100000) == 99999) {
      cerr << "Warning: Num Singular Matrices = " << nSingular << endl;
    }
    //cerr << "Warning: Singular Matrix at slice: " << slice << endl;
    nSingular++;
    det = -1.0;
    gradient(0)[0] = sqrt(-1.0);
    return;
  }
  Transpose (cofactors);
  cofactors = det * cofactors;

  // Now compute gradient of determinant
  for (int ptcl=first; ptcl<=last; ptcl++) {
    gradient(ptcl-first) = 0.0;
    dVec gradPhi;
    for (int refPtcl=first; refPtcl<=last; refPtcl++) {
      GetActionDeriv(slice,sliceDiff,refPtcl,ptcl,gradPhi,detMatrix);
      gradient(ptcl-first) = gradient(ptcl-first) + gradPhi*cofactors(refPtcl-first, ptcl-first);
    }
  }
}


void NodalActionClass::GradientDetFD (int slice, double &det, Array<dVec,1> &gradient){}
//void NodalActionClass::GradientDetFD (int slice, double &det, Array<dVec,1> &gradient)
//{
//  SpeciesClass &species = Path.Species(SpeciesNum);
//  int first = species.FirstPtcl;
//  int last = species.LastPtcl;
//  // Fill up determinant matrix
//  int myStartSlice, myEndSlice;
//  int myProc = PathData.Path.Communicator.MyProc();
//  Path.SliceRange (myProc, myStartSlice, myEndSlice);
//  int refSlice = Path.GetRefSlice()-myStartSlice;
//  double t = abs(refSlice-slice) * PathData.Path.tau;
//  double beta = PathData.Path.TotalNumSlices * PathData.Path.tau;
//  t = min (t, fabs(beta-t));
//  double lambda = species.lambda;
//  double C = 1.0/(4.0*M_PI * lambda * t);
//
//  int sliceDiff = abs(slice-refSlice);
//  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
//  assert (sliceDiff <= Path.TotalNumSlices);
//
//  // HACK HACK HACK for now;  should work for serial mode.
////   if (Path.GetRefSlice() < Path.NumTimeSlices())
////     for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
////       Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);
//
//  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
//    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
//      dVec diff;
//      double dist;
//      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
//      double action = 0.0;
//      for (int dim=0; dim<NDIM; dim++)
//        action += ActionSplines(sliceDiff)[dim](diff[dim]);
//      DetMatrix(refPtcl-first, ptcl-first) = exp(-action);
//    }
//  }
//
//  // Compute determinant
//  det = Determinant (DetMatrix);
//
//  dVec disp;
//  double dist;
//
//  double eps = 1.0e-5;
//  for (int ptcl=first; ptcl <= last; ptcl++) {
//    dVec delta = 0.0;
//    dVec &r = Path(slice, ptcl);
//    double dplus, dminus;
//    for (int dim=0; dim<NDIM; dim++) {
//      delta = 0.0;
//      delta[dim] = 0.5*eps;
//      for (int ref=first; ref <= last; ref++) {
//        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
//        dVec diff = disp + delta;
//        double action = 0.0;
//        for (int dim=0; dim<NDIM; dim++)
//          action += ActionSplines(sliceDiff)[dim](diff[dim]);
//        DetMatrix(ref-first, ptcl-first) = exp(-action);
//      }
//      dplus = Determinant (DetMatrix);
//      for (int ref=first; ref <= last; ref++) {
//        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
//        dVec diff = disp - delta;
//        double action = 0.0;
//        for (int dim=0; dim<NDIM; dim++)
//          action += ActionSplines(sliceDiff)[dim](diff[dim]);
//        DetMatrix(ref-first, ptcl-first) = exp(-action);
//      }
//      dminus = Determinant(DetMatrix);
//      for (int ref=first; ref <= last; ref++) {
//        Path.RefDistDisp (slice, ref, ptcl, dist, disp);
//        dVec diff = disp;
//        double action = 0.0;
//        for (int dim=0; dim<NDIM; dim++)
//          action += ActionSplines(sliceDiff)[dim](diff[dim]);
//        DetMatrix(ref-first, ptcl-first) = exp(-action);
//      }
//      gradient(ptcl-first)[dim] = (dplus-dminus)/eps;
//    }
//  }
//}


bool NodalActionClass::IsGroundState()
{
  return (false);
}

void NodalActionClass::WriteInfo (IOSectionClass &out)
{
  // do nothing in base class.
}
