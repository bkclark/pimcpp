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

#include "PathClass.h"
#include "Actions/ActionsClass.h"

double PathClass::MinImageDistance(dVec v1, dVec v2)
{
  dVec disp = MinImageDisp(v1, v2);
  double dist = sqrt(dot(disp,disp));
  return dist;
}

dVec PathClass::MinImageDisp(dVec v1, dVec v2)
{
  dVec disp = v2-v1;
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  return disp;
}

const dVec& PathClass::GetOpenTail()
{
  return Path(OpenLink,NumParticles());
}

const dVec PathClass::ReturnOpenHead()
{
  return Path(OpenLink,OpenPtcl);
}

void PathClass::SetHead(const dVec &r,int join)
{
  (*this)(OpenLink, OpenPtcl) = r;
  if (OpenLink==NumTimeSlices()-1){
    if (join==NumTimeSlices()-1)
      (*this)(0,Permutation(OpenPtcl))=r;
    else 
      (*this)(0,OpenPtcl)=r;
  }
}

void PathClass::SetTail(const dVec &r)
{
  (*this)(OpenLink, NumParticles()) = r;
}

void PathClass::BackflowRefDistDisp (int slice, int refPtcl, int ptcl, double &dist, dVec &disp, Array<dVec,1> &tempPath)
{
  dVec pathPos = Path(slice,ptcl);
  dVec refPos = RefPath(refPtcl);
  for (int i=0; i<NumParticles(); i++) {
    // add in backflow coordinates
  }
  disp = pathPos - refPos;

  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));
}

void PathClass::RefDistDisp (int slice, int refPtcl, int ptcl, double &dist, dVec &disp)
{
  disp = Path(slice, ptcl)- RefPath(refPtcl);

  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl) -RefPath(refPtcl);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i)) 
      DBdisp(i) += Box(i);
    if (fabs(DBdisp(i)-disp(i)) > 1.0e-12){ 
      perr<<DBdisp(i)<<" "<<disp(i)<<endl;
    }
  }
#endif
}

void PathClass::RefDistDisp (int slice, int refPtcl, int ptcl, double &dist, dVec &disp, Array<dVec,1> &tempPath)
{
  disp = tempPath(ptcl) - RefPath(refPtcl);

  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl) -RefPath(refPtcl);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i)) 
      DBdisp(i) += Box(i);
    if (fabs(DBdisp(i)-disp(i)) > 1.0e-12){ 
      perr<<DBdisp(i)<<" "<<disp(i)<<endl;
    }
  }
#endif
}

void PathClass::SetIonConfig(int config)
{
  ConfigNum = config;
  int first = Species(IonSpecies).FirstPtcl;

  SetMode(OLDMODE);
  for (int ptcl=0; ptcl<IonConfigs[config].size(); ptcl++) 
    for (int slice=0; slice<NumTimeSlices(); slice++)
      SetPos(slice, ptcl+first, IonConfigs[config](ptcl));

  SetMode(NEWMODE);
  for (int ptcl=0; ptcl<IonConfigs[config].size(); ptcl++) 
    for (int slice=0; slice<NumTimeSlices(); slice++)
      SetPos(slice, ptcl+first, IonConfigs[config](ptcl));
}

void PathClass::Read (IOSectionClass &inSection)
{
  int myProc = Communicator.MyProc();

  /// Defaults
  Equilibrate = 1;
  CenterOfMass = 0.0;
  StoreNodeDist = false;
  StoreNodeDet = false;
  UseNodeImportance = false;
  kVecsSetup = false;
  LongRange = false;
  DavidLongRange = false;
  IlkkaLongRange = false;
  SetupPermFirstTime = true;

  SetMode(OLDMODE);
  SignWeight = 1;
  Sign = 1;
  NowOpen = false;
  HeadSlice = 0;

  SetMode(NEWMODE);
  SignWeight = 1;
  Sign = 1;
  NowOpen = false;
  HeadSlice = 0;

  // Coupling
  if (!inSection.ReadVar("FunnyCoupling",FunnyCoupling))
    FunnyCoupling=false;
  double tempExistsCoupling;
  if (!inSection.ReadVar("ExistsCoupling",tempExistsCoupling))
    ExistsCoupling=-1.0;
  else
    ExistsCoupling=tempExistsCoupling;

  // Order N
  OrderN=false;
  inSection.ReadVar("OrderN",OrderN);

  // Open Loops / Worm
  if (!inSection.ReadVar("OpenLoops",OpenPaths))
    OpenPaths=false;
  if (!inSection.ReadVar("WormOn",WormOn))
    WormOn=false;
#ifndef OPEN_LOOPS
  if (OpenPaths) {
    cerr << "OpenPaths are not enabled in the code!\n" << "Reconfigure with --enable-open and \n" << "  make clean; make ...\nAborting.\n";
    abort();
  }
#endif

  // Read system information
  assert(inSection.ReadVar ("NumTimeSlices", TotalNumSlices));
  assert(inSection.ReadVar ("tau", tau));
  Array<double,1> tempBox;
  Array<bool,1> tempPeriodic;
  TinyVector<bool,NDIM> periodic;
  assert (inSection.ReadVar ("IsPeriodic", tempPeriodic));
  assert (tempPeriodic.size() == NDIM);
  bool needBox = false;
  for (int i=0; i<NDIM; i++) {
    needBox = needBox || tempPeriodic(i);
    periodic[i] = tempPeriodic(i);
  }
  SetPeriodic (periodic);
  double desiredDensity;
  bool useDensity=inSection.ReadVar("Density",desiredDensity);
  bool useBox=inSection.ReadVar("Box",tempBox);
  assert(useBox || useDensity);
  if (needBox && !useDensity) {
    verr << "Using periodic boundary conditions.\n";
    assert(tempBox.size()==NDIM);
    for (int counter=0;counter<tempBox.size();counter++)
      Box(counter)=tempBox(counter);
    SetBox (Box);
  } else if (useDensity) {
    inSection.OpenSection("Particles");
    int numSpecSections=inSection.CountSections("Species");
    int totalNumParticles=0;
    for (int specNumber=0;specNumber<numSpecSections;specNumber++){
      inSection.OpenSection("Species",specNumber);
      int specNumParticles;
      assert(inSection.ReadVar("NumParticles",specNumParticles));
      totalNumParticles += specNumParticles;
      inSection.CloseSection();
    }
    inSection.CloseSection();
    double currentBoxVolume=tempBox(0)*tempBox(1)*tempBox(2);
    double desiredBoxVolume=1.0/desiredDensity*totalNumParticles;
    double scaleBox=pow(desiredBoxVolume/currentBoxVolume,1.0/3.0);
    assert(tempBox.size()==NDIM);
    for (int counter=0;counter<tempBox.size();counter++)
      Box(counter)=tempBox(counter)*scaleBox;
    verr<<"The BOX I am setting is "<<Box<<endl;
    SetBox (Box);
  } else {
    assert(tempBox.size()==NDIM);
    for (int counter=0;counter<tempBox.size();counter++)
      Box(counter)=tempBox(counter);
    SetBox (Box);
    verr << "Using free boundary conditions.\n";
  }

  // Read particle information
  assert(inSection.OpenSection("Particles"));
  int numSpecies = inSection.CountSections("Species");
  for (int Species=0; Species < numSpecies; Species++) {
    inSection.OpenSection("Species", Species);
    SpeciesClass *newSpecies = ReadSpecies(inSection);
    inSection.CloseSection();
    AddSpecies (newSpecies);
  }
  inSection.CloseSection(); // Particles

  // Now actually allocate the path
  Allocate();

  /// Read to see if we are using correlated sampling
  string ionSpecies;
  if (inSection.ReadVar ("IonSpecies", ionSpecies)) {
    CorrelatedSampling = true;
    IonSpecies = SpeciesNum(ionSpecies);
    assert(IonSpecies!=-1);
    Array<double,2> ionConfigA, ionConfigB;
    assert(inSection.ReadVar("IonConfigA", ionConfigA));
    assert(inSection.ReadVar("IonConfigB", ionConfigB));
    assert(ionConfigA.extent(0) == ionConfigB.extent(0));
    assert(ionConfigA.extent(1) == NDIM);
    assert(ionConfigB.extent(1) == NDIM);
    IonConfigs[0].resize(ionConfigA.extent(0));
    IonConfigs[1].resize(ionConfigB.extent(0));
    for (int i=0; i<ionConfigA.extent(0); i++)
      for (int j=0; j<3; j++) {
        IonConfigs[0](i)[j] = ionConfigA(i,j);
        IonConfigs[1](i)[j] = ionConfigB(i,j);
      }
    cerr << "IonConfigs[0] = " << IonConfigs[0] << endl;
    cerr << "IonConfigs[1] = " << IonConfigs[1] << endl;
    SetIonConfig(0);
  }

  /// Checking rounding mode
  bool roundOkay = true;
  for (int i=0; i<1000; i++) {
    double x = 10.0*drand48()-5.0;
    if (nearbyint(x) != round(x))
      roundOkay = false;
  }
  if (!roundOkay) {
    cerr << "Rounding mode is not set to ""round"".  Aborting!\n";
    abort();
  }

}


inline bool Include(dVec k)
{
  //  assert (NDIM == 3);
  if (k[0] > 0.0)
    return true;
  else if ((k[0]==0.0) && (k[1]>0.0))
    return true;
  else if ((NDIM==3) && ((k[0]==0.0) && (k[1]==0.0) && (k[2] > 0.0)))
    return true;
  else
    return false;
}


void PathClass::Allocate()
{

  SetMode(NEWMODE);

  assert(TotalNumSlices>0);
  int myProc = Communicator.MyProc();
  int numProcs = Communicator.NumProcs();
  /// Everybody gets the same number of time slices if possible.
  /// Otherwise the earlier processors get the extra one slice  
  /// until we run out of extra slices. The last slice on processor i
  /// is the first slice on processor i+1.
  MyNumSlices=TotalNumSlices/numProcs+1+(myProc<(TotalNumSlices % numProcs));

  // Initialize reference slice position to be 0 on processor 0.
  RefSlice = 0;

  /// Set the particle range for the new species
  int numParticles = 0;
  stringstream SysStr;
  SysStr << "N: ";
  for (int speciesNum=0; speciesNum<SpeciesArray.size(); speciesNum++) {
    int N = SpeciesArray(speciesNum)->NumParticles;
    int first = numParticles;
    string name = SpeciesArray(speciesNum)->Name;
    SpeciesArray(speciesNum)->FirstPtcl = first;
    numParticles=numParticles + N;
    SpeciesArray(speciesNum)->LastPtcl = first + N-1;;
    SpeciesArray(speciesNum)->Ptcls.resize(N);

    for (int i=0; i < N; i++)
      SpeciesArray(speciesNum)->Ptcls(i) = i+first;

    SysStr<<N<<" ("<<name<<"), ";
  }
  SysStr << "M: " << TotalNumSlices << ", tau: " << tau << endl;
  if (COMM::WorldProc() == 0)
    cout << SysStr.str();

  if (WormOn)
    numParticles += 2;
  Path.resize(MyNumSlices,numParticles+OpenPaths);
  ParticleExist.resize(MyNumSlices,numParticles+OpenPaths);
  RefPath.resize(numParticles+OpenPaths);
  Permutation.resize(numParticles+OpenPaths);

  SpeciesNumber.resize(numParticles+OpenPaths);
  if (WormOn)
    SpeciesNumber = 1;///HACK! HACK HACK!
  DoPtcl.resize(numParticles+OpenPaths);

  /// Assign the species number to the SpeciesNumber array
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++) {
    for (int i=SpeciesArray(speciesNum)->FirstPtcl; i<= SpeciesArray(speciesNum)->LastPtcl; i++)
      SpeciesNumber(i) = speciesNum;
  }

  // Set to the identity permutaiton
  for (int ptcl=0;ptcl<Permutation.size();ptcl++)
    Permutation(ptcl) = ptcl;

  //  InitializeJosephsonCode();
}


void PathClass::SetupkVecs(IOSectionClass& in)
{
  // Only setup kVecs if not done so already
  if (!kVecsSetup) {

    // Read in the k-space radius.
    assert(in.ReadVar("kCutoff", kCutoff));

    // Do setup
#if NDIM==2
    SetupkVecs2D();
#endif
#if NDIM==3
    SetupkVecs3D();
#endif
    if (DavidLongRange || IlkkaLongRange)
      SortRhoK();
    Rho_k.resize(MyNumSlices, NumSpecies(), kVecs.size());
  }
  UpdateRho_ks();
  bool kVecsSetup = true;
}


void PathClass::SetupkVecs2D()
{
  assert (NDIM == 2);
  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box[i];

  int numVecs=0;
  for (int i=0;i<NDIM;i++){
    MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    // perr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
  }

  dVec k;
  TinyVector<int,NDIM> ki;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      if ((dot(k,k)<kCutoff*kCutoff) && Include(k))
        numVecs++;
    }
  }
  kIndices.resize(numVecs);
  verr << "kCutoff = " << kCutoff << endl;
  verr << "Number of kVecs = " << numVecs << endl;
  kVecs.resize(numVecs);
  for (int i=0; i<NDIM; i++)
    C[i].resize(2*MaxkIndex[i]+1);
  numVecs = 0;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    ki[0]= ix+MaxkIndex[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      ki[1]= iy+MaxkIndex[1];
      if ((dot(k,k)<kCutoff*kCutoff) && Include(k)) {
        //perr<<"This k vec is "<<k[0]<<" "<<k[1]<<endl;
        kVecs(numVecs) = k;
        kIndices(numVecs)=ki;
        numVecs++;
      }
    }
  }
  SortRhoK();
  ///Now it's conceivable that we might want to add some k vectors to be additionally computed
}


void PathClass::SetupkVecs3D()
{
  assert (NDIM == 3);
  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box[i];

  int numVecs=0;
  for (int i=0;i<NDIM;i++){
    MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    // perr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
  }

  dVec k;
  TinyVector<int,NDIM> ki;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
        k[2] = iz*kBox[2];
        if ((dot(k,k)<kCutoff*kCutoff) && Include(k))
          numVecs++;
      }
    }
  }
  kIndices.resize(numVecs);
  verr << "kCutoff = " << kCutoff << endl;
  verr << "Number of kVecs = " << numVecs << endl;
  verr << "MaxkIndex = " << MaxkIndex << endl;
  kVecs.resize(numVecs);
  MagK.resize(numVecs);
  for (int i=0; i<NDIM; i++)
    C[i].resize(2*MaxkIndex[i]+1);
  numVecs = 0;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    ki[0]= ix+MaxkIndex[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      ki[1]= iy+MaxkIndex[1];
      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
        k[2] = iz*kBox[2];
        ki[2]= iz+MaxkIndex[2];
        if ((dot(k,k)<kCutoff*kCutoff) && Include(k)) {
          kVecs(numVecs) = k;
          kIndices(numVecs) = ki;
          MagK(numVecs) = sqrt(dot(k,k));
          numVecs++;
        }
      }
    }
  }
}


///Puts in the vector MagKInt a number that corresponds to the sorted
///order of the magnitude of the k vectors. Doesn't actually change
///what's in kVec. Currently only used for the reading of David's long
///range class.
void PathClass::SortRhoK()
{
  MagK.resize(kVecs.size());
  MagKint.resize(kVecs.size());
  for (int kVec=0;kVec<kVecs.size();kVec++){
    MagK(kVec)=sqrt(dot(kVecs(kVec),kVecs(kVec)));
  }
  int smallIndex=0;

  for (int counter=0;counter<MagKint.size();counter++){
    MagKint(counter)=-1;
  }

  for (int currentNum=0;currentNum<kVecs.size();currentNum++){
    for (int counter=0;counter<MagK.size();counter++){
      if (MagKint(counter)==-1 && MagK(counter)<MagK(smallIndex)){
        smallIndex=counter;
      }
    }
    for (int counter=0;counter<MagK.size();counter++){
      if (abs(MagK(smallIndex)-MagK(counter))<1e-4){
        MagKint(counter)=currentNum;
      }
    }
    smallIndex=-1;
    for (int counter=0;counter<MagK.size();counter++){
      if (MagKint(counter)==-1)
        smallIndex=counter;
    }
    if (smallIndex==-1)
      currentNum=kVecs.size()+1;
  }

}


void PathClass::CalcRho_ks_Slow(int slice, int species,Array<dVec,1> &thekvecs, Array<complex<double> ,3> rho_k)
{
  for (int ki=0; ki<thekvecs.size(); ki++) {
    complex<double> rho;
    rho = 0.0;
    for (int ptcl=Species(species).FirstPtcl; ptcl <= Species(species).LastPtcl; ptcl++) {
      const dVec &r = (*this)(slice, ptcl);
      double phase = dot(r, thekvecs(ki));
      rho += complex<double> (cos(phase), sin(phase));
    }
    rho_k(slice, species, ki) = rho;
  }
}


void PathClass::CalcRho_ks_Slow(int slice, int species)
{
  for (int ki=0; ki<kVecs.size(); ki++) {
    complex<double> rho;
    rho = 0.0;
    for (int ptcl=Species(species).FirstPtcl; ptcl <= Species(species).LastPtcl; ptcl++) {
      const dVec &r = (*this)(slice, ptcl);
      double phase = dot(r, kVecs(ki));
      rho += complex<double> (cos(phase), sin(phase));
    }
    Rho_k(slice, species, ki) = rho;
  }
}


void PathClass::CalcRho_ks_Fast(int slice,int species)
{
  //  perr<<"Beginning the calcrhok stuff"<<endl;
  // Zero out Rho_k array
  for (int ki=0;ki<kIndices.size();ki++)
    Rho_k(slice,species,ki)=0.0;

  for (int ptcl=Species(species).FirstPtcl;  ptcl <= Species(species).LastPtcl; ptcl++) {
    // First, compute C arrays
    const dVec &r = (*this)(slice, ptcl);
    for (int dim=0;dim<NDIM;dim++){
      complex<double> tempC;
      double phi = r[dim] * kBox[dim];
      tempC=complex<double>(cos(phi), sin(phi));
      C[dim](MaxkIndex[dim]) = 1.0;
      for (int n=1; n<=MaxkIndex[dim]; n++) {
        C[dim](MaxkIndex[dim]+n) = tempC * C[dim](MaxkIndex[dim]+n-1);
        C[dim](MaxkIndex[dim]-n) = conj(C[dim](MaxkIndex[dim]+n));
      }
    }
    // Now, loop over k-vector indices;
    for (int ki=0; ki<kIndices.size(); ki++) {
      const TinyVector<int,NDIM> &kIndex = kIndices(ki);
#if NDIM==3
      Rho_k(slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
#if NDIM==2
      Rho_k(slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1]);
#endif
    }
  }
}


void PathClass::UpdateRho_ks(int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  int skip = 1<<level;
  /// First, reset new to old;
  for (int slice=slice1; slice<=slice2; slice+=skip)
    for (int species=0; species<NumSpecies(); species++)
      for (int ki=0; ki<kIndices.size(); ki++)
        Rho_k[NEWMODE](slice,species,ki) = Rho_k[OLDMODE](slice,species,ki);

  /// Now update the new depending on the changed positions.
  for (int pi=0; pi<changedParticles.size(); pi++) {
    int ptcl = changedParticles(pi);
    int species = ParticleSpeciesNum(ptcl);
    for (int slice=slice1; slice <= slice2; slice+=skip) {
      dVec &rnew = Path[NEWMODE](slice,ptcl);
      dVec &rold = Path[OLDMODE](slice,ptcl);
      for (int dim=0;dim<NDIM;dim++){
        complex<double> tempC;
        double phi = rnew[dim] * kBox[dim];
        tempC=complex<double>(cos(phi), sin(phi));
        C[dim](MaxkIndex[dim]) = 1.0;
        for (int n=1; n<=MaxkIndex[dim]; n++) {
          C[dim](MaxkIndex[dim]+n) = tempC * C[dim](MaxkIndex[dim]+n-1);
          C[dim](MaxkIndex[dim]-n) = conj(C[dim](MaxkIndex[dim]+n));
        }
      }
      // Now, loop over k-vector indices;
      for (int ki=0; ki<kIndices.size(); ki++) {
        const TinyVector<int,NDIM> &kIndex = kIndices(ki);
#if NDIM==3
        Rho_k[NEWMODE](slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
#if NDIM==2
        Rho_k[NEWMODE](slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1]);
#endif
      }
      for (int dim=0;dim<NDIM;dim++){
        complex<double> tempC;
        double phi = rold[dim] * kBox[dim];
        tempC=complex<double>(cos(phi), sin(phi));
        C[dim](MaxkIndex[dim]) = 1.0;
        for (int n=1; n<=MaxkIndex[dim]; n++) {
          C[dim](MaxkIndex[dim]+n) = tempC * C[dim](MaxkIndex[dim]+n-1);
          C[dim](MaxkIndex[dim]-n) = conj(C[dim](MaxkIndex[dim]+n));
        }
      }
      // Now, loop over k-vector indices;
      for (int ki=0; ki<kIndices.size(); ki++) {
        const TinyVector<int,NDIM> &kIndex = kIndices(ki);
#if NDIM==3
        Rho_k[NEWMODE](slice,species,ki) -= C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
#if NDIM==2
        Rho_k[NEWMODE](slice,species,ki) -= C[0](kIndex[0])*C[1](kIndex[1]);
#endif
      }
    }
  }
}


void PathClass::UpdateRho_ks()
{
  ModeType mode = GetMode();
  SetMode(OLDMODE);
  for (int slice=0; slice<NumTimeSlices(); slice++)
    for (int species=0; species<NumSpecies(); species++)
      CalcRho_ks_Fast(slice,species);
  Rho_k[1] = Rho_k[0];
//   SetMode(NEWMODE);
//   for (int slice=0; slice<NumTimeSlices(); slice++)
//     for (int species=0; species<NumSpecies(); species++)
//       CalcRho_ks_Fast(slice,species);
  SetMode(mode);
}


void PathClass::MoveJoin(int oldJoin, int newJoin)
{
  //  cerr<<"IN move join"<<endl;
  //  perr<<"My time slices is "<<NumTimeSlices()<<endl;
  //  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    //    perr<<Permutation(ptcl)<<endl;
  //  }
  //  perr<<oldJoin<<" "<<newJoin<<" "<<endl;
  //  perr<<"Starting"<<endl;
  //  perr<<Path(OpenLink,OpenPtcl)<<endl;
  bool swappedAlready=false;
  int numSpecies=NumSpecies();
  //  cerr<<"A MOve join"<<endl;
  if (newJoin>oldJoin){
    //    cerr<<"Here I am "<<oldJoin<<" "<<newJoin<<" "<<" "<<endl;
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      //      cerr<<"Slice "<<timeSlice<<endl;
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
        // cerr<<"ptcl "<<ptcl<<" "<<Permutation(ptcl)<<endl;
        Path[OLDMODE](timeSlice,ptcl)=Path[NEWMODE](timeSlice,Permutation(ptcl));
        if (OpenPaths) {
          // cerr<<"IN OPEN PATHS"<<endl;
          if (timeSlice==(int)OpenLink && Permutation(ptcl)==(int)OpenPtcl && !swappedAlready){
            OpenPtcl[OLDMODE]=ptcl;
            OpenPtcl[NEWMODE]=ptcl;
            swappedAlready=true;
          }
        }
      }
    }
    //    cerr<<"B MOve join"<<endl;
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++) {
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
        Path[NEWMODE](timeSlice,ptcl)=Path[OLDMODE](timeSlice,ptcl);
    }
  }
  else if (oldJoin>newJoin){
    //  cerr<<"C MOve join"<<endl;
    //  else if (oldJoin>=newJoin){//CHANGED!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
        if (OpenPaths)
          if (timeSlice==(int)OpenLink && ptcl==(int)OpenPtcl && !swappedAlready){
            OpenPtcl[OLDMODE]=Permutation(ptcl);
            OpenPtcl[NEWMODE]=Permutation(ptcl);
            swappedAlready=true;
          }
        Path[OLDMODE](timeSlice,Permutation(ptcl))=Path[NEWMODE](timeSlice,ptcl);
      }
    }
    //    cerr<<"D MOve join"<<endl;
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
        Path[NEWMODE](timeSlice,ptcl)=Path[OLDMODE](timeSlice,ptcl);
    }
  }
  //  perr<<Path(OpenLink,OpenPtcl)<<endl;
  //  perr<<"Ending"<<endl;
  //  cerr<<"E MOve join"<<endl;
  if (WormOn)
    MoveJoinParticleExist(oldJoin,newJoin);
  //  cerr<<"out move join"<<endl;
}


void PathClass::AcceptCopy(int startSlice,int endSlice, const Array <int,1> &activeParticles)
{
  // Accept Weights
  Sign.AcceptCopy();
  ExistsCoupling.AcceptCopy();
  SignWeight.AcceptCopy();

  cm2 = 0.0;
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    // HACK!! Though not really, I commented this out b/c I believe it's not necessary
    // since we make sure to PutInBox everytime we measure distances. Additionally
    // this messes up centroid calculations. -ewb
    //for (int slice=startSlice;slice<=endSlice;slice++)
    //  PutInBoxFast(Path(slice,ptcl));

    // Calculate center of mass movement
    dVec cmPart = 0.0;
    for (int slice=startSlice;slice<endSlice;slice++){
      dVec newLoc = Path[NEWMODE](slice,ptcl);
      dVec oldLoc = Path[OLDMODE](slice,ptcl);
      dVec minDisp = MinImageDisp(oldLoc,newLoc);
      CenterOfMass = CenterOfMass + minDisp;
      cmPart = cmPart + minDisp;
    }
    cm2 += dot(cmPart,cmPart)/TotalNumSlices;

    // Accept new coordinates and permutation
    Path[OLDMODE](Range(startSlice, endSlice), ptcl) = Path[NEWMODE](Range(startSlice, endSlice), ptcl);
    if (WormOn)
      ParticleExist[OLDMODE](Range(startSlice, endSlice), ptcl) = ParticleExist[NEWMODE](Range(startSlice, endSlice), ptcl);
    Permutation.AcceptCopy(ptcl);
  }

  // Accept new rho_k
  Rho_k[OLDMODE](Range(startSlice,endSlice), Range::all(), Range::all()) =
    Rho_k[NEWMODE](Range(startSlice,endSlice), Range::all(), Range::all());

  // Accept new nodal distances
  if (StoreNodeDist)
    for (int species=0; species<NumSpecies(); species++)
      NodeDist[OLDMODE](Range(startSlice,endSlice), species) = NodeDist[NEWMODE](Range(startSlice,endSlice), species);
  if (StoreNodeDet)
    for (int species=0; species<NumSpecies(); species++)
      NodeDet[OLDMODE](Range(startSlice,endSlice), species) = NodeDet[NEWMODE](Range(startSlice,endSlice), species);

  // Accept open paths
  if (OpenPaths) {
    NowOpen.AcceptCopy();
    OpenPtcl.AcceptCopy();
    OpenLink.AcceptCopy();

    for (int counter=0; counter<NumTimeSlices(); counter++)
      Path[OLDMODE](counter,NumParticles()) = Path[NEWMODE](counter,NumParticles());

    for (int ptcl=0;ptcl<NumParticles();ptcl++) {
      Path[OLDMODE](OpenLink[NEWMODE],ptcl) = Path[NEWMODE](OpenLink[NEWMODE],ptcl);
      Path[OLDMODE](0,ptcl) = Path[NEWMODE](0,ptcl);
      Permutation.AcceptCopy(ptcl);
    }
  }

  // Regrid for O(N)
  if (OrderN) {
    for (int slice=startSlice; slice<=endSlice; slice++) {
      for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
        int ptcl = activeParticles(ptclIndex);
        Cell.ReGrid(slice,ptcl);
        //      Cell.BinParticles(slice);
      }
    }
  }
}


void PathClass::RejectCopy(int startSlice,int endSlice, const Array <int,1> &activeParticles)
{
  // Reject Weights
  Sign.RejectCopy();
  ExistsCoupling.RejectCopy();
  SignWeight.RejectCopy();

  // Reject particle coordinates and permutation
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[NEWMODE](Range(startSlice, endSlice), ptcl) = Path[OLDMODE](Range(startSlice, endSlice), ptcl);
    if (WormOn)
      ParticleExist[NEWMODE](Range(startSlice, endSlice), ptcl) = ParticleExist[OLDMODE](Range(startSlice, endSlice), ptcl);
    Permutation.RejectCopy(ptcl);
  }

  // Reject rho_k
  Rho_k[NEWMODE](Range(startSlice,endSlice), Range::all(), Range::all()) =
    Rho_k[OLDMODE](Range(startSlice,endSlice), Range::all(), Range::all());

  // Reject nodal distances
  if (StoreNodeDist)
    for (int species=0; species<NumSpecies(); species++)
      NodeDist[NEWMODE](Range(startSlice,endSlice), species) = NodeDist[OLDMODE](Range(startSlice,endSlice), species);
  if (StoreNodeDet)
    for (int species=0; species<NumSpecies(); species++)
      NodeDet[NEWMODE](Range(startSlice,endSlice), species) = NodeDet[OLDMODE](Range(startSlice,endSlice), species);

  // Reject open paths
  if (OpenPaths) {
    NowOpen.RejectCopy();
    OpenPtcl.RejectCopy();
    OpenLink.RejectCopy();

    Path[NEWMODE](Range(startSlice,endSlice),NumParticles()) = Path[OLDMODE](Range(startSlice,endSlice),NumParticles());
    for (int ptcl=0;ptcl<NumParticles();ptcl++) {
      Path[NEWMODE](OpenLink[NEWMODE],ptcl) = Path[OLDMODE](OpenLink[NEWMODE],ptcl);
      Path[NEWMODE](0,ptcl) = Path[OLDMODE](0,ptcl);
      Permutation.RejectCopy(ptcl);
    }

  }

}


void PathClass::ShiftData(int slicesToShift)
{

  if (WormOn)
    ShiftParticleExist(slicesToShift);
  ShiftPathData(slicesToShift);
  if (StoreNodeDist)
    ShiftNodeDist(slicesToShift);
  if (StoreNodeDet)
    ShiftNodeDet(slicesToShift);
  if (LongRange) {
    // ShiftRho_kData(slicesToShift);
    UpdateRho_ks();
  }
  OpenLink.AcceptCopy(); ///the open link has changed and you want to accept it
  RefSlice += slicesToShift;
  while (RefSlice >= TotalNumSlices)
    RefSlice -= TotalNumSlices;
  while (RefSlice < 0)
    RefSlice += TotalNumSlices;
  //  perr<<"My ref slice at particle 0 is "<<Path(RefSlice,0)<<endl;
  if (OpenPaths){

    //    perr<<"Here my open link is "<<OpenLink<<endl;
    int openLinkOld=(int)OpenLink;
    int startSliceOpenLinkProc;
    int endSliceOpenLinkProc;
    SliceRange(SliceOwner(RefSlice),startSliceOpenLinkProc,endSliceOpenLinkProc);
    OpenLink=RefSlice-startSliceOpenLinkProc;
    if ((int)OpenLink==0){
      OpenLink=NumTimeSlices()-1;
    }
    //  perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    ////    perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl;
    //    if (openLinkOld!=(int)OpenLink)
    //      cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl<<Communicator.MyProc()<<" "<<SliceOwner(RefSlice)<<endl;
    //    assert(openLinkOld==(int)OpenLink);
    OpenLink[OLDMODE]=OpenLink[NEWMODE];
  }

}


void PathClass::ShiftRho_kData(int slicesToShift)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numSlices  = Rho_k.extent(0);
  int numSpecies = Rho_k.extent(1);
  int numk       = Rho_k.extent(2);
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  /// First shifts the data in the A copy left or right by the
  /// appropriate amount 
  if (slicesToShift>0) {
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int species=0; species<numSpecies; species++)
        for (int ki=0; ki<numk; ki++)
          Rho_k[0](slice,species,ki)=Rho_k[0](slice-slicesToShift,species,ki);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int species=0; species<numSpecies; species++)
        for (int ki=0; ki<numk; ki++)
          Rho_k[0](slice,species,ki)=Rho_k[0](slice-slicesToShift,species,ki);
  }

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numSpecies*numk;
  Array<complex<double>,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0) {
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
      for (int species=0; species<numSpecies; species++)
        for (int ki=0; ki<numk; ki++) {
          /// If shifting forward, don't send the last time slice (so always)
          /// send slice-1
          sendBuffer(buffIndex)=Rho_k[1](slice-1,species,ki);
          buffIndex++;
        }
  } else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
      for (int species=0; species<numSpecies; species++)
        for (int ki=0; ki<numk; ki++) {
          /// If shifting backward, don't send the first time slice (so always)
          /// send slice+1
          sendBuffer(buffIndex)=Rho_k[1](slice+1,species,ki);
          buffIndex++;
        }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);

  if (slicesToShift>0)
    startSlice=0;
  else
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
    for (int species=0; species<numSpecies; species++)
      for (int ki=0; ki<numk; ki++) {
        Rho_k[0](slice,species,ki)=receiveBuffer(buffIndex);
        buffIndex++;
      }

  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int species=0; species<numSpecies; species++)
      for (int ki=0; ki<numk; ki++)
        Rho_k[1](slice,species,ki) = Rho_k[0](slice,species,ki);

  // And we're done!
}


void PathClass::ShiftNodeDet(int slicesToShift)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numSpecies=NumSpecies();
  int numSlices=NumTimeSlices();
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int species=0;species<numSpecies;species++)
        NodeDet[NEWMODE](slice,species) = NodeDet[NEWMODE](slice-slicesToShift,species);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int species=0;species<numSpecies;species++)
        NodeDet[NEWMODE](slice,species) = NodeDet[NEWMODE](slice-slicesToShift,species);
  }

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numSpecies;
  Array<double,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0){
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int species=0;species<numSpecies;species++) {
        ///If shifting forward, don't send the last time slice (so always)
        ///send slice-1
        sendBuffer(buffIndex)=NodeDet[OLDMODE](slice-1,species);
        buffIndex++;
      }
    }
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int species=0;species<numSpecies;species++) {
        ///If shifting backward, don't send the first time slice (so always)
        ///send slice+1
        sendBuffer(buffIndex)=NodeDet[OLDMODE](slice+1,species);
        buffIndex++;
      }
    }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);

  if (slicesToShift>0)
    startSlice=0;
  else
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
    for (int species=0;species<numSpecies;species++){
      NodeDet[NEWMODE](slice,species)=receiveBuffer(buffIndex);
      buffIndex++;
    }
  }

  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int species=0;species<numSpecies;species++)
      NodeDet[OLDMODE](slice,species) = NodeDet[NEWMODE](slice,species);

}


void PathClass::ShiftNodeDist(int slicesToShift)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numSpecies=NumSpecies();
  int numSlices=NumTimeSlices();
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int species=0;species<numSpecies;species++)
        NodeDist[NEWMODE](slice,species) = NodeDist[NEWMODE](slice-slicesToShift,species);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int species=0;species<numSpecies;species++)
        NodeDist[NEWMODE](slice,species) = NodeDist[NEWMODE](slice-slicesToShift,species);
  }

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numSpecies;
  Array<double,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0){
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int species=0;species<numSpecies;species++) {
        ///If shifting forward, don't send the last time slice (so always)
        ///send slice-1
        sendBuffer(buffIndex)=NodeDist[OLDMODE](slice-1,species);
        buffIndex++;
      }
    }
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int species=0;species<numSpecies;species++) {
        ///If shifting backward, don't send the first time slice (so always)
        ///send slice+1
        sendBuffer(buffIndex)=NodeDist[OLDMODE](slice+1,species);
        buffIndex++;
      }
    }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);

  if (slicesToShift>0)
    startSlice=0;
  else
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
    for (int species=0;species<numSpecies;species++){
      NodeDist[NEWMODE](slice,species)=receiveBuffer(buffIndex);
      buffIndex++;
    }
  }

  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int species=0;species<numSpecies;species++)
      NodeDist[OLDMODE](slice,species) = NodeDist[NEWMODE](slice,species);

}


void PathClass::ShiftPathData(int slicesToShift)
{
  //  perr<<"Slices to shift are "<<slicesToShift<<endl;
  //  perr<<"I'm in shiftpathdata with numparticles "<<NumParticles()+OpenPaths<<"and Openpaths being "<<OpenPaths<<"and slicesToShift is "<<slicesToShift<<endl;
  //  sleep(10);
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numPtcls=NumParticles()+OpenPaths;
  int numSlices=NumTimeSlices();
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
        Path[NEWMODE](slice,ptcl) = Path[NEWMODE](slice-slicesToShift,ptcl);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
        Path[NEWMODE](slice,ptcl) = Path[NEWMODE](slice-slicesToShift,ptcl);
  }

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numPtcls;
  Array<dVec,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0){
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
        ///If shifting forward, don't send the last time slice (so always)
        ///send slice-1
        sendBuffer(buffIndex)=Path[OLDMODE](slice-1,ptcl);
        buffIndex++;
      }
    }
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
        ///If shifting backward, don't send the first time slice (so always)
        ///send slice+1
        sendBuffer(buffIndex)=Path[OLDMODE](slice+1,ptcl);
        buffIndex++;
      }
    }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);

  if (slicesToShift>0)
    startSlice=0;
  else
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
    for (int ptcl=0;ptcl<numPtcls;ptcl++){
      Path[NEWMODE](slice,ptcl)=receiveBuffer(buffIndex);
      buffIndex++;
    }
  }

  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int ptcl=0; ptcl<numPtcls; ptcl++)
      Path[OLDMODE](slice,ptcl) = Path[NEWMODE](slice,ptcl);

  // And we're done!
  if (OpenPaths){ //only works for serialo positive slices
    if ((int)OpenLink+slicesToShift >= NumTimeSlices())
      OpenLink=OpenLink+1;
    if ((int)OpenLink+slicesToShift <0)
      OpenLink=OpenLink-1;
    OpenLink=((int)OpenLink+slicesToShift+NumTimeSlices()) % NumTimeSlices();
    if ((int)OpenLink==0){
      perr<<"equal to 0"<<endl;
      OpenLink=NumTimeSlices()-1;
    }

    //    OpenLink=RefSlice;
    //    if ((int)OpenLink==0){
    //      OpenLink=NumTimeSlices()-1;
    //    }
    //    perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    if (openLinkOld!=(int)OpenLink)
    //      perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    assert(openLinkOld==(int)OpenLink);
    OpenLink[OLDMODE]=OpenLink[NEWMODE];
  }
  //  perr<<Path(OpenLink,NumParticles())<<endl;
  //  perr<<"I leave shiftpathdata"<<endl;
  //  sleep(10);

}


/// This function sends the reference slice from whatever processor is
/// currently holding it to all the other processors.  Note that only
/// the NEW copy is updated.  If we want the OLD copy to be updated,
/// we must use RefPath.AcceptCopy().
void PathClass::BroadcastRefPath()
{
  SetMode (NEWMODE);
  // Figure which processor has the reference slice
  int procWithRefSlice = SliceOwner (RefSlice);

  Array<dVec,1> buffer(NumParticles());
  if (procWithRefSlice == Communicator.MyProc()) {
    int myStart, myEnd;
    SliceRange (Communicator.MyProc(), myStart, myEnd);
    for (int ptcl=0; ptcl<NumParticles(); ptcl++)
      buffer(ptcl) = Path(RefSlice-myStart, ptcl);
  }

  //  perr << "procWithRefSlice = " << procWithRefSlice << endl;

  // Do broadcast
  Communicator.Broadcast(procWithRefSlice, buffer);

  // Now, all processors have reference slice.  Note that the 
  // labeling of particles may not be right since with haven't
  // propogated the permuation matrices.

  // Now, copy into path
  for (int ptcl=0; ptcl<NumParticles(); ptcl++)
    RefPath(ptcl) = buffer(ptcl);
}


// Determine sign of a given path. Total sign will have to multiply
// this number for all paths.
int PathClass::GetSign()
{
  int N = NumParticles();
  int ptcl = 0;
  Array<bool,1> CountedAlready(N);
  CountedAlready = false;
  int sign = 1;
  while (ptcl < N) {
    int species = ParticleSpeciesNum(ptcl);
    if (Species(species).GetParticleType() == 0 && !CountedAlready(ptcl)) {
      int startPtcl = ptcl;
      int roamingPtcl = ptcl;
      int cycleLength = 0;
      roamingPtcl = Permutation(roamingPtcl);
      while (roamingPtcl != startPtcl) {
        CountedAlready(roamingPtcl) = true;
        cycleLength++;
        roamingPtcl = Permutation(roamingPtcl);
      }
      sign *= pow(-1,cycleLength);
    }
    ptcl++;
  }
  return sign;
}


// Combines the permutation vectors of all the processors.  Note:
// only processor 0 gets the result.  All others get junk.
void PathClass::TotalPermutation(Array<int,1> &permVec)
{
  int myProc = Communicator.MyProc();
  int numProcs = Communicator.NumProcs();
  int numPtcls = NumParticles();
  if (permVec.size()!=Permutation.size())
    permVec.resize(Permutation.size());
  if (permMat.extent(0) !=numProcs || permMat.extent(1)!=Permutation.size()){
    permMat.resize(numProcs,Permutation.size());
    //    Array<int,2> permMat(numProcs, permVec.size());
  }
  //  cerr<<"Perm sizes: "<<Permutation.size()<<" "<<numProcs<<" "<<permMat.extent(0)<<" "<<permMat.extent(1)<<" "<<permVec.size()<<endl;
  for (int i=0; i<permVec.size(); i++)
    permVec(i) = Permutation(i);

  // First, collect all the individual permutation vectors
  Communicator.Gather (permVec, permMat, 0);
  // Now apply them in sequence.
  if (Communicator.MyProc() == 0) {
    for (int pi=0; pi < numPtcls; pi++) {
      int ptcl = pi;
      for (int proc=0; proc<numProcs; proc++)
        ptcl = permMat(proc, ptcl);
      permVec(pi) = ptcl;
    }
  }
}

///Must start diagonal or this is goign to break because of where
///you are starting the openlink at 0
void PathClass::InitOpenPaths()
{
  perr<<"Starting to initialize"<<endl;
  if (OpenPaths){

    perr<<"openpaths"<<endl;
    SetMode(OLDMODE);

    if (Communicator.MyProc()==0)
      OpenLink=NumTimeSlices()-1;
    else
      OpenLink=-1;
    perr<<"Set up the open link"<<endl;
    OpenPtcl=Species(OpenSpeciesNum).FirstPtcl;
    perr<<"set up the open particle"<<endl;
    SetMode(NEWMODE);

    if (Communicator.MyProc()==0)
      OpenLink=NumTimeSlices()-1;
    else 
      OpenLink=-1;
    OpenPtcl=0;
    if (Communicator.MyProc()==0){
      perr<<"Preparing for moving things in path around"<<endl;
      Path[OLDMODE]((int)OpenLink,NumParticles())=Path[OLDMODE]((int)OpenLink,(int)OpenPtcl);
      perr<<"Moved the first thing"<<endl;
      Path[NEWMODE]((int)OpenLink,NumParticles())=Path[NEWMODE]((int)OpenLink,(int)OpenPtcl);
      perr<<"Moved the second thing"<<endl;
      Path[OLDMODE](0,NumParticles())=Path[OLDMODE](0,(int)OpenPtcl);
      perr<<"Moved the first thing"<<endl;
      Path[NEWMODE](0,NumParticles())=Path[NEWMODE](0,(int)OpenPtcl);
    }
  }
  perr<<"Initialized the open paths"<<endl;
}


void PathClass::WarpAtoB(dVec &pos)
{
  double weightSum = 0.0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[0](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    weightSum += 1.0/(dist2*dist2);
  }
  dVec shift;
  for(int i=0; i<NDIM; i++)
    shift[i] = 0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[0](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    double weight = 1.0/(dist2*dist2*weightSum);
    shift = shift + weight*(IonConfigs[1](i)-IonConfigs[0](i));
  }
  pos = pos + shift;
}


void PathClass::WarpBtoA(dVec &pos)
{
  double weightSum = 0.0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[1](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    weightSum += 1.0/(dist2*dist2);
  }
  dVec shift;
  for(int i=0; i<NDIM; i++)
    shift[i] = 0;
  for (int i=0; i<IonConfigs[1].size(); i++) {
    dVec disp = pos - IonConfigs[1](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    double weight = 1.0/(dist2*dist2*weightSum);
    shift = shift + weight*(IonConfigs[0](i)-IonConfigs[1](i));
  }
  pos = pos + shift;
}


void PathClass::WarpPaths (int ionSpecies)
{
  SpeciesClass &ions = Species(ionSpecies);
  int N = ions.NumParticles;
  Array<dVec,1> Rold(N), Rnew(N), Rdelta(N);
  for (int i=0; i<N; i++) {
    Rold(i) = Path[0](0,i+ions.FirstPtcl);
    Rnew(i) = Path[1](0,i+ions.FirstPtcl);
    Rdelta(i) = Rnew(i)-Rold(i);
    PutInBox(Rdelta(i));
  }
  SetMode(OLDMODE);
  for (int elec=0; elec<NumParticles(); elec++) {
    if (ParticleSpeciesNum(elec) != ionSpecies) {
      for (int slice=0; slice<NumTimeSlices(); slice++) {
        dVec disp;
        double dist;
        double totalWeight = 0.0;
        for (int ion=0; ion<N; ion++) {
          DistDisp(slice, elec, ion+ions.FirstPtcl, dist, disp);
          totalWeight += 1.0/(dist*dist*dist*dist);
        }
        for (int ion=0; ion<N; ion++) {
          DistDisp(slice, elec, ion+ions.FirstPtcl, dist, disp);
          double weight = 1.0/(dist*dist*dist*dist*totalWeight);
          Path[1](slice, elec) += weight*Rdelta(ion);
        }
        Path[0](slice, elec) = Path[1](slice, elec);
      }
    }
  }
}


void PathClass::DistDispFast (int sliceA, int sliceB, int ptcl1, int ptcl2, double &distA, double &distB,dVec &dispA, dVec &dispB)
{
  dispA = Path(sliceA, ptcl2) - Path(sliceA,ptcl1);
  dispB = Path(sliceB, ptcl2) - Path(sliceB,ptcl1);
  dVec preA=dispA;
  dVec preB=dispB;
  dVec boxOver2=GetBox()/2;
  dVec box=GetBox();
  for (int dim=0;dim<NDIM;dim++){
    if (dispA[dim]>boxOver2[dim])
      dispA[dim]-=box[dim];
    if (dispA[dim]<-boxOver2[dim])
      dispA[dim]+=box[dim];
  }
  dVec dispC=dispA-dispB;
  dVec m;
  dVec temptempdispB;
  m[0] = nearbyint((dispA[0]-dispB[0])*BoxInv[0]);
  m[1] = nearbyint((dispA[1]-dispB[1])*BoxInv[1]);
  temptempdispB[0] = dispB[0]+m[0]*IsPeriodic[0]*Box[0];
  temptempdispB[1] = dispB[1]+m[1]*IsPeriodic[1]*Box[1];
  for (int dim=0;dim<NDIM;dim++){
    if (dispC[dim]>boxOver2[dim])
      dispB[dim]+=box[dim];
    if (dispC[dim]<-boxOver2[dim])
      dispB[dim]-=box[dim];
  }
  bool debug=false;
  if (debug){
    dVec tempDispA;
    dVec tempDispB;
    double tempDistA;
    double tempDistB;
    DistDisp (sliceA, sliceB, ptcl1, ptcl2, tempDistA, tempDistB,tempDispA, tempDispB);
    //  if (!((tempDispA==dispA) && (tempDispB==dispB))){
    if ((abs(tempDispA(0)-dispA(0))>1e-8) ||
        (abs(tempDispA(1)-dispA(1))>1e-8) ||
        (abs(tempDispB(0)-dispB(0))>1e-8) ||
        (abs(tempDispB(1)-dispB(1))>1e-8)) {
      cerr<<tempDispA<<" "<<dispA<<" "<<tempDispB<<" "<<dispB<<endl;
      cerr<<temptempdispB<<" "<<m<<endl;
      cerr<<dispC<<endl;
      cerr<<preA<<" "<<preB<<endl;
    }
    //  assert(tempDispA==dispA);
    //  assert(tempDispB==dispB);
    assert(abs(tempDispA(0)-dispA(0))<1e-8);
    assert(abs(tempDispA(1)-dispA(1))<1e-8);
    assert(abs(tempDispB(0)-dispB(0))<1e-8);
    assert(abs(tempDispB(1)-dispB(1))<1e-8);
  }
  distA = sqrt(dot(dispA,dispA));
  distB = sqrt(dot(dispB,dispB));
  return;

  //do nothing for now

}


void PathClass::PutInBox(dVec &v,dVec &box)
{
  dVec boxInv;
  for (int dim=0;dim<NDIM;dim++)
    boxInv(dim)=1.0/box(dim);
  for (int i=0; i<NDIM; i++) {
    double n = -floor(v(i)*boxInv(i)+0.5);
    v(i) += n*IsPeriodic(i)*box(i);
  }
}


void PathClass::PutInBoxFast (dVec &v)
{
  bool debug=false;
  dVec oldV;
  if (debug)
    oldV=v;
  dVec box=GetBox();
  dVec halfBox=GetBox()/2;
  for (int dim=0;dim<NDIM;dim++){
    if (v[dim]>halfBox[dim])
      v[dim]-=box[dim];
    if (v[dim]<-halfBox[dim])
      v[dim]+=box[dim];
  }
  if (debug){
    PutInBox(oldV);
    assert(v==oldV);
  }
}


// Setup permutation sectors
void PathClass::SetupPermSectors(int n, int MaxNSectors)
{
  if (Communicator.MyProc() == 0 && SetupPermFirstTime) {
    SetupPermFirstTime = false;
    cout << CloneStr << " Setting up Permutation Sectors" << endl;
    vector<int> a;
    a.resize(n);
    for (int i=0; i<n; i++) {
      a[i] = 0;
    }
    int k = 1;
    int y = n-1;
    vector< vector<int> > tmpPossPerms;
    while (k != 0 && (MaxNSectors > PossPerms.size() || !MaxNSectors)) {
      int x = a[k-1] + 1;
      k -= 1;
      while (2*x <= y) {
        a[k] = x;
        y -= x;
        k += 1;
      }
      int l = k+1;
      while (x <= y && (MaxNSectors > PossPerms.size() || !MaxNSectors)) {
        a[k] = x;
        a[l] = y;
        vector<int> b;
        for (vector<int>::size_type j=0; j!=k+2; j++)
          b.push_back(a[j]);
        tmpPossPerms.push_back(b);
        x += 1;
        y -= 1;
      }
      a[k] = x+y;
      y = x+y-1;
      vector<int> c;
      for (vector<int>::size_type j=0; j!=k+1; j++)
        c.push_back(a[j]);
      tmpPossPerms.push_back(c);
    }


    //cout << CloneStr << " Putting sectors in map" << endl;
    int nSectors = tmpPossPerms.size();
    for (vector<int>::size_type j=0; j != nSectors; j++) {
      //sort(tmpPossPerms[j].begin(),tmpPossPerms[j].end());
      //for (int i = 0; i < tmpPossPerms[j].size(); i++)
      //  cout << tmpPossPerms[j][i] << " ";
      //cout << endl;
      PossPerms[tmpPossPerms[j]] = j;
    }
  }

}


void PathClass::GetPermInfo(vector<int> &Cycles, int &PermSector)
{
  int N = NumParticles();
  if (CountedAlready.size() != N) {
    CountedAlready.resize(N);
    TotalPerm.resize(N);
  }
  TotalPermutation(TotalPerm); /// Only proc 0 gets TotalPerm
  CountedAlready = false;
  int ptcl = 0;
  int PermNumber = 0;
  if (Communicator.MyProc() == 0) {
    while (ptcl < N) {
      if (!CountedAlready(ptcl)) {
        int startPtcl = ptcl;
        int roamingPtcl = ptcl;
        int cycleLength = 0;
        roamingPtcl = TotalPerm(roamingPtcl);
        while (roamingPtcl != startPtcl) {
          CountedAlready(roamingPtcl) = true;
          cycleLength++;
          roamingPtcl = TotalPerm(roamingPtcl);
        }
        Cycles.push_back(cycleLength+1);
        PermNumber += cycleLength;
      }
      ptcl++;
    }
  } else
    return;

  //while (!Cycles.empty()) {
  //  int cycle = Cycles.back()-1;
  //  if (ThisPerm.find(cycle) == ThisPerm.end())
  //    ThisPerm.insert(pair<int,int>(cycle,1));
  //  else
  //    ThisPerm[cycle] += 1;
  //  Cycles.pop_back();
  //}

  sort(Cycles.begin(),Cycles.end());
  PossPermsIterator = PossPerms.find(Cycles);
  if (PossPermsIterator == PossPerms.end()) {
    cerr << "Broken Permutation: " << endl;
    for (vector<int>::size_type i=0; i != Cycles.size(); i++)
      cerr << Cycles[i] << " ";
    cerr << endl;
    exit(1);
  } else {
    PermSector = PossPermsIterator->second;
    return;
  }
}


