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

#include "Communication/Communication.h"
#include "Blitz.h"
#include "Random/Random.h"
#include "IO/IO.h"

using namespace IO;
// #define NDIM 3

typedef TinyVector<int,NDIM> State;
typedef TinyVector<double,NDIM> dVec;

inline bool
equal (const State& st1, const State& st2)
{
  bool eq = true;
  for (int i=0; i<NDIM; i++)
    eq = eq && (st1[i] == st2[i]);
  return eq;
}

class ParticleClass
{
protected:
  int NumParticles;
  dVec Box;
  dVec kPrim;
  Array<State,1> AllowedStates;
  Array<complex<double>,1> Rhok;
  double Esum, E2sum;
  int BlockSize, NumBlocks;
  double lambda, beta;
  CommunicatorClass Communicator;
  RandomClass Random;
  double Energy();
public: 
  Array<State,1> OccupiedStates;
  virtual double AcceptProb (int ptclNum, State &newState) = 0;
  // Returns sampling probablilty ratio oldProb/NewProb;
  virtual double Sample(int &ptclToMove, State &newState);
  void MCStep();
  void Observe();
  void DumpResults();
  void Read (IOSectionClass &in);
  virtual void FillStates();
  void Run();
  void PrintStates();
  ParticleClass() :
    Random(Communicator)
  {
    //    Random.Init();
  }
};


class BoltzmannonClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
};

class BosonClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
};

class FermionClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
  void FillStates();
  void FillFermiSea();
};


////////////////////////////////////
// ParticleClass Member Functions //
////////////////////////////////////
void ParticleClass::Read(IOSectionClass &in)
{
  Array<double,1> box;
  assert(in.ReadVar ("Box", box));
  for (int i=0; i<NDIM; i++) {
    Box[i] = box(i);
    kPrim[i] = 2.0*M_PI/box(i);
  }
  assert(in.ReadVar ("NumParticles", NumParticles));
  assert(in.ReadVar ("beta", beta));
  assert(in.ReadVar ("BlockSize", BlockSize));
  assert(in.ReadVar ("NumBlocks", NumBlocks));
  assert(in.ReadVar ("lambda", lambda));

  int seed;
  bool haveSeed = in.ReadVar ("Seed", seed);
  // Now, set up random number generator
  if (haveSeed)
    Random.Init (seed);
  else
    Random.InitWithRandomSeed();


  Esum = 0.0;
  E2sum = 0.0;
  FillStates();
}

void ParticleClass::FillStates()
{
  OccupiedStates.resize(NumParticles);
  for (int ptcl=0; ptcl<NumParticles; ptcl++)
    OccupiedStates(ptcl) = 0;
}

void ParticleClass::MCStep()
{
  State newState;
  int ptclToMove;
 
  double sampleRatio = Sample (ptclToMove, newState);

  double toAccept = AcceptProb (ptclToMove, newState);
  if (toAccept*sampleRatio>Random.Local())
    OccupiedStates(ptclToMove)=newState;
}

double ParticleClass::Sample(int &ptclToMove, State &newState)
{
  ptclToMove=Random.LocalInt(OccupiedStates.size());

  if (Random.Local()>0.5){ //make a energy changing move
    int toChange;
    if (Random.Local()>0.5)
      toChange=1;
    else
      toChange=-1;
    int dimToChange=Random.LocalInt(NDIM);
    newState=OccupiedStates(ptclToMove);
    newState[dimToChange]+=toChange;

  }
  else {
    newState=OccupiedStates(ptclToMove);
    int dimToChange=Random.LocalInt(NDIM);
    newState[dimToChange]=newState[dimToChange]*-1;
    //    cerr<<"This is "<<ptclToMove<<" "<<OccupiedStates(ptclToMove)<<" "<<newState<<endl;
    
  }
  return 1.0;
}


double ParticleClass::Energy()
{
  double E=0.0;
  for (int ptcl=0; ptcl<NumParticles; ptcl++) {
    //    cerr << "State: " << OccupiedStates(ptcl) << endl;
    for (int i=0; i<NDIM; i++)
      E += kPrim[i]*kPrim[i]*OccupiedStates(ptcl)[i]*OccupiedStates(ptcl)[i];
  }
  E *= lambda;
  return (E);
}


void ParticleClass::Run()
{
  double meanSum, mean2Sum;
  meanSum = mean2Sum = 0.0;
  for (int block=0; block<NumBlocks; block++) {
    // Reset block sums
    Esum = 0.0;
    for (int step=0; step<BlockSize; step++) {
      MCStep();
      double E = Energy();
      Esum += E;
    }
    double E = Esum / BlockSize;
    meanSum += E;
    mean2Sum += E*E;
    //    PrintStates();
  }
  // Write block data
  double mean = meanSum / NumBlocks;
  double mean2 = mean2Sum / NumBlocks;
  double var = sqrt (mean2 - mean*mean);
  fprintf (stderr, "Energy:\n%1.12e +/- %1.12e\n", mean, sqrt(var/NumBlocks));
}


///////////////////////////////////
// FermionClass Member Functions //
///////////////////////////////////

void ParticleClass::PrintStates()
{
  cerr<<"Beginning to Print States"<<endl;
  for (int counter=0;counter<OccupiedStates.size();counter++){
    cerr<<OccupiedStates(counter)<<endl;
  }
  cerr<<"Done Printing States";


}


void FermionClass::FillStates()
{
  OccupiedStates.resize (NumParticles);
  double kmin = min(kPrim[0], kPrim[1]);
#if NDIM==3
  kmin = min (kmin, kPrim[2]);
#endif
  double deltaEmin = lambda * kmin*kmin;

  int n = 0;
  double Ecut = deltaEmin;
  while (n < NumParticles) {
    int ixmax = (int)ceil(sqrt(Ecut/lambda)/kPrim[0]);
    int iymax = (int)ceil(sqrt(Ecut/lambda)/kPrim[1]);
#if NDIM==3
    int izmax = (int)ceil(sqrt(Ecut/lambda)/kPrim[2]);
#endif
    dVec k;
    State tryState;
    for (int ix=-ixmax; ix<=ixmax; ix++) {
      tryState[0] = ix;
      k[0] = kPrim[0]*ix;
      for (int iy=-iymax; iy<=iymax; iy++) {
        tryState[1] = iy;
        k[1] = kPrim[1]*iy;
#if NDIM==3
        for (int iz=-izmax; iz<=izmax; iz++) {
          tryState[2] = iz;
          k[2] = kPrim[2]*iz;
#endif
          double E = lambda * dot(k,k);
          if (n < NumParticles)
            if (E < Ecut) {
              bool occupied = false;
              for (int i=0; i<n; i++)
                occupied = occupied || equal(tryState,OccupiedStates(i));
              if (!occupied) {
                OccupiedStates(n) = tryState;
                n++;
              }
            }
        }
        Ecut += deltaEmin;
      }
    }
  }
}


double FermionClass::AcceptProb (int ptclToMove, State &newState)
{
  // First, reject if new state is already occupied -- Pauli exclusion.
  for (int ptcl=0; ptcl<NumParticles; ptcl++)
    if (ptcl != ptclToMove)
      if (equal(newState,OccupiedStates(ptcl)))
	return 0.0;

  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  return (exp (-beta * (newEnergy-oldEnergy)));
}


////////////////////////////////////////
// Boltzmannon Class Member Functions //
////////////////////////////////////////
double BoltzmannonClass::AcceptProb (int ptclToMove, State &newState)
{
  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  double acceptProb = exp (-beta * (newEnergy-oldEnergy));
  return acceptProb;
}



/////////////////////////////////
// BosonClass Member Functions //
/////////////////////////////////

double BosonClass::AcceptProb (int ptclToMove, State &newState)
{
  int Nold = 0;
  int Nnew = 0;
  for (int i=0; i<NumParticles; i++) {
    if (equal(OccupiedStates(i),newState))
      Nnew++;
    if (equal(OccupiedStates(i),OccupiedStates(ptclToMove)))
      Nold++;
  }

  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  double acceptProb = exp (-beta * (newEnergy-oldEnergy));

  // Boson statistics:
  acceptProb *= (double)(Nnew+1)/(double)Nold;
  return (acceptProb);
}


//////////
// Main //
//////////
main(int argc, char **argv)
{

  if (argc < 2) {
    cout << "Usage:\n";
    cout << "FreeParticles myfile.in\n"; 
  }
  else {
    COMM::Init(argc, argv);
    ParticleClass *particle;
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    string particleType;
    assert (in.ReadVar ("Type", particleType));
    if (particleType == "FERMION")
      particle = new FermionClass;
    else if (particleType == "BOSON")
      particle = new BosonClass;
    else if (particleType == "BOLTZMANNON")
      particle = new BoltzmannonClass;
    else {
      cerr << "Unrecognized particle type " << particleType << ". Exitting.\n";
      exit(1);
    }
      
    particle->Read (in);
    particle->Run();
    delete particle;
  }
}

