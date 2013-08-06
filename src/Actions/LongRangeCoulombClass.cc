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


#include "LongRangeCoulombClass.h"
#include "../PathDataClass.h"

// hack
//ofstream out;
bool isEnergy;
double halfbox;
int allPairs, excludedPairs;

LongRangeCoulombClass::LongRangeCoulombClass(PathDataClass &pathData,
			       Array<PairActionFitClass* ,2> &pairMatrix,
			       Array<PairActionFitClass*, 1> &pairArray) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix),
  PairArray(pairArray)
{
  elementary_charge = 1.602*pow(10.0,-19);
  N_Avogadro = 6.022*pow(10.0,23.0);
  kcal_to_joule = 4184;
  epsilon_not = 8.85*pow(10.0,-12);
  angstrom_to_m = pow(10.0,10);
  SI = 1/(4*M_PI*epsilon_not);
  prefactor = SI*angstrom_to_m*elementary_charge*elementary_charge*N_Avogadro/kcal_to_joule;
  initPhi = false;
  //out.open("LRCEnergies.dat");
  //out << "##kspace self real realTest" << endl;
  allPairs = 0;
  //excludedPairs = 0;
}

void LongRangeCoulombClass::Read(IOSectionClass& in)
{
  assert(in.ReadVar("Alpha",alpha));
  double conversion = 1.0;
  in.ReadVar("Conversion",conversion);
  prefactor = prefactor*conversion;
  in.ReadVar("Prefactor",prefactor);
  doRealSpace = true;
  in.ReadVar("CalcRealSpace",doRealSpace);
  cerr << "LRC: Overall coulomb prefactor " << prefactor << endl;
  sq_alpha = sqrt(alpha);
  Array<string, 1> specNames;
  assert(in.ReadVar("ActiveSpecies",specNames));
  activeSpecies.resize(Path.NumSpecies());
  activeSpecies = false;
  for(int s=0; s<specNames.size(); s++)
    activeSpecies(Path.SpeciesNum(specNames(s))) = true;
  cerr << "LongRangeCoulomb active species are " << specNames << " " << activeSpecies << endl;
  PtclCharge.resize(PathData.NumSpecies());
  for(int s=0; s<PathData.NumSpecies(); s++)
    PtclCharge(s) = PathData.Species(s).pseudoCharge;
  cerr << "Original Species charges:" << PtclCharge << endl;
  Array<double,1> setCharge;
  assert(in.ReadVar("ActiveCharges",setCharge));
  int chargeIndex = 0;
	for(int s=0; s<specNames.size(); s++) {
	  int myS = Path.SpeciesNum(specNames(s));
    PtclCharge(myS) = setCharge(chargeIndex);
    chargeIndex++;
  }
  cerr << "Assigned charges for each species:" << PtclCharge << endl;
  int kPts = Path.kVecs.size();
  cerr << "LongRangeCoulomb: using " << kPts << " kPts from Path.kVecs" << endl;
  phi.resize(kPts);
  k2.resize(kPts);
  for(int i=0; i<kPts; i++){
    double kSq = dot(Path.kVecs(i), Path.kVecs(i));
    k2(i) = kSq;
    phi(i) = exp(-kSq/(4*alpha*alpha));
  }
  cerr << "LRC Read loaded phi of size " << phi.size() << endl;

  self = 0.0;
  double a = alpha/sqrt(M_PI);
  for(int n=0; n<Path.NumParticles(); n++){
    int spec = Path.ParticleSpeciesNum(n);
    if(activeSpecies(spec)) {
      //double q = PathData.Species(spec).Charge;
      double q = PtclCharge(spec);
      self += prefactor*a*q*q;
    }
  }
  cerr << "Computed self-interaction correction " << self << endl;

  volume = 1.0;
  for(int i=0; i<3; i++)
    volume *= Path.GetBox()(i);
  cerr << "I have box volume " << volume << endl;
  halfbox = 0.5*Path.GetBox()(0);
  cerr << "initPhi " << initPhi << endl;
  initPhi = false;
  LRfirstTime = true;
}


inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}

double LongRangeCoulombClass::SingleAction (int slice1, int slice2,
			      const Array<int,1> &changedParticles,
			      int level)
{
  isEnergy = false;
  if (slice1 == 0 && slice2 == PathData.Path.TotalNumSlices) {
    slice2 -= 1;
  }
  double U = ComputeEnergy(slice1, slice2, changedParticles, level);
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  //cerr << "LRC returning U " << U << " * tau " << levelTau << " = " << levelTau*U << endl;
  return(levelTau*U);
}


double LongRangeCoulombClass::d_dBeta (int slice1, int slice2,  int level)
{
  isEnergy = true;
  Array<int,1> changedParticles(PathData.Path.NumParticles());
  for(int i=0; i<PathData.Path.NumParticles(); i++)
    changedParticles(i) = i;
  double U = ComputeEnergy(slice1, slice2-1, changedParticles, level);
  return(U);
}

double LongRangeCoulombClass::ComputeEnergy (int slice1, int slice2,
			      const Array<int,1> &changedParticles, int level)
{
  //cerr << "in LRCE " << firstTime << endl;
  if(LRfirstTime){
    int kPts = Path.kVecs.size();
    k2.resize(kPts);
    phi.resize(kPts);
    for(int i=0; i<kPts; i++){
      double kSq = dot(Path.kVecs(i), Path.kVecs(i));
      k2(i) = kSq;
      phi(i) = exp(-k2(i)/(4*alpha*alpha));
      //cerr << i << " " << kSq << " " << k2(i) << " " << phi(i) << endl;
    }
    cerr << "LRC initialized phi of size " << phi.size() << " inside action!" << endl;
    LRfirstTime= false;
    cerr << "activespec " << activeSpecies << endl;

    // initialize Rho_k
    for (int slice=0; slice<=Path.TotalNumSlices; slice+=1) {
      for (int species=0; species<Path.NumSpecies(); species++) {
        if(activeSpecies(species)) {
          Path.CalcRho_ks_Fast(slice, species);
        }
      }
    }
  }

  int skip = (1<<level);
  if (GetMode() == NEWMODE)
  {
    //Path.UpdateRho_ks(slice1, slice2, changedParticles, level);
    for (int slice=slice1; slice<=slice2; slice+=skip) {
      for (int species=0; species<Path.NumSpecies(); species++) {
        if(activeSpecies(species)) {
          Path.CalcRho_ks_Fast(slice, species);
        }
      }
    }
  }
  //cerr << "RHO CALC FAST " << endl;
  //for(int s=0; s<Path.TotalNumSlices; s++)
  //  for(int spec=0; spec<Path.NumSpecies(); spec++)
  //    for(int k=0; k<Path.kVecs.size(); k++)
  //      cerr << s << " " << spec << " " << k << " " << Path.Rho_k(s, spec, k) << endl;

  //for (int slice=slice1; slice<=slice2; slice+=skip) {
  //  for (int species=0; species<Path.NumSpecies(); species++) {
  //    Path.CalcRho_ks_Slow(slice, species);
  //  }
  //}
  //cerr << "RHO CALC SLOW" << endl;
  //for(int s=0; s<Path.TotalNumSlices; s++)
  //  for(int spec=0; spec<Path.NumSpecies(); spec++)
  //    for(int k=0; k<Path.kVecs.size(); k++)
  //      cerr << s << " " << spec << " " << k << " " << Path.Rho_k(s, spec, k) << endl;
  
  double kspace = 0.0;
  double real = 0.0;
  double realTest = 0.0;
  double correction = 0;
  //double vacuum = 0.0;
  // charge-neutral only!!
  double factor = 1.0/(volume);

  for (int slice=slice1; slice<=slice2; slice+=skip) {
    correction += self;


    // First, do the Real-space terms 
    if(doRealSpace) {
      for (int counter=0; counter<Path.DoPtcl.size(); counter++)
        Path.DoPtcl(counter)=true;
      for(int iIndex=0; iIndex<changedParticles.size(); iIndex++){
        int i = changedParticles(iIndex);
        Path.DoPtcl(i) = false;
        int speci = Path.ParticleSpeciesNum(i);
        if(activeSpecies(speci)) {
          //double qi = PathData.Species(speci).Charge;
          double qi = PtclCharge(speci);
          //cerr << "real: spec " << speci << " with q " << qi << endl;

          for(int j=0; j<PathData.Path.NumParticles(); j++){
            if(Path.DoPtcl(j)){
              int specj = Path.ParticleSpeciesNum(j);
              if(activeSpecies(specj)) {
                double qj = PtclCharge(specj);
                //double qj = PathData.Species(specj).Charge;
                if((qi*qj) != 0.0){
                  dVec Rij;
                  double rijmag;
	        	      PathData.Path.DistDisp(slice, i, j, rijmag, Rij);

                  if(PathData.Mol(i) != PathData.Mol(j)){
                    allPairs++;
                    dVec Oij;
                    double Oijmag;
	        	        PathData.Path.DistDisp(slice, PathData.Mol(i), PathData.Mol(j), Oijmag, Oij);
                    if(Oijmag < halfbox){
                      // real-space term
                      real += prefactor* qi*qj * erfc(alpha*rijmag)/rijmag;
                      //real += prefactor* qi*qj * (1 - erf(alpha*rijmag))/rijmag;
                      realTest += prefactor*qi*qj/rijmag;
                    }
                    //else{
                    //  excludedPairs++;
                    //  if(allPairs%500000 == 0){
                    //    cerr << "excluded " << excludedPairs << " pairs with ratio " << double(excludedPairs)/allPairs << endl;
                    //  }
                    //}
                  }
                }
              }
            }
          }
        }
        //dVec Ri = PathData.Path(slice,i);
        //vacuum += 2*M_PI/(3*volume) * qi*qi * dot(Ri,Ri);
      }
    }

    // Now k-space terms over all ptcl pairs!
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      //int speci = Path.ParticleSpeciesNum(species);
      if(activeSpecies(species)) {
        double qi = PtclCharge(species);
        //double qi = PathData.Species(species).Charge;
        //cerr << "homo: spec " << species << " with q " << qi << " Path.NumSpec is " << Path.NumSpecies() << endl;
        if(qi != 0.0){
          for (int ki=0; ki<Path.kVecs.size(); ki++) {
	          double rhok2 = mag2(Path.Rho_k(slice,species,ki));
            double myKterm = prefactor*factor * 2*M_PI * qi*qi * rhok2 * phi(ki)/k2(ki);
	          kspace += myKterm;
	          //kspace += prefactor*factor * 2*M_PI * qi*qi * rhok2 * phi(ki)/k2(ki);
            //cerr << "homo " << ki << " qi " << qi << " rhok2 " << rhok2 << " phi_k " << phi(ki) << " k^-2 " << 1.0/k2(ki) << " total " << myKterm << endl;
	        }
        }
      }
    }

    //cerr << "khomo " << kspace;
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++){
      if(activeSpecies(species1)) {
        for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
          if(activeSpecies(species2)) {
            //int speci = Path.ParticleSpeciesNum(species1);
            //int specj = Path.ParticleSpeciesNum(species2);
            //double qi = PathData.Species(species1).Charge;
            //double qj = PathData.Species(species2).Charge;
            double qi = PtclCharge(species1);
            double qj = PtclCharge(species2);
            if(qi != 0.0 && qj != 0.0){
	            for (int ki=0; ki<Path.kVecs.size(); ki++) {
	              double rhorho = 
	                Path.Rho_k(slice, species1, ki).real() *
	                Path.Rho_k(slice, species2, ki).real() + 
	                Path.Rho_k(slice, species1, ki).imag() *
	                Path.Rho_k(slice, species2, ki).imag();
	              double myKterm = prefactor*factor * 2*M_PI * qi*qj * rhorho * phi(ki)/k2(ki);
	              //kspace += prefactor*factor * 2*M_PI * qi*qj * rhorho * phi(ki)/k2(ki);
	              kspace += myKterm;
                //cerr << "hetero " << ki << " qi " << qi << " qj " << qj << " rhorho " << rhorho << " phi_k " << phi(ki) << " k^-2 " << 1.0/k2(ki) << " total " << myKterm << endl;
	            }
            }
          }
        }
      }
    }

    //cerr << " khetero " << kspace << endl;
    // assemble self-energy correction: need to include terms for other sites on the same molecule!!
    // if this is right, it can be done more efficiently (only need to update terms from molecules that were moved!
    double a = alpha/sqrt(M_PI);
    for(int i=0; i<PathData.Path.NumParticles(); i++){
      int speci = Path.ParticleSpeciesNum(i);
      if(activeSpecies(speci)) {
        double qi = PtclCharge(speci);
        //double qi = PathData.Species(speci).Charge;
        //correction += prefactor*a*qi*qi;
        //cerr << "ptcl " << i << " has friends";
        for(int jIndex=0; jIndex<PathData.Mol.MembersOf(PathData.Mol(i)).size(); jIndex++){
          int j = PathData.Mol.MembersOf(PathData.Mol(i))(jIndex);
          //cerr << " " << j;
          if(i!=j){
            int specj = Path.ParticleSpeciesNum(j);
            if(activeSpecies(specj)) {
              //double qj = PathData.Species(specj).Charge;
              double qj = PtclCharge(specj);
              dVec Rij;
              double rijmag;
	        	  PathData.Path.DistDisp(slice, i, j, rijmag, Rij);
              // not sure about factor of 0.5; don't think it's right
              correction += 0.5*prefactor* qi*qj * erf(alpha*rijmag)/rijmag;
            }
          }
        }
      }
    }
  }

  //if(isEnergy)
    //out << kspace << " " << correction << " " << real << " " << realTest << endl;

  //if(isEnergy)
  //  cerr << "LR ONLY KSPACE " << kspace << " " << self << " " << correction << " " << real << " " << realTest << endl;
  //cerr << "LR terms " << doRealSpace << " " << kspace << " " << self << " " << correction << " " << real << " " << realTest << endl;
  double U = (kspace + real - correction);
  //double U = kspace;
  return (U);
}

void LongRangeCoulombClass::GradAction(int slice1, int slice2, 
			   const Array<int,1> &ptcls, int level,
			   Array<dVec,1> &gradVec)
{
  
  int skip = 1<<level;
  
  for (int slice = slice1; slice<=slice2; slice+=skip) {
    double factor = ((slice==slice1)||(slice==slice2)) ? 1.0 : 2.0;
    for (int pi=0; pi<ptcls.size(); pi++) {
      int ptcl = ptcls(pi);
      int species1 = Path.ParticleSpeciesNum(ptcl);
      if(activeSpecies(species1)) {
        for (int ki=0; ki<Path.kVecs.size(); ki++) {
          dVec r = Path(slice,ptcl);
          dVec k = Path.kVecs(ki);
          double phi = dot(r,k);
          double re, im;
	  //          sincos(phi, &im, &re);
	  re=cos(phi);
	  im=sin(phi);
          complex<double> z(re, im);
          complex<double> rho_uSum(0.0, 0.0);
          for (int si=0; si<Path.NumSpecies(); si++) {
            if(activeSpecies(si)) {
              PairActionFitClass &PA = *PairMatrix(species1,si);
              rho_uSum += PA.Ulong_k(level,ki) * Path.Rho_k(slice, si, ki);
            }
          }
          rho_uSum = conj (rho_uSum);
          // Now, compute the imaginary part of the product, and
          // multiply by k
          gradVec(pi) -= 
            factor*(z.real()*rho_uSum.imag() + z.imag()*rho_uSum.real())*k;
        }
      }
    } // end pi loop
  } // end slice loop

}

string
LongRangeCoulombClass::GetName()
{
  return "LongRange";
}
