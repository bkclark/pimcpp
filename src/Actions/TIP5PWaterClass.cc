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

#include "../PathDataClass.h"
#include "TIP5PWaterClass.h"

TIP5PWaterClass::TIP5PWaterClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  //Do  nothing for now
}

bool TIP5PWaterClass::dVecsEqual(dVec u,dVec v){
  bool equal = true;
  for (int i = 0; i < 3; i++){
    if (u(i) != v(i))
      equal = false;
  }
  return equal;
}


double 
TIP5PWaterClass::SingleAction (int startSlice, int endSlice, 
			       const Array<int,1> &activeParticles, int level)
{
//  cerr << "I'm calculating the Action" << endl;
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }

// obnoxious physical constants
  double elementary_charge = 1.602*pow(10.0,-19);
  double N_Avogadro = 6.022*pow(10.0,23.0);
  double kcal_to_joule = 4184;
  double epsilon_not = 8.85*pow(10.0,-12);
  double angstrom_to_m = pow(10.0,10);
  double SI = 1/(4*M_PI*epsilon_not);
  double k_B = 1.3807*pow(10.0,-23);
  double erg_to_eV = 1.6*pow(10.0,12);
  double joule_to_eV = pow(10.0,19)/1.602;

  double TotalU = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  int speciesO=Path.SpeciesNum("O");
  int speciesp=Path.SpeciesNum("p");
  int speciese=Path.SpeciesNum("e");
  double CUTOFF = 7.75; // setcutoff: spherical cutoff in angstroms

  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
/*   need to exclude intramolecular interactions - compare with line commented out below
      for(int i = 0;i<=4;++i){
        Path.DoPtcl(ptcl1+4*i) = false;
      }
*/
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
//    cerr << "Choosing particle " << ptcl1 << " which is of species " << species1 << endl;

    if (species1==speciesO){
      int counter = 0;
//    Calculate o-o interaction
      for (int ptcl2=Path.Species(speciesO).FirstPtcl;ptcl2<=Path.Species(speciesO).LastPtcl;ptcl2++) {///loop over oxygen
	if (Path.DoPtcl(ptcl2)){
//          int PairIndex = PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2));
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    dVec r, rp;
	    double rmag, rpmag;
	    PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
//  implement spherical cutoff  
            if (rmag <= CUTOFF){
              counter ++;
              double sigma_over_r = PathData.Species(species1).Sigma/rmag;
              double sigma_over_cutoff = PathData.Species(species1).Sigma/CUTOFF;
              double offset = pow(sigma_over_cutoff,12) - pow(sigma_over_cutoff,6);
	      double lj = 4*PathData.Species(species1).Epsilon*(pow(sigma_over_r,12)-pow(sigma_over_r,6) - offset); // this is in kcal/mol 
	      TotalU += lj;
//	      cerr << "lj  " << lj << " at distance " << rmag << endl;
            }
            else{
//              cerr << rmag << " is outside cutoff" << endl;
            }
	  }
	}
      }
      ///end calculating o-o interactions
//      cerr << "I calculated Lennard-Jones interactions with: " << counter << " particles" << endl;
    }
    else{
      /// calculating coulomb interactions
      for (int ptcl2=Path.Species(speciesp).FirstPtcl;ptcl2<=Path.Species(speciesp).LastPtcl;ptcl2++) {
        ///loop over protons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciesp).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const*(1.0/rmag - 1.0/CUTOFF);
	      TotalU += coulomb;
//	      cerr << "protons " << coulomb << " at distance " << rmag  << " with offset " << coulomb_const/CUTOFF << endl;
            }
          }
	}
      }
      for (int ptcl2=Path.Species(speciese).FirstPtcl;ptcl2<=Path.Species(speciese).LastPtcl;ptcl2++) {
        ///loop over electrons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciese).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const*(1.0/rmag - 1.0/CUTOFF);
	      TotalU += coulomb;
//	      cerr << "electrons " << coulomb << " at distance " << rmag << " with offset " << coulomb_const/CUTOFF << endl;
            }
	  }
	}
      }
    
	
    }
   /// end calculating coulomb interactions 
  
  
  }
  double TotalU_times_tau = TotalU*PathData.Path.tau;
//  cerr << TotalU << " and times tau " << TotalU_times_tau << " at temp " << 1.0/PathData.Path.tau << endl;
//  cerr << "I'm returning TIP5P action " << TotalU_times_tau << endl;
  return (TotalU_times_tau);
}

double TIP5PWaterClass::d_dBeta (int startSlice, int endSlice,  int level)
{
  double thermal = 6/PathData.Path.tau;
  Array<int,1> activeParticles(PathData.Path.NumParticles());
  for (int i=0;i<PathData.Path.NumParticles();i++){
    activeParticles(i)=i;
  }
 
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }

// obnoxious physical constants
  double elementary_charge = 1.602*pow(10.0,-19);
  double N_Avogadro = 6.022*pow(10.0,23.0);
  double kcal_to_joule = 4184;
  double epsilon_not = 8.85*pow(10.0,-12);
  double angstrom_to_m = pow(10.0,10);
  double SI = 1/(4*M_PI*epsilon_not);
  double k_B = 1.3807*pow(10.0,-23);
  double erg_to_eV = 1.6*pow(10.0,12);
  double joule_to_eV = pow(10.0,19)/1.602;

  double TotalU = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  int speciesO=Path.SpeciesNum("O");
  int speciesp=Path.SpeciesNum("p");
  int speciese=Path.SpeciesNum("e");
  double CUTOFF = 7.75; // setcutoff: spherical cutoff in angstroms

  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
/*   need to exclude intramolecular interactions - compare with line commented out below
      for(int i = 0;i<=4;++i){
        Path.DoPtcl(ptcl1+4*i) = false;
      }
*/
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
//    cerr << "Choosing particle " << ptcl1 << " which is of species " << species1 << endl;

    if (species1==speciesO){
      int counter = 0;
//    Calculate o-o interaction
      for (int ptcl2=Path.Species(speciesO).FirstPtcl;ptcl2<=Path.Species(speciesO).LastPtcl;ptcl2++) {///loop over oxygen
	if (Path.DoPtcl(ptcl2)){
//          int PairIndex = PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2));
	  
	  for (int slice=startSlice;slice<endSlice;slice+=skip){
//cerr << "slice " << slice << endl;
	    dVec r, rp;
	    double rmag, rpmag;
	    PathData.Path.DistDisp(slice, ptcl1, ptcl2,
				   rmag, r);
//  implement spherical cutoff  
            if (rmag <= CUTOFF){
              counter ++;
              double sigma_over_r = PathData.Species(species1).Sigma/rmag;
	      double lj = 4*PathData.Species(species1).Epsilon*(pow(sigma_over_r,12)-pow(sigma_over_r,6)); // this is in kcal/mol 
	      TotalU += lj;
//	      cerr << "lj  " << lj << " at distance " << rmag << endl;
            }
	  }
	}
      }
      ///end calculating o-o interactions
//      cerr << "I calculated Lennard-Jones interactions with: " << counter << " particles" << endl;
    }
    else{
      /// calculating coulomb interactions
      for (int ptcl2=Path.Species(speciesp).FirstPtcl;ptcl2<=Path.Species(speciesp).LastPtcl;ptcl2++) {///loop over protons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  for (int slice=startSlice;slice<endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciesp).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const/rmag;
	      TotalU += coulomb;
//	      cerr << "protons " << coulomb << " at distance " << rmag  << " with offset " << coulomb_const/CUTOFF << endl;
            }
          }
	}
      }
      for (int ptcl2=Path.Species(speciese).FirstPtcl;ptcl2<=Path.Species(speciese).LastPtcl;ptcl2++) {///loop over electrons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  for (int slice=startSlice;slice<endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciese).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const/rmag;
	      TotalU += coulomb;
//	      cerr << "electrons " << coulomb << " at distance " << rmag << " with offset " << coulomb_const/CUTOFF << endl;
            }
	  }
	}
      }
    
	
    }
   /// end calculating coulomb interactions 
  
  
  }
//  double TotalU_times_beta = TotalU*kcal_to_joule/(N_Avogadro*k_B)*PathData.Path.tau;
//  cerr << TotalU << " and times beta " << TotalU_times_beta << " at temp " << 1.0/PathData.Path.tau << endl;
//  cerr << "Energy function is returning " << TotalU << endl;

  double energy_per_molecule = TotalU/PathData.Mol.NumMol();
//  return energy_per_molecule; // + thermal;
  return TotalU;
}

double TIP5PWaterClass::OOSeparation (int slice,int ptcl1,int ptcl2)
{
  int speciesO=Path.SpeciesNum("O");
  int Optcl1, Optcl2;
  dVec Or;
  double Ormag;

  Optcl1 = Path.Species(speciesO).FirstPtcl + PathData.Mol(ptcl1);
  Optcl2 = Path.Species(speciesO).FirstPtcl + PathData.Mol(ptcl2);
  PathData.Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Or);
//  cerr << "We have particles " << ptcl1 << " and " << ptcl2 << " belonging to molecules " << PathData.Mol(ptcl1) << " and " << PathData.Mol(ptcl2) << " and calculate the distance between O " << Optcl1 << " and " << Optcl2 << " from " << PathData.Mol(Optcl1) << " and " << PathData.Mol(Optcl2) << endl;
  return Ormag;
}

void TIP5PWaterClass::Read (IOSectionClass &in)
{
  //do nothing for now
}

double TIP5PWaterClass::RotationalKinetic(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level)
{
  double RotK = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = activeParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double FourLambdaTauInv=1.0/(4.0*lambda_p*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      // Load the coordinates we'll be manipulating; these are coords for adjacent time slices of a proton.
      dVec coord1 = PathData.Path(slice,ptcl);
      dVec coord2 = PathData.Path(slice+skip,ptcl);
//cerr << "loaded " << coord1 << " and " << coord2 << " for slice " << slice << " and ptcl " << ptcl << endl;
      // Identify the index for the corresponding oxygen, i.e. COM
      int Optcl = FindCOM(ptcl);
      dVec Ocoord1 = PathData.Path(slice,Optcl);
      dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << Ocoord1 << " and " << Ocoord2 << endl;
      // We want to measure rotations relative to the COM ONLY, so we need to redefine coords WRT their respective COMs so that we can calculate polar and azimuthal angles WRT the COM.
      // Sync COMs for coord1 and coord2
//      dVec COMdisp = Displacement(slice,slice+skip,Optcl,Optcl);
//      coord2 -= COMdisp;
      // Redefine coordinates WRT COM
      coord1 -= Ocoord1;
      coord2 -= Ocoord2;
//cerr << "now WRT COM coords are " << coord1 << " and " << coord2 << endl;
      // Get the corresponding polar and azimuthal angles
      double theta1,theta2,phi1,phi2;
      GetAngles(coord1,theta1,phi1);
      GetAngles(coord2,theta2,phi2); 
//cerr << "corresponding angles are " << theta1 << ", " << phi1 << ", " << theta2 << ", " << phi2 << endl;
      // Calculate the rotational kinetic energy
      double vel_squared = CalcEnergy(theta1,theta2-theta1,phi2-phi1);
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

      double GaussProd = 1.0;
//    for (int dim=0; dim<NDIM; dim++) {
//  	int NumImage=1;
      double GaussSum=0.0;
//	for (int image=-NumImage; image<=NumImage; image++) {
//	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
      GaussSum += exp(-vel_squared*FourLambdaTauInv);
//      }

      GaussProd *= GaussSum;
//      }
      RotK -= log(GaussProd);    
    //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//  //cerr << "I'm returning kinetic action " << RotK << endl;
  return (RotK);
}

double TIP5PWaterClass::RotationalEnergy(int startSlice, int endSlice, int level)
{
  double spring=0.0;
  // ldexp(double x, int n) = x*2^n
  double levelTau=ldexp(Path.tau, level);
//   for (int i=0; i<level; i++) 
//     levelTau *= 2.0;
  spring  = 0.0;  
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double FourLambdaTauInv = 1.0/(4.0*lambda_p*levelTau);
  int TotalNumParticles = Path.NumParticles();
  int startparticle = 3*TotalNumParticles/5;
  int endparticle = 4*TotalNumParticles/5;
  for (int ptcl=startparticle; ptcl<endparticle; ptcl++) {
    // Do free-particle part
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    if (speciesNum == PathData.Path.SpeciesNum("p")){
      SpeciesClass &species = Path.Species(speciesNum);
      for (int slice=startSlice; slice<endSlice; slice+=skip) {
	spring += (0.5*2)/levelTau;

      // Load the coordinates we'll be manipulating; these are coords for adjacent time slices of a proton.
      dVec coord1 = PathData.Path(slice,ptcl);
      dVec coord2 = PathData.Path(slice+skip,ptcl);
//cerr << "loaded " << coord1 << " and " << coord2 << " for slice " << slice << " and ptcl " << ptcl << endl;
      // Identify the index for the corresponding oxygen, i.e. COM
      int Optcl = FindCOM(ptcl);
      dVec Ocoord1 = PathData.Path(slice,Optcl);
      dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << Ocoord1 << " and " << Ocoord2 << endl;
      // We want to measure rotations relative to the COM ONLY, so we need to redefine coords WRT their respective COMs so that we can calculate polar and azimuthal angles WRT the COM.
      // Sync COMs for coord1 and coord2
//      dVec COMdisp = Displacement(slice,slice+skip,Optcl,Optcl);
//      coord2 -= COMdisp;
      // Redefine coordinates WRT COM
      coord1 -= Ocoord1;
      coord2 -= Ocoord2;
//cerr << "now WRT COM coords are " << coord1 << " and " << coord2 << endl;
      // Get the corresponding polar and azimuthal angles
      double theta1,theta2,phi1,phi2;
      GetAngles(coord1,theta1,phi1);
      GetAngles(coord2,theta2,phi2); 
//cerr << "corresponding angles are " << theta1 << ", " << phi1 << ", " << theta2 << ", " << phi2 << endl;
      // Calculate the rotational kinetic energy
      double vel_squared = CalcEnergy(theta1,theta2-theta1,phi2-phi1);
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

        double GaussSum;
        double numSum;
//      dVec GaussSum=0.0;
//      dVec numSum=0.0;
//      for (int dim=0; dim<NDIM; dim++) {
//        for (int image=-NumImage; image<=NumImage; image++) {
//	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	double d2overFLT = vel_squared*FourLambdaTauInv;
	double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = -d2overFLT/levelTau* expPart; 
//	GaussSum[dim] += expPart;
//	numSum[dim] += -d2overFLT/levelTau* expPart;
//	}
//	Z *= GaussSum[dim];
        Z += GaussSum;
//      }
/*      double scalarnumSum=0.0;
        for (int dim=0;dim<NDIM;dim++) {
          dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
        }
*/
	//cerr << "Z = " << Z << endl;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  cerr << "spring = " << spring << endl;
  return spring;
}

void TIP5PWaterClass::GetAngles(dVec disp, double &theta, double &phi)
{
  double x = disp(0);
  double y = disp(1);
  double z = disp(2);
  double R = O_H_moment_arm;
  theta = acos(z/R);
  double SineTheta = sin(theta);  
  phi = acos(x/(R*SineTheta));
  double phicheck = asin(y/(R*SineTheta));
  if(y<0){
    phi = 2*M_PI - phi;
  }
}

double TIP5PWaterClass::CalcEnergy(double reftheta,double dtheta, double dphi)
{
  double R = O_H_moment_arm;
  double SineTheta = sin(reftheta);
  double omega_squared = dtheta*dtheta;// + dphi*dphi*SineTheta*SineTheta;
  double vel_squared = omega_squared*R*R;
  return vel_squared;
}

double TIP5PWaterClass::SineOfPolar(dVec coords)
{
  double x = coords(0);
  double z = coords(2);
  double R = O_H_moment_arm;
  double theta = acos(z/R);  
cerr << "I find theta = " << theta << " for coords " << x << ", " << coords(1) << ", " << z << ";";
  double SineTheta = sin(theta);
cerr << "returning " << SineTheta << endl;
  return SineTheta;
}

dVec TIP5PWaterClass::COMVelocity (int slice1,int slice2,int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  dVec Ovel;

  Optcl = Path.Species(speciesO).FirstPtcl + PathData.Mol(ptcl);
  Ovel = PathData.Path.Velocity(slice1, slice2, Optcl);
//cerr << "I'm correcting velocity for ptcl " << ptcl << " of species " << Path.ParticleSpeciesNum(ptcl) << ".  Found COM oxygen at " << Optcl;
//cerr << "Returning COM velocity " << Ovel << endl;
  return Ovel;
}

dVec TIP5PWaterClass::COMCoords (int slice, int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  dVec relative_coords;

  Optcl = Path.Species(speciesO).FirstPtcl + PathData.Mol(ptcl);
  relative_coords = PathData.Path(slice,ptcl) - PathData.Path(slice,Optcl);
//cerr << "I'm correcting velocity for ptcl " << ptcl << " of species " << Path.ParticleSpeciesNum(ptcl) << ".  Found COM oxygen at " << Optcl;
//cerr << "Returning COM velocity " << Ovel << endl;
  return relative_coords;
}

dVec TIP5PWaterClass::Displacement(int slice1, int slice2, int ptcl1, int ptcl2)
{
  dVec disp;
  disp = PathData.Path(slice1,ptcl1) - PathData.Path(slice2,ptcl2);
  return disp;
}

int TIP5PWaterClass::FindCOM(int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  Optcl = Path.Species(speciesO).FirstPtcl + PathData.Mol(ptcl);
  return Optcl;
}

int TIP5PWaterClass::FindOtherProton(int ptcl)
{
  int speciesp=PathData.Path.SpeciesNum("p");
  int otherptcl;
  otherptcl = Path.Species(speciesp).FirstPtcl + PathData.Mol(ptcl);
  if (otherptcl == ptcl){
    otherptcl += Path.NumParticles()/5;
  }
  return otherptcl;
}

double TIP5PWaterClass::dotprod(dVec vec1, dVec vec2, double mag)
{
  double total = 0;
  for(int i = 0; i<3; i++){
    total += vec1[i]*vec2[i];
  }
  double norm = 1.0/mag;
  return total*norm;
}

/*
double TIP5PWaterClass::ProtonKineticAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  double R = O_H_moment_arm;
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = lambda_p;
//    double lambda = Path.Species(species).lambda;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
// now calculate the dotprod between these UNIT vectors
        double dot = dotprod(coord1,coord2,R*R);
// theta = arccos[dotprod of unit vectors]
        double theta = acos(dot);
        double vel_squared = R*R*theta*theta;
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

        double GaussProd = 1.0;
        double GaussSum=0.0;
        GaussSum += exp(-vel_squared*FourLambdaTauInv);
        GaussProd *= GaussSum;
        TotalK -= log(GaussProd);    
      }
    }
  }
  //cerr << "I'm returning kinetic action " << RotK << endl;
  return (TotalK);
}

double TIP5PWaterClass::ProtonKineticEnergy (int slice1, int slice2, int level)
{
  double R = O_H_moment_arm;
  double spring=0.0;
  double Z = 0.0;
  double levelTau=ldexp(Path.tau, level);
  spring  = 0.0;  
  int skip = 1<<level;
  const int NumImage=1;  
  int TotalNumParticles = Path.NumParticles();
  int startparticle = 3*TotalNumParticles/5;
  int endparticle = 4*TotalNumParticles/5;
  for (int ptcl=startparticle; ptcl<endparticle; ptcl++) {
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = lambda_p;
//    double lambda = species.lambda;
    if (speciesNum == PathData.Path.SpeciesNum("p")) {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
	spring += (0.5*2)/levelTau;
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
// now calculate the dotprod between these UNIT vectors
        double dot = dotprod(coord1,coord2,R*R);
// theta = arccos[dotprod of unit vectors]
        double theta = acos(dot);
        double vel_squared = R*R*theta*theta;
//cerr << "from which I calculate vel_squared " << vel_squared << endl;
      
        double GaussSum;
        double numSum;
	double d2overFLT = vel_squared*FourLambdaTauInv;
	double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = -d2overFLT/levelTau* expPart; 
        Z += GaussSum;
	//cerr << "Z = " << Z << endl;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  cerr << "spring = " << spring << endl;
  return spring;
}

*/
//  Original versions of these functions

double TIP5PWaterClass::ProtonKineticAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = lambda_p;
//    double lambda = Path.Species(species).lambda;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
        dVec vel;
        vel = coord1 - coord2;
        double GaussProd = 1.0;
        for (int dim=0; dim<NDIM; dim++) {
  	  int NumImage=1;
	  double GaussSum=0.0;
	  for (int image=-NumImage; image<=NumImage; image++) {
	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	    GaussSum += exp(-dist*dist*FourLambdaTauInv);
	  }
	  GaussProd *= GaussSum;
        }
        TotalK -= log(GaussProd);    
      }
    }
  }
  return (TotalK);
}

double TIP5PWaterClass::ProtonKineticEnergy (int slice1, int slice2, int level)
{
  double spring=0.0;
  double levelTau=ldexp(Path.tau, level);
  spring  = 0.0;  
  int skip = 1<<level;
  const int NumImage=1;  
  int TotalNumParticles = Path.NumParticles();
  int startparticle = 3*TotalNumParticles/5;
  int endparticle = 4*TotalNumParticles/5;
  for (int ptcl=startparticle; ptcl<endparticle; ptcl++) {
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = lambda_p;
//    double lambda = species.lambda;
    if ((speciesNum == PathData.Path.SpeciesNum("p")))
//    if ((speciesNum == PathData.Path.SpeciesNum("p")) || (speciesNum == PathData.Path.SpeciesNum("e")))
    {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
//	spring += (0.5*NDIM)/levelTau;
	spring += (0.5*2)/levelTau;
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
        dVec vel;
        vel = coord1 - coord2;
	double Z = 1.0;
	dVec GaussSum=0.0;
	dVec numSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImage; image<=NumImage; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    double d2overFLT = dist*dist*FourLambdaTauInv;
	    double expPart = exp(-d2overFLT);
	    GaussSum[dim] += expPart;
	    numSum[dim] += -d2overFLT/levelTau* expPart;
	  }
	  Z *= GaussSum[dim];
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++) {
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
	} //cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
	spring += scalarnumSum/Z; 
      }
    }
  }
  //  cerr << "spring = " << spring << endl;
  return spring;
}

double TIP5PWaterClass::SecondProtonKineticAction(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level)
{
  double RotK = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = activeParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double FourLambdaTauInv=1.0/(4.0*lambda_p*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      dVec X1 = PathData.Path(slice,ptcl);
      dVec X2 = PathData.Path(slice+skip,ptcl);
//cerr << "loaded " << coord1 << " and " << coord2 << " for slice " << slice << " and ptcl " << ptcl << endl;
      int Optcl = FindCOM(ptcl);
      int otherptcl = FindOtherProton(ptcl);
//cerr << "I found other proton " << otherptcl << " and oxygen " << Optcl << " for ptcl " << ptcl << endl;
      dVec O1 = PathData.Path(slice,Optcl);
      dVec O2 = PathData.Path(slice+skip,Optcl);
      dVec P1 = PathData.Path(slice,otherptcl);
//cerr << "identified COM " << Optcl << " with coords " << Ocoord1 << " and " << Ocoord2 << endl;
      X1 -= O1;
      X2 -= O2;
      P1 -= O1;
//cerr << "now WRT COM coords are " << coord1 << " and " << coord2 << endl;
      dVec L = O1 - P1;
      dVec p1 = X1 - P1;
      dVec p2 = X2 - P1;
      double Beta = M_PI - HOH_angle;
      double CosBeta = cos(Beta);
      double SinBeta = sin(Beta);
      double scale = 1 + CosBeta;
      L = L*scale;
      dVec u1 = p1 - L;
      dVec u2 = p2 - L;
      double norm = dotprod(u1,u1,1.0);
      norm *= dotprod(u2,u2,1.0);
      double dot = dotprod(u1,u2,norm);
      double psi = acos(dot);
      double lever_arm = O_H_moment_arm*SinBeta;      
      // Calculate the rotational kinetic energy
      double vel_squared = psi*psi*lever_arm*lever_arm;
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

      double GaussProd = 1.0;
//    for (int dim=0; dim<NDIM; dim++) {
//  	int NumImage=1;
      double GaussSum=0.0;
//	for (int image=-NumImage; image<=NumImage; image++) {
//	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
      GaussSum += exp(-vel_squared*FourLambdaTauInv);
//      }

      GaussProd *= GaussSum;
//      }
      RotK -= log(GaussProd);    
    //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//  //cerr << "I'm returning kinetic action " << RotK << endl;
  return (RotK);
}

double TIP5PWaterClass::SecondProtonKineticEnergy(int startSlice, int endSlice, int level)
{
  double spring=0.0;
  double levelTau=ldexp(Path.tau, level);
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double FourLambdaTauInv = 1.0/(4.0*lambda_p*levelTau);
  int TotalNumParticles = Path.NumParticles();
  int startparticle = 4*TotalNumParticles/5;
  int endparticle = 5*TotalNumParticles/5;
  for (int ptcl=startparticle; ptcl<endparticle; ptcl++) {
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    if (speciesNum == PathData.Path.SpeciesNum("p")){
//      SpeciesClass &species = Path.Species(speciesNum);
      for (int slice=startSlice; slice<endSlice; slice+=skip) {
	spring += (0.5*1)/levelTau;
        dVec X1 = PathData.Path(slice,ptcl);
        dVec X2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        int otherptcl = FindOtherProton(ptcl);
//cerr << "p2 coords are " << X1 << " and " << X2 << ".  For ptcl " << ptcl << " found other proton " << otherptcl << " and oxygen " << Optcl << endl;
        dVec O1 = PathData.Path(slice,Optcl);
        dVec O2 = PathData.Path(slice+skip,Optcl);
        dVec P1 = PathData.Path(slice,otherptcl);
        dVec P2 = PathData.Path(slice,otherptcl);
//cerr << "Coords for oxygens and other proton at slice 1 are ";
//cerr << O1 << endl;
//cerr << O2 << endl;
//cerr << P1 << endl;

//        X1 -= O1;
//        X2 -= O2;
//        P1 -= O1;
//cerr << "NOT correcting for COM we have " << endl;
//cerr << "X1 " << X1 << endl;
//cerr << "X2 " << X2 << endl;
//cerr << "P1 " << P1 << endl;
        dVec L = O1 - P1;
        dVec p1 = X1 - P1;
        dVec p2 = X2 - P2;
//cerr << "L " << L << endl;
//cerr << "p1 " << p1 << endl;
//cerr << "p2 " << p2 << endl;
        double Beta = M_PI - HOH_angle;
//cerr << "beta is " << Beta << endl;
        double CosBeta = cos(Beta);
        double SinBeta = sin(Beta);
        double scale = 1 + CosBeta;
        L = L*scale;
//cerr << "so with scale " << scale << " L becomes " << L << endl;
        dVec u1 = p1 - L;
        dVec u2 = p2 - L;

//cerr << "u1 " << u1 << endl;
//cerr << "u2 " << u2 << endl;
        double norm = sqrt(dotprod(u1,u1,1.0));
//cerr << "normalization is " << norm;
        norm *= sqrt(dotprod(u2,u2,1.0));
//cerr << " and now it's " << norm << endl;
        double dot = dotprod(u1,u2,norm);
        double psi = acos(dot);
//cerr << "so i get the angle psi " << psi << endl;
        double lever_arm = O_H_moment_arm*SinBeta;      
        double vel_squared = psi*psi*lever_arm*lever_arm;
//cerr << "and vel_squared is " << vel_squared << endl;
        double GaussSum;
        double numSum;
        double d2overFLT = vel_squared*FourLambdaTauInv;
        double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = -d2overFLT/levelTau* expPart; 
        Z += GaussSum;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  //cerr << "spring = " << spring << endl;
  return spring;
}

double TIP5PWaterClass::FixedAxisAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level)
{
  double R = O_H_moment_arm;
//cerr << "ok here we go: R is " << R << endl;
  double RotK = 0.0;
//cerr << "active is " << activeParticles << endl;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double dt = levelTau*hbar;
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl1 = activeParticles(ptclIndex);
    int ptcl2 = ptcl1 + numMol;
    double FourLambdaTauInv=1.0/(4.0*lambda_p*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
//cerr << "ptcl1 is " << ptcl1 << " and ptcl2 is " << ptcl2 << " at slice " << slice << endl;
      // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << O << " and " << Oprime << endl;
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
//cerr << "now WRT COM coords are " << endl;
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
//cerr << "n " << n << endl;
//cerr << "nprime " << nprime << endl;
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
//cerr << "Well I calculate that the angle between n and nprime is " << acos(dotprod(n,nprime,1.0)) << endl;
     // if (n == nprime){
      if (dVecsEqual(n,nprime)){
        r = z1;
        theta = 0.0;
//cerr << "axis of rotation is z (theta = 0)" << r << endl;
      }
      else{
        r = Normalize(CrossProd(n,nprime));
        theta = GetAngle(n,nprime);
//cerr << "axis of rotation is r" << r << endl;
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
      double alpha = HOH_half_angle;
      double SinAlpha = sin(alpha);
      double CosAlpha = cos(alpha);
      double phi = GetAngle(z1,r);
      double psi = GetAngle(z1prime,r);
//cerr << "P1 " << P1 << endl;
//cerr << "n " << n << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "phi " << phi << endl;
//cerr << "psi " << psi << endl;
/*      if (phi = 3.14159)
        phi = 0.0;
      if (psi = 3.14159)
        psi = 0.0;
*/
//cerr << "total PHI " << GetAngle(z1,z1prime) << " vs. sum phi+psi " << phi+psi << endl;
// This getup selects the smaller angle, "mechanically" allowing proton exchange, although it doesn't do anything about parity i.e. sign flips so it's kind of wrong.
/*    if(phi1<phi2){
        phi = phi1;
      }
      else{
        phi = phi2;
      }
//cerr << "I chose psi " << psi << endl;
      if(psi1<psi2){
        psi = psi1;
      }
      else{
        psi = psi2;
      }
*/


      // Calculate angle of rotation (psi)

      if (phi == 0 && psi == 0 && theta == 0){
//cerr << "___________________________________________________zeroing___________" << endl;
        vel_squared = 0;
        prefactor = 1;
      }
      else{
//cerr << "FIXED AXIS TIP" << endl;
        // Calculate quaternions
        double wx = theta/dt;
        double wy = phi/dt;
        double wz = psi/dt;
        double w2 = wx*wx + wy*wy + wz*wz;
/* THESE PARAMTERS ASSUME CONSTANT ANGULAR VELOCITY
        double chi = 1 - w2*dt*dt/8; // not sure about this one: omega dot omega_not??
        double xi = -0.5*wy*dt + 0.25*(0.5*w2*wy)*(dt*dt*dt/6);
        double eta = 0.5*wx*dt + 0.25*(-0.5*w2*wx)*(dt*dt*dt/6);
        double zeta = 0.5*wz*dt + 0.25*(-0.5*w2*wz)*(dt*dt*dt/6);
// THESE ARE THE PARAMETERS AS PRINTED
 
       double chi = 1 - w2*dt*dt/24; // not sure about this one: omega dot omega_not??
        double xi = -0.75*wy*dt + 0.25*(0.5*w2*wy - 2*wy/(dt*dt))*(dt*dt*dt/6);
        double eta = 0.75*wx*dt + 0.25*(-0.5*w2*wx + 2*wx/(dt*dt))*(dt*dt*dt/6);
        double zeta = 0.25*wz*dt + 0.25*(wx*wz/dt - wx*wy/dt - 0.5*w2*wz + 2*wz/(dt*dt))*(dt*dt*dt/6);

        // Calculate Fixed-Axis Approximation parameters
        double gamma = 2*acos(chi);
        double SinHalfGamma = sin(gamma/2);
//if (chi < -1){
cerr << "EXPANSION GIVES " << endl;
  cerr << "theta " << theta << endl;
  cerr << "phi " << phi << endl;
  cerr << "psi " << psi << endl;
  cerr << "dt " << dt << " wx " << wx << " wy " << wy << " wz " << wz << " w2 " << w2*dt*dt << endl;
  cerr << "w2*dt*dt/4 " << w2*dt*dt/4 << endl;
  cerr << "chi " << chi << " eta " << eta << " xi " << xi << " zeta " << zeta << endl;
  cerr << "gamma " << gamma << " and SinHalfGamma " << SinHalfGamma << endl;
//}
*/
// THIS IS TESTING WHETHER OR NOT WE CAN USE EXACT EXPRESSIONS ITO EULER ANGLES?
        double theta1 = 0;
	double theta2 = theta;
	double phi1= -phi;
	double phi2 = 0;
	double psi1 = 0;
	double psi2 = psi;
        double chi1 = cos(theta1/2)*cos((psi1 + phi1)/2);
        double eta1 = sin(theta1/2)*cos((psi1 - phi1)/2);
        double xi1 = sin(theta1/2)*sin((psi1 - phi1)/2);
        double zeta1 = cos(theta1/2)*sin((psi1 + phi1)/2);
        double chi2 = cos(theta2/2)*cos((psi2 + phi2)/2);
        double eta2 = sin(theta2/2)*cos((psi2 - phi2)/2);
        double xi2 = sin(theta2/2)*sin((psi2 - phi2)/2);
        double zeta2 = cos(theta2/2)*sin((psi2 + phi2)/2);
	double chi = chi2 - chi1;
	double xi = xi2-xi1;
	double zeta = zeta2 - zeta1;
	double eta = eta2 - eta1;
        double gamma = 2*acos(chi);
        double SinHalfGamma = sin(gamma/2);
//cerr << "EXPLICIT FORM GIVES " << endl;
//  cerr << "chi " << chi << " eta " << eta << " xi " << xi << " zeta " << zeta << endl;
//  cerr << "gamma " << gamma << " and SinHalfGamma " << SinHalfGamma << endl;

if (chi > 1 || chi < -1){
  cerr << "OH CRAP CHI IS " << chi << endl;
}
        // Rotation Moment Matrix
        double Lx = 2*R*R;
        double Ly = 2*R*R*SinAlpha*SinAlpha;
        double Lz = 2*R*R*SinAlpha*SinAlpha;
        // Axis of Rotation;
        double ex = eta/SinHalfGamma;
        double ey = -xi/SinHalfGamma;
        double ez = zeta/SinHalfGamma;
        // Rotation Matrix Elements <e|L|e>
        double Mx = ex*Lx*ex;
        double My = ey*Ly*ey;
        double Mz = ez*Lz*ez;
        double MSum = Mx + My + Mz;
        vel_squared = MSum*gamma*gamma;
        prefactor = gamma/(2*SinHalfGamma);
      }
//cerr << "from which I calculate vel_squared                      " << vel_squared << endl;
      double GaussSum = prefactor*exp(-vel_squared*FourLambdaTauInv);
//cerr << "-vel^2/4LT is " << -vel_squared*FourLambdaTauInv << endl;
      RotK -= log(GaussSum);    
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//cerr << "I'm returning kinetic action " << RotK << endl;
//cerr << "*************************************" << endl;
  return (RotK);
}

// ORIGINAL
/*
double TIP5PWaterClass::FixedAxisAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level)
{
  double R = O_H_moment_arm;
//cerr << "ok here we go: R is " << R << endl;
  double RotK = 0.0;
//cerr << "active is " << activeParticles << endl;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl1 = activeParticles(ptclIndex);
    int ptcl2 = ptcl1 + numMol;
    double FourLambdaTauInv=1.0/(4.0*lambda_p*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
//cerr << "ptcl1 is " << ptcl1 << " and ptcl2 is " << ptcl2 << " at slice " << slice << endl;
      // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << O << " and " << Oprime << endl;
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
//cerr << "now WRT COM coords are " << endl;
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
//cerr << "n " << n << endl;
//cerr << "nprime " << nprime << endl;
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
//cerr << "Well I calculate that the angle between n and nprime is " << acos(dotprod(n,nprime,1.0)) << endl;
      //if (n == nprime){
      if (dVecsEqual(n,nprime)){
        theta = 0.0;
        r = z1; 
//cerr << "axis of rotation is z (theta = 0)" << r << endl;
      }
      else{
        theta = GetAngle(n,nprime);
        r = Normalize(CrossProd(n,nprime));
//cerr << "axis of rotation is r" << r << endl;
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
//cerr << "theta " << theta << endl;
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
      double alpha = HOH_half_angle;
      double SinAlpha = sin(alpha);
      double CosAlpha = cos(alpha);
      double phi1 = GetAngle(z1,r);
      double psi1 = GetAngle(z1prime,r);
//cerr << "P1 " << P1 << endl;
//cerr << "n " << n << endl;
//cerr << "P2 " << P2 << endl;
      double phi = phi1;
      double psi = psi1;
//cerr << "phi " << phi << endl;
//cerr << "psi " << psi << endl;
//      if (phi = 3.14159)
        phi = 0.0;
      if (psi = 3.14159)
        psi = 0.0;
//
//cerr << "phi " << phi << endl;
//cerr << "psi " << psi << endl;

// This getup selects the smaller angle, "mechanically" allowing proton exchange, although it doesn't do anything about parity i.e. sign flips so it's kind of wrong.
//    if(phi1<phi2){
        phi = phi1;
      }
      else{
        phi = phi2;
      }
//cerr << "I chose psi " << psi << endl;
      if(psi1<psi2){
        psi = psi1;
      }
      else{
        psi = psi2;
      }
//


      // Calculate angle of rotation (psi)

      if (phi == 0 && psi == 0 && theta == 0){
//cerr << "___________________________________________________zeroing___________" << endl;
        vel_squared = 0;
        prefactor = 1;
      }
      else{
//cerr << "FIXED AXIS TIP" << endl;
        // Calculate quaternions
        double chi = cos(theta/2)*cos((psi + phi)/2);
        double eta = sin(theta/2)*cos((psi - phi)/2);
        double xi = sin(theta/2)*sin((psi - phi)/2);
        double zeta = cos(theta/2)*sin((psi + phi)/2);
//cerr << "chi " << chi << " eta " << eta << " xi " << xi << " zeta " << zeta << endl;
        // Calculate Fixed-Axis Approximation parameters
        double gamma = 2*acos(chi);
        double SinHalfGamma = sin(gamma/2);
//cerr << "gamma " << gamma << " and SinHalfGamma " << SinHalfGamma << endl;
        // Rotation Moment Matrix
        double Lx = R*R;
        double Ly = R*R*CosAlpha*CosAlpha;
        double Lz = R*R*SinAlpha*SinAlpha;
        // Axis of Rotation;
        double ex = eta/SinHalfGamma;
        double ey = -xi/SinHalfGamma;
        double ez = zeta/SinHalfGamma;
        // Rotation Matrix Elements <e|L|e>
        double Mx = Lx*Lx*ex;
        double My = Ly*Ly*ey;
        double Mz = Lz*Lz*ez;
        double MSum = Mx + My + Mz;
        vel_squared = MSum*gamma*gamma;
        prefactor = gamma/(2*SinHalfGamma);
      }
//cerr << "from which I calculate vel_squared                      " << vel_squared << endl;

      double GaussProd = 1.0;
//    for (int dim=0; dim<NDIM; dim++) {
//  	int NumImage=1;
      double GaussSum=0.0;
//	for (int image=-NumImage; image<=NumImage; image++) {
//	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
      GaussSum += prefactor*exp(-vel_squared*FourLambdaTauInv);
//      }

      GaussProd *= GaussSum;
//      }
      RotK -= log(GaussProd);    
    //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//cerr << "I'm returning kinetic action " << RotK << endl;
//cerr << "*************************************" << endl;
  return (RotK);
}
*/


double TIP5PWaterClass::FixedAxisEnergy(int startSlice, int endSlice, int level)
{
  double R = O_H_moment_arm;
  double spring = 0.0;
  double levelTau=ldexp(Path.tau, level);
  double dt = levelTau*hbar;
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double lambda = lambda_p;
  double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  int startparticle = 3*numMol;
  int endparticle = 4*numMol;
  for (int slice=startSlice; slice<endSlice; slice+=skip) {
    double prod = 1.0;
    for (int ptcl1=startparticle; ptcl1<endparticle; ptcl1++) {
      int ptcl2 = ptcl1 + numMol;
      int speciesNum  = Path.ParticleSpeciesNum(ptcl1);
      if (speciesNum == PathData.Path.SpeciesNum("p")){
        // Load coords and their corresponding oxygens (COMs)
        dVec P1 = PathData.Path(slice,ptcl1);
        dVec P2 = PathData.Path(slice,ptcl2);
        dVec P1prime = PathData.Path(slice+skip,ptcl1);
        dVec P2prime = PathData.Path(slice+skip,ptcl2);
        int Optcl = FindCOM(ptcl1);
        dVec O = PathData.Path(slice,Optcl);
        dVec Oprime = PathData.Path(slice+skip,Optcl);
        // Redefine coordinates WRT COM
        P1 -= O;
        P2 -= O;
        P1prime -= Oprime;
        P2prime -= Oprime;
        P1 = Normalize(P1);
        P2 = Normalize(P2);
        P1prime = Normalize(P1prime);
        P2prime = Normalize(P2prime);
        // Calculate bisectors for each configuration
        dVec n = GetBisector(P1,P2);
        dVec nprime = GetBisector(P1prime,P2prime);
        n = Normalize(n);
        nprime = Normalize(nprime);
        double vel_squared;
        double prefactor;
        double theta;
        dVec r;
        double CDsqrt;
      //  if (n == nprime){
        if (dVecsEqual(n,nprime)){
          r = Normalize(CrossProd(P1,P2));
        }
        else{
          r = Normalize(CrossProd(n,nprime));
        }
        theta = GetAngle(n,nprime);
        double alpha = HOH_half_angle;
        double SinAlpha = sin(alpha);
        double CosAlpha = cos(alpha);
        dVec u1 = Normalize(P1 - Scale(n,CosAlpha));
        dVec z1 = Normalize(CrossProd(P1,n));
        dVec u2 = Normalize(P2 - Scale(n,CosAlpha));
        dVec z2 = Normalize(CrossProd(P2,n));
        dVec u1prime = Normalize(P1prime - Scale(nprime,CosAlpha));
        dVec z1prime = Normalize(CrossProd(P1prime,nprime));
        dVec u2prime = Normalize(P2prime - Scale(nprime,CosAlpha));
        dVec z2prime = Normalize(CrossProd(P2prime,nprime));
        double phi1 = GetAngle(z1,r);
        double phi2 = GetAngle(z2,r);
        double psi1 = GetAngle(z1prime,r);
        double psi2 = GetAngle(z2prime,r);
        double phi = phi1;
        double psi = psi1;
/*        if (phi == 3.14159)
          phi = 0.0;
        if (psi == 3.14159)
          psi = 0.0;
*/
// This getup selects the smaller angle, "mechanically" allowing proton exchange, although it doesn't do anything about parity i.e. sign flips so it's kind of wrong.
/*    if(phi1<phi2){
        phi = phi1;
      }
      else{
        phi = phi2;
      }
//cerr << "I chose psi " << psi << endl;
      if(psi1<psi2){
        psi = psi1;
      }
      else{
        psi = psi2;
      }
*/
        // Calculate angle of rotation (psi)
        if (phi == 0 && psi == 0 && theta ==0){
          vel_squared = 0;
          prefactor = 1;
          CDsqrt = 1;
        }
        else{
          double wx = theta/dt;
          double wy = phi/dt;
          double wz = psi/dt;
          double w2 = wx*wx + wy*wy + wz*wz;
/* THESE PARAMTERS ASSUME CONSTANT ANGULAR VELOCITY
        double chi = 1 - w2*dt*dt/8; // not sure about this one: omega dot omega_not??
        double xi = -0.5*wy*dt + 0.25*(0.5*w2*wy)*(dt*dt*dt/6);
        double eta = 0.5*wx*dt + 0.25*(-0.5*w2*wx)*(dt*dt*dt/6);
        double zeta = 0.5*wz*dt + 0.25*(-0.5*w2*wz)*(dt*dt*dt/6);

          double chi = 1 - w2*dt*dt/24; // not sure about this one: omega dot omega_not??
          double xi = -0.75*wy*dt + 0.25*(0.5*w2*wy - 2*wy/(dt*dt))*(dt*dt*dt/6);
          double eta = 0.75*wx*dt + 0.25*(-0.5*w2*wx + 2*wx/(dt*dt))*(dt*dt*dt/6);
          double zeta = 0.25*wz*dt + 0.25*(wx*wz/dt - wx*wy/dt - 0.5*w2*wz + 2*wz/(dt*dt))*(dt*dt*dt/6);
//cerr << "chi " << chi << " eta " << eta << " xi " << xi << " zeta " << zeta << endl;
          // Calculate Fixed-Axis Approximation parameters
*/
        double theta1 = 0;
	double theta2 = theta;
	double phi1= -phi;
	double phi2 = 0;
	double psi1 = 0;
	double psi2 = psi;
        double chi1 = cos(theta1/2)*cos((psi1 + phi1)/2);
        double eta1 = sin(theta1/2)*cos((psi1 - phi1)/2);
        double xi1 = sin(theta1/2)*sin((psi1 - phi1)/2);
        double zeta1 = cos(theta1/2)*sin((psi1 + phi1)/2);
        double chi2 = cos(theta2/2)*cos((psi2 + phi2)/2);
        double eta2 = sin(theta2/2)*cos((psi2 - phi2)/2);
        double xi2 = sin(theta2/2)*sin((psi2 - phi2)/2);
        double zeta2 = cos(theta2/2)*sin((psi2 + phi2)/2);
	double chi = chi2 - chi1;
	double xi = xi2-xi1;
	double zeta = zeta2 - zeta1;
	double eta = eta2 - eta1;
          double gamma = 2*acos(chi);
          double SinHalfGamma = sin(gamma/2);
//cerr << "gamma " << gamma << " and SinHalfGamma " << SinHalfGamma << endl;
          // Rotation Moment Matrix
          double Lx = 2*R*R;
          double Ly = 2*R*R*SinAlpha*SinAlpha;
          double Lz = 2*R*R*SinAlpha*SinAlpha;
          // Axis of Rotation;
          double ex = eta/SinHalfGamma;
          double ey = -xi/SinHalfGamma;
          double ez = zeta/SinHalfGamma;
          // Rotation Matrix Elements <e|L|e>
          double Mx = ex*Lx*ex;
          double My = ey*Ly*ey;
          double Mz = ez*Lz*ez;
          double MSum = Mx + My + Mz;
          vel_squared = MSum*gamma*gamma;
          prefactor = gamma/(2*SinHalfGamma);
          CDsqrt = sqrt(Lx*Ly*Lz)*gamma/(4*sqrt(2.0)*pow((lambda*levelTau),1.5)*SinHalfGamma);

// ORIGINAL
/*
          // Calculate quaternions
          double chi = cos(theta/2)*cos((psi + phi)/2);
          double eta = sin(theta/2)*cos((psi - phi)/2);
          double xi = sin(theta/2)*sin((psi - phi)/2);
          double zeta = cos(theta/2)*sin((psi + phi)/2);
          // Calculate Fixed-Axis Approximation parameters
          double gamma = 2*acos(chi);
          double SinHalfGamma = sin(gamma/2);
          // Rotation Moment Matrix
          double Lx = R*R;
          double Ly = R*R*CosAlpha*CosAlpha;
          double Lz = R*R*SinAlpha*SinAlpha;
          // Axis of Rotation;
          double ex = eta/SinHalfGamma;
          double ey = -xi/SinHalfGamma;
          double ez = zeta/SinHalfGamma;
          // Rotation Matrix Elements <e|L|e>
          double Mx = Lx*Lx*ex;
          double My = Ly*Ly*ey;
          double Mz = Lz*Lz*ez;
          double MSum = Mx + My + Mz;
          vel_squared = MSum*gamma*gamma;
          prefactor = gamma/(2*SinHalfGamma);
          CDsqrt = sqrt(Lx*Ly*Lz)*gamma/(4*sqrt(2.0)*pow((lambda*levelTau),1.5)*SinHalfGamma);
*/
        }

	double d2overFLT = vel_squared*FourLambdaTauInv;
	double expPart = exp(-d2overFLT);
        prod = (1.5 - d2overFLT)*CDsqrt/levelTau*expPart; 
//cerr << "d2overFLT " << d2overFLT << endl;
//cerr << "expPart " << expPart << endl;
//cerr << "prod " << prod << endl;
        Z += expPart;
//cerr << "now Z is " << Z << endl;
    spring += prod; 
//cerr << "spring " << spring << endl;
      }//conditional
    }//loop over particle
  }//loop over slice
  spring = spring/(endSlice*Z);
  //cerr << "returning spring ENERGY = " << spring << endl;
  return spring;
}

// This is the newly modifed version from staticPIMC++ -- see log for 20050711
double TIP5PWaterClass::NewRotKinAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level)
{
  double R = O_H_moment_arm;
  double lambda = lambda_p;
  double RotK = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl1 = activeParticles(ptclIndex);
    int ptcl2 = ptcl1 + numMol;
    double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
//cerr << "n " << n << endl;
//cerr << "nprime " << nprime << endl;
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
      //if (n == nprime){
      if (dVecsEqual(n,nprime)){
        r = z1; 
        theta = 0.0;
      }
      else{
        r = Normalize(CrossProd(n,nprime));
        theta = GetAngle(n,nprime);
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
      double alpha = HOH_half_angle;
      double SinAlpha = sin(alpha);
      double CosAlpha = cos(alpha);
      double phi = GetAngle(z1,r);
      double psi = GetAngle(z1prime,r);

      double lpsi = R*SinAlpha;
      double I_theta = R*R*(cos(phi)*cos(phi) + CosAlpha*CosAlpha*sin(phi)*sin(phi));
      double vel_theta_squared = R*R*theta*theta;
//cerr << "lpsi " << lpsi << endl;
        double vel_psi_squared = lpsi*lpsi*(psi*psi + phi*phi);
//cerr << "vel_psi_sq " << vel_psi_squared << endl;
        vel_squared = vel_psi_squared + vel_theta_squared;
        //vel_squared = vel_psi_squared;//vel_theta_squared;
//cerr << "from which I calculate vel_squared                      " << vel_squared << endl;

      double GaussProd = 1.0;
      double GaussSum=0.0;
      GaussSum += exp(-vel_squared*FourLambdaTauInv);

      GaussProd *= GaussSum;
      RotK -= log(GaussProd);    
      //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//cerr << "I'm returning kinetic action " << RotK << endl;
//cerr << "*************************************" << endl;
  return (RotK);
}

double TIP5PWaterClass::NewRotKinEnergy(int startSlice, int endSlice, int level)
{
  double R = O_H_moment_arm;
  double spring=0.0;
  double levelTau=ldexp(Path.tau, level);
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double lambda = lambda_p;
  double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  int startparticle = 3*numMol;
  int endparticle = 4*numMol;
  for (int ptcl1=startparticle; ptcl1<endparticle; ptcl1++) {
    int ptcl2 = ptcl1 + numMol;
    int speciesNum = Path.ParticleSpeciesNum(ptcl1);
    if (speciesNum == PathData.Path.SpeciesNum("p")){
      for (int slice=startSlice; slice<endSlice; slice+=skip) {
        // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
      //if (n == nprime){
      if (dVecsEqual(n,nprime)){
        r = z1; 
        theta = 0.0;
      }
      else{
        r = Normalize(CrossProd(n,nprime));
        theta = GetAngle(n,nprime);
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
        double alpha = HOH_half_angle;
        double SinAlpha = sin(alpha);
        double CosAlpha = cos(alpha);
        double phi = GetAngle(z1,r);
        double psi = GetAngle(z1prime,r);

        double lpsi = R*SinAlpha;
        double I_theta = R*R*(cos(phi)*cos(phi) + CosAlpha*CosAlpha*sin(phi)*sin(phi));
        double vel_theta_squared = R*R*theta*theta;
//cerr << "lpsi " << lpsi << endl;
        double vel_psi_squared = lpsi*lpsi*(psi*psi + phi*phi);
//cerr << "vel_psi_sq " << vel_psi_squared << endl;
        vel_squared = vel_psi_squared + vel_theta_squared;
        double CDsqrt = R*lpsi*lpsi/(pow((lambda*levelTau),1.5));

        double GaussSum;
        double numSum;
        double d2overFLT = vel_squared*FourLambdaTauInv;
        double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = (1.5 - d2overFLT)/levelTau*CDsqrt*expPart; 
        Z += GaussSum;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  cerr << "returning spring = " << spring << endl;
  return spring;
}

// The old version; not really sure what it's supposed to do.
/*
double TIP5PWaterClass::NewRotKinAction(int startSlice, int endSlice, const Array<int,1> &activeParticles, int level)
{
  double R = O_H_moment_arm;
  double lambda = lambda_p;
//cerr << "ok here we go: R is " << R << endl;
  double RotK = 0.0;
//cerr << "active is " << activeParticles << endl;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl1 = activeParticles(ptclIndex);
    int ptcl2 = ptcl1 + numMol;
    double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
//cerr << "ptcl1 is " << ptcl1 << " and ptcl2 is " << ptcl2 << " at slice " << slice << endl;
      // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << O << " and " << Oprime << endl;
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
//cerr << "now WRT COM coords are " << endl;
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
//cerr << "n " << n << endl;
//cerr << "nprime " << nprime << endl;
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
//cerr << "Well I calculate that the angle between n and nprime is " << acos(dotprod(n,nprime,1.0)) << endl;
      //if (n == nprime){
      if (dVecsEqual(n,nprime)){
        r = z1; 
        theta = 0.0;
//cerr << "axis of rotation is z (theta = 0)" << r << endl;
      }
      else{
        r = Normalize(CrossProd(n,nprime));
        theta = GetAngle(n,nprime);
//cerr << "axis of rotation is r" << r << endl;
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
//cerr << "theta " << theta << endl;
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
      double alpha = HOH_half_angle;
      double SinAlpha = sin(alpha);
      double CosAlpha = cos(alpha);
      double phi = GetAngle(z1,r);
      double psi = GetAngle(z1prime,r);

      double lpsi = R*SinAlpha;
      double I_theta = R*R*(cos(phi)*cos(phi) + CosAlpha*CosAlpha*sin(phi)*sin(phi));
      double vel_theta_squared = I_theta*theta*theta;
//cerr << "lpsi " << lpsi << endl;
        double vel_psi_squared = lpsi*lpsi*(psi*psi);
//cerr << "vel_psi_sq " << vel_psi_squared << endl;
        vel_squared = vel_psi_squared + vel_theta_squared;
        //vel_squared = vel_psi_squared;//vel_theta_squared;
//cerr << "from which I calculate vel_squared                      " << vel_squared << endl;

      double GaussProd = 1.0;
      double GaussSum=0.0;
      GaussSum += exp(-vel_squared*FourLambdaTauInv);

      GaussProd *= GaussSum;
      RotK -= log(GaussProd);    
      //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//cerr << "I'm returning kinetic action " << RotK << endl;
//cerr << "*************************************" << endl;
  return (RotK);
}

double TIP5PWaterClass::NewRotKinEnergy(int startSlice, int endSlice, int level)
{
  double R = O_H_moment_arm;
  double spring=0.0;
  double levelTau=ldexp(Path.tau, level);
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double lambda = lambda_p;
  double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
  int TotalNumParticles = Path.NumParticles();
  int numMol = TotalNumParticles/5;
  int startparticle = 3*numMol;
  int endparticle = 4*numMol;
  for (int ptcl1=startparticle; ptcl1<endparticle; ptcl1++) {
    int ptcl2 = ptcl1 + numMol;
//cerr << "I'm working with particles " << ptcl1 << " and " << ptcl2 << endl;
    int speciesNum = Path.ParticleSpeciesNum(ptcl1);
    if (speciesNum == PathData.Path.SpeciesNum("p")){
      for (int slice=startSlice; slice<endSlice; slice+=skip) {
//cerr << "slice " << slice << " -- " << slice+skip << endl;
        // Load coords and their corresponding oxygens (COMs)
      dVec P1 = PathData.Path(slice,ptcl1);
      dVec P2 = PathData.Path(slice,ptcl2);
      dVec P1prime = PathData.Path(slice+skip,ptcl1);
      dVec P2prime = PathData.Path(slice+skip,ptcl2);
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      int Optcl = FindCOM(ptcl1);
      dVec O = PathData.Path(slice,Optcl);
      dVec Oprime = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << O << " and " << Oprime << endl;
      // Redefine coordinates WRT COM
      P1 -= O;
      P2 -= O;
      P1prime -= Oprime;
      P2prime -= Oprime;
      P1 = Normalize(P1);
      P2 = Normalize(P2);
      P1prime = Normalize(P1prime);
      P2prime = Normalize(P2prime);
//cerr << "now WRT COM coords are " << endl;
//cerr << "P1 " << P1 << endl;
//cerr << "P2 " << P2 << endl;
//cerr << "P1prime " << P1prime << endl;
//cerr << "P2prime " << P2prime << endl;
      // Calculate bisectors for each configuration
      dVec n = Normalize(GetBisector(P1,P2));
      dVec nprime = Normalize(GetBisector(P1prime,P2prime));
//cerr << "n " << n << endl;
//cerr << "nprime " << nprime << endl;
      double vel_squared;
      double prefactor;
      double theta;
      dVec z1 = Normalize(CrossProd(P1,P2));
      dVec z1prime = Normalize(CrossProd(P1prime,P2prime));
      dVec r;
//cerr << "Well I calculate that the angle between n and nprime is " << acos(dotprod(n,nprime,1.0)) << endl;
      //if (n == nprime){
      if (dVecsEqual(n,nprime)){
        r = z1; 
        theta = 0.0;
//cerr << "axis of rotation is z (theta = 0)" << r << endl;
      }
      else{
        r = Normalize(CrossProd(n,nprime));
        theta = GetAngle(n,nprime);
//cerr << "axis of rotation is r" << r << endl;
      }
      // Calculate the cross product - the axis of rotation
      // Calculate polar angles and trig functions
      // Calculate azimuthal angle
//cerr << "theta " << theta << endl;
      // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
        double alpha = HOH_half_angle;
        double SinAlpha = sin(alpha);
        double CosAlpha = cos(alpha);
        double phi = GetAngle(z1,r);
        double psi = GetAngle(z1prime,r);

        double lpsi = R*SinAlpha;
        double I_theta = R*R*(cos(phi)*cos(phi) + CosAlpha*CosAlpha*sin(phi)*sin(phi));
        double vel_theta_squared = I_theta*theta*theta;
//cerr << "lpsi " << lpsi << endl;
        double vel_psi_squared = lpsi*lpsi*(psi*psi);
//cerr << "vel_psi_sq " << vel_psi_squared << endl;
        vel_squared = vel_psi_squared + vel_theta_squared;
        double CDsqrt = R*lpsi*lpsi/(pow((lambda*levelTau),1.5));

        double GaussSum;
        double numSum;
        double d2overFLT = vel_squared*FourLambdaTauInv;
        double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = (1.5 - d2overFLT)/levelTau*CDsqrt*expPart; 
        Z += GaussSum;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  cerr << "returning spring = " << spring << endl;
  return spring;
}
*/


dVec TIP5PWaterClass::CrossProd(dVec v1, dVec v2)
{
  dVec cross;
  cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
  cross[1] = -v1[0]*v2[2] + v1[2]*v2[0];
  cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return cross;
}

double TIP5PWaterClass::Mag(dVec v)
{
  double mag = sqrt(dotprod(v,v,1.0));
  return mag;
}

double TIP5PWaterClass::GetAngle(dVec v1, dVec v2)
{
  double mag = Mag(v1);
  mag *= Mag(v2);
  double dot = dotprod(v1,v2,mag);
  double angle = acos(dot);
//  if (dot > 1.0){
//    cerr << "OH CRAP: DOT PRODUCT IS " << dot << " between " << v1 << " and " << v2 << "; I used mag " << mag << " and I'm going to return " << angle << endl;
//  }
  if (dot-1 < 0.0001 && dot-1 > 0){
//    cerr << "correcting angle" << endl;
    angle = 0.0;
  }
  return angle;
}

dVec TIP5PWaterClass::GetBisector(dVec v1, dVec v2)
{
  dVec bisector = v1 + v2;
  return bisector;
}

dVec TIP5PWaterClass::Normalize(dVec v)
{
  double mag = Mag(v);
  dVec norm = v/mag;
//cerr << "Test normalization: mag of v is " << mag << " and normalized it's " << Mag(norm) << endl;
  return norm;
}

dVec TIP5PWaterClass::Strip(dVec R, dVec u){
  double phi = GetAngle(R,u);
  dVec n = Scale(R,cos(phi));
  return u - n; 
}

dVec TIP5PWaterClass::Scale(dVec v, double scale)
{
  double mag = Mag(v);
  dVec norm = v/mag;
  norm *= scale;
//cerr << "Test scaling: mag of v is " << mag << " and scaled it's" << Mag(norm) << endl;
  return norm;
}

double TIP5PWaterClass::CalcPsi(double theta)
{
  double Alpha = HOH_half_angle;
  double SinAlpha = sin(Alpha);
  double TROUBLE1 = M_PI/2 - Alpha;
  double TROUBLE2 = TROUBLE1 + HOH_angle;
  double psi;
  if (theta < TROUBLE1)
    psi = 0.0;
  else if (theta > TROUBLE2)
    psi = M_PI;
  else
    psi = acos(cos(theta)/SinAlpha);  
  return psi;
}

double TIP5PWaterClass::CalcCutoff(int ptcl1, int ptcl2, int slice, double Rcmag){
  // get oxygen particle ids
  int Optcl1 = FindCOM(ptcl1);
  int Optcl2 = FindCOM(ptcl2);
  // get vectors of oxygens
  dVec O1 = PathData.Path(slice,Optcl1);
  dVec O2 = PathData.Path(slice,Optcl2);
  // get vector between oxygens
  dVec Roo;
  double Ormag;
  PathData.Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Roo);
  // calculate vector that would place O2 at radial cutoff Rcmag
  dVec Rc = Scale(Roo,Rcmag);
  // get constituent coordinates WRT oxygen COM
  dVec P1 = PathData.Path(slice,ptcl1);
  dVec P2 = PathData.Path(slice,ptcl2);
  P1 -= O1;
  P2 -= O2;
  // solve for vector between constituents (if molecule 2 is at cutoff radius)
  dVec R12 = Rc + P2 - P1;
  // obtain cutoff magnitude
  double r12 = Mag(R12);
  return r12;
}


string
TIP5PWaterClass::GetName()
{
  return "TIP5PWater";
}
