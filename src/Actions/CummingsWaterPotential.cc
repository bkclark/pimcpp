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
#include "CummingsWaterPotential.h"
#include "../Moves/MoveUtils.h"

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}

CummingsWaterPotentialClass::CummingsWaterPotentialClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  elementary_charge = 1.602*pow(10.0,-19);
  N_Avogadro = 6.022*pow(10.0,23.0);
  kcal_to_joule = 4184;
  epsilon_not = 8.85*pow(10.0,-12);
  angstrom_to_m = pow(10.0,10);
  SI = 1/(4*M_PI*epsilon_not);
  k_B = 1.3807*pow(10.0,-23);
  erg_to_eV = 1.6*pow(10.0,12);
  joule_to_eV = pow(10.0,19)/1.602;
  bohr_per_angstrom = 1.890359;

	ReadComplete = false;
}

string CummingsWaterPotentialClass::GetName(){
	return("CummingsWaterPotentialClass");
}

// global update of tensor T_ij for all molecules i!=j
// called by Solve Dipole (must be called before any call to Update_Ep)
void CummingsWaterPotentialClass::Update_T(int slice) {
  for(int i=0; i<PathData.Mol.NumMol(); i++) {
    //// check integrity of molecules ///////
    //int numMol = PathData.Mol.NumMol();
    //dVec rcheck;
    //double rmagcheck;
	  //PathData.Path.DistDisp(slice, i, i+numMol, rmagcheck, rcheck);
    //if(abs(rmagcheck - 0.27) > 1e-4)
    //  cerr << slice << " mol " << i << " m " << rmagcheck;
	  //PathData.Path.DistDisp(slice, i+numMol, i+2*numMol, rmagcheck, rcheck);
    //if(abs(rmagcheck - 0.9572) > 1e-4)
    //  cerr << slice << " mol " << i << " p1 " << rmagcheck;
	  //PathData.Path.DistDisp(slice, i+numMol, i+3*numMol, rmagcheck, rcheck);
    //if(abs(rmagcheck - 0.9572) > 1e-4)
    //  cerr << slice << " mol " << i << " p2 " << rmagcheck;
    /////////////////////////////////////////
    
    //for(int j=i+1; j<PathData.Mol.NumMol(); j++) {
    for(int j=0; j<PathData.Mol.NumMol(); j++) {
      if(i!=j) {
        dVec r;
        double rmag;
	      PathData.Path.DistDisp(slice, j, i, rmag, r);
        PathData.Path.PutInBox(r);
        double rmag3 = rmag*rmag*rmag;
        Array<double,2> r_direct(3,3);
        for(int m=0; m<3; m++)
          for(int n=0; n<3; n++)
            r_direct(m,n) = r(m) * r(n);
        //cerr << "Update_T " << i <<" " << j << " rdirect " << r_direct << endl;
        double expPart = exp(-1 * rmag*rmag/(4*sigma_COM*sigma_COM));
        double g = erf(rmag/(2*sigma_COM)) - rmag/(sqrt(M_PI) * sigma_COM) * expPart;
        double f = g - rmag3/(6*sqrt(M_PI)*sigma_COM*sigma_COM*sigma_COM) * expPart;
        Array<double,2> G(3,3);
        G = 0.0;
        G(0,0) = g; G(1,1) = g; G(2,2) = g;
        for(int x=0; x<3; x++) {
          for(int y=0; y<3; y++) {
            T(slice,i,j)(x,y) = (f*3*r_direct(x,y)/(rmag*rmag) - G(x,y))/rmag3;
            //T(slice,j,i)(x,y) = (f*3*r_direct(x,y)/(rmag*rmag) - G(x,y))/rmag3;
          }
        }
        //cerr << "Tij " << i << " " << j << " " << T(slice,i,j) << endl;
      }
    }
  }
}

// global update to E field due to induced dipoles
// called by Solve Dipole
// Update_T must be called prior to call to this function
void CummingsWaterPotentialClass::Update_Ep(int slice) {
  Update_RF_p(slice);
  for (int i=0; i<PathData.Mol.NumMol(); i++)
    for(int x=0; x<3; x++)
      //Ep(slice,i)(x) = 0.;
      Ep(slice,i)(x) = RF_p(slice,i)(x)/(conversion*prefactor_Efield); // this unit business is hokey but I think this consistent
  //for (int i=0; i<PathData.Mol.NumMol()-1; i++){
  //  for (int j=i+1; j<PathData.Mol.NumMol(); j++){
  //    for(int m=0; m<3; m++) {
  //      for(int n=0; n<3; n++) {
  //        Ep(slice,i)(m) += T(slice,i,j)(m,n) * p(slice,j)(n);
  //        Ep(slice,j)(m) += T(slice,j,i)(m,n) * p(slice,i)(n);
  //      }
  //    }
  //  }
  //}
  for (int i=0; i<PathData.Mol.NumMol(); i++){
    for (int j=0; j<PathData.Mol.NumMol(); j++){
      if(i != j) {
        //cerr << "Update_Ep T_" << i << j << T(slice,i,j) << endl;
        for(int m=0; m<3; m++) {
          for(int n=0; n<3; n++) {
            Ep(slice,i)(m) += T(slice,i,j)(m,n) * p(slice,j)(n);
          }
        }
      }
    }
  }
  //for (int i=0; i<PathData.Mol.NumMol(); i++)
  //  cerr << "Ep " << i << " " << Ep(slice,i) << endl;
}

void CummingsWaterPotentialClass::Update_RF_q(int slice) {
  // update reaction field permanent dipole, add to Eq
  int numMol = PathData.Mol.NumMol();
  for(int m=0; m<numMol; m++) {
    dVec b;
    b = PathData.Path(slice,m+numMol) - PathData.Path(slice,m);
    b = Normalize(b);
    p_permanent(slice,m) = mu_0 * b;
    RF_q(slice, m) = conversion * prefactor_RF * p_permanent(slice, m);
  }
  for(int l=0; l<numMol; l++) {
    for(int m=l+1; m<numMol; m++) {
      dVec R_lm;
      double rmag;
      PathData.Path.DistDisp(slice,l, m, rmag, R_lm);
      if(rmag < Rcut_RF) {
        RF_q(slice, m) += conversion * prefactor_RF * p_permanent(slice,l);
        RF_q(slice, l) += conversion * prefactor_RF * p_permanent(slice,m);
      }
    }
  }
}

void CummingsWaterPotentialClass::Update_RF_p(int slice) {
  // update reaction field induced dipole
  int numMol = PathData.Mol.NumMol();
  for(int l=0; l<numMol; l++) {
    RF_p(slice, l) = prefactor_RF * p(slice, l);
    for(int m=l+1; m<numMol; m++) {
      dVec R_lm;
      double rmag;
      PathData.Path.DistDisp(slice,l, m, rmag, R_lm);
      if(rmag < Rcut_RF) {
        RF_p(slice, m) += prefactor_RF * p_permanent(slice,l);
        RF_p(slice, l) += prefactor_RF * p_permanent(slice,m);
      }
    }
  }
}

// called by Solve Dipole
void CummingsWaterPotentialClass::Update_Eq(int slice) {
  for(int mol1=0; mol1<PathData.Mol.NumMol(); mol1++)
    Update_Eq(slice, mol1);
}

// particle-wise update to Eq
void CummingsWaterPotentialClass::Update_Eq(int slice, int i) {
  for(int x=0; x<3; x++) {
    // initalize to reaction field contribution
    Eq(slice,i)(x) = RF_q(slice, i)(x);
    //  skip Reaction Field
    //Eq(slice,i)(x) = 0.;
  }
  //int species1=Path.ParticleSpeciesNum(i);
  for(int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++) {
   	int species2=Path.ParticleSpeciesNum(ptcl2);
    // make sure it's a charge-carrying ptcl
    if (Interacting(species2,1) && PathData.Mol(i) != PathData.Mol(ptcl2)) {
      double rmag;
      dVec r;
      PathData.Path.DistDisp(slice,ptcl2,i,rmag,r);
      if(rmag < CUTOFF) {
	    //double Ormag;
	    //dVec Or;
	    //PathData.Path.DistDisp(slice,PathData.Mol(i),PathData.Mol(ptcl2),Ormag,Or);
      //if (Ormag <= CUTOFF){
        //PathData.Path.PutInBox(r);
        double E_2_scalar = conversion * prefactor_Efield  / (rmag*rmag*rmag)
          * PathData.Species(species2).pseudoCharge
          * (erf(rmag / sqrt(2 * (sigma_COM*sigma_COM + PathData.Species(species2).chargeSpread)))
              - sqrt(2) * rmag / sqrt(M_PI * (sigma_COM*sigma_COM + PathData.Species(species2).chargeSpread))
              * exp(-1 * rmag*rmag / (2 * (sigma_COM*sigma_COM + PathData.Species(species2).chargeSpread))));
        dVec E_2;
        for(int x=0; x<3; x++) {
          E_2(x) = E_2_scalar * r(x);
          Eq(slice, i)(x) += E_2(x);
        }
        //cerr << slice << " Eq " << i << ", " << ptcl2 << " " << E_2 << endl;;
      }
    }
  }
  //cerr << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2 " << slice << ", " << i << " Eq " << Eq(slice,i) << endl;
  //cerr << "by the way prefactor " << prefactor_Efield*conversion << endl;
}

// call this function for global updates to E fields and 
// induced dipole moments
// then sum quantities in SingleAction and dBeta calls
void CummingsWaterPotentialClass::SolveDipole(int slice) {
  //cerr << "hack generating com coords" << endl;
  //ofstream out("com.dat");
  //int numMoly = PathData.Mol.NumMol();
  //for(int m=0; m<numMoly; m++) {
  //  dVec O = PathData.Path(0,m);
  //  dVec M = PathData.Path(0,m+numMoly);
  //  dVec om = M-O;
  //  cerr << m << " " << O << M << om << endl;
  //  //PathData.Path.PutInBox(om);
  //  dVec c = Renormalize(om,0.032549015367683058);
  //  cerr << c;
  //  c += O;
  //  PathData.Path.PutInBox(c);
  //  cerr << c << endl;
  //  out << c[0] << ", " << c[1] << ", " << c[2] << ",\n";
  //}
  //out.close();
  //exit(0);

  //cerr << "SolveDipole taking ";
  Update_T(slice);
  // maybe this doesn't have to be a global update?
  // probably doesn't matter compared to Ep update
  Update_RF_q(slice);
  Update_Eq(slice);

  //if(firstTime) {
  //  firstTime = false;
    int numSlices = PathData.Path.NumTimeSlices();
    int numMol = PathData.Mol.NumMol();
  //  for(int s=0; s<numSlices; s++) {
      int s = slice;
      for(int m=0; m<numMol; m++) {
      // initialize p... should there be a better initial value?
        //dVec b;
        //b = PathData.Path(s,m+numMol) - PathData.Path(s,m);
        //PathData.Path.PutInBox(b);
        //b = Normalize(b);
        p(s,m) = alpha*Eq(slice,m)*(1.0/(conversion*prefactor_Efield));
        //p(s,m) = mu_0 * b;
        //p(s,m) = ((0., 0., 0.));
        //cerr << m << " " << b << " " << p(s,m) << endl;
        //cerr << m << " " << PathData.Path(s,m+numMol) << " " << PathData.Path(s,m) << " " << b << " " << p(s,m) << endl;
      }
  //  }
  //}
  
  Array<dVec, 1> p0(PathData.Mol.NumMol());
  Array<dVec, 1> p_last(PathData.Mol.NumMol());
  // zero out dipoles
  //for(int m=0; m<PathData.Mol.NumMol(); m++)
  //  p(slice,m) = ((0.0,0.0,0.0));
  
  // store initial dipoles
  for(int i=0; i<PathData.Mol.NumMol(); i++)
    for(int x=0; x<3; x++)
      p0(i)(x) = p(slice,i)(x);

  double dp_max = 2*tolerance;
  double dp;

  // iterate dipole update and check for convergence
  // go thru loop at least once
  int tries = 0;
  double netP = 0.;
  while(dp_max > tolerance) {
    Update_Ep(slice);
    dp_max = 0;
    for(int i=0; i<PathData.Mol.NumMol(); i++) {
      double pmag = sqrt(dot(p(slice,i),p(slice,i)));
      netP += pmag;
      //if(i==0)
      //cerr << i << " " << tries << " " << p(slice,i) << " " << Eq(slice,i)/prefactor_Efield << " " << Ep(slice,i) << endl;
      for(int x=0; x<3; x++) {
        p_last(i)(x) = p(slice,i)(x);
        p(slice,i)(x) = alpha*(Eq(slice,i)(x)/(conversion*prefactor_Efield) + Ep(slice,i)(x));
        //p(slice,i)(x) = p0(slice,i)(x) + alpha*(Eq(slice,i)(x)/prefactor_Efield + Ep(slice,i)(x));
      }
      dp = sqrt(dot((p(slice,i)-p_last(i)), (p(slice,i)-p_last(i))));
      if(dp > dp_max)
        dp_max = dp;
    }
    tries++;
    //if(tries>100)
    //cerr << tries << " " << dp_max << " " << tolerance << endl;
  }
  if(tries>20)
    cerr << tries << " attempts to converge dipole moments avg netP " << netP/(tries*PathData.Mol.NumMol()) << endl;
  double avg_netP = netP/(tries*PathData.Mol.NumMol());
  if(avg_netP>0.9)
    cerr <<  "netP " << avg_netP << endl;
}

double CummingsWaterPotentialClass::SingleAction (int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level)
{
  //cerr<<"CUMMING SINGLE ACTION BEGINS"<<endl;
  //cerr << "MoleculeInteractions::Action__________________ for slices " << startSlice << " to " << endSlice;// << endl;
  //bCount++;
  //if(bCount%64 == 0)
  //yout << bCount << endl;
	assert(ReadComplete);
  if(startSlice == 0 && endSlice == 0){
    startSlice -= 1;
    endSlice += 1;
  }
  else if (startSlice == 0 && endSlice == PathData.Path.TotalNumSlices-1) {
    startSlice -= 1;
    endSlice += 1;
  }
  //cerr << "MolInAct over slices " << startSlice+1 << " " << endSlice-1 << endl;
	bool IsAction = true;
	double TotalU = ComputeEnergy(startSlice+1, endSlice-1, activeParticles, level, IsAction);
	//cerr << " RETURNING " << TotalU*PathData.Path.tau << " tau " << PathData.Path.tau << endl;
	//cerr<<"SINGLE ACTION ENDS"<<endl;
  return(TotalU*PathData.Path.tau);
}

double CummingsWaterPotentialClass::d_dBeta (int startSlice, int endSlice,  int level)
{
  //cerr << "MoleculeInteractions::d_dBeta__________________ for slices " << startSlice << " to " << endSlice<<endl;// << endl;
	assert(ReadComplete);
  Array<int,1> activeParticles(PathData.Path.NumParticles());
  for (int i=0;i<PathData.Path.NumParticles();i++)
    activeParticles(i)=i;
	
	bool IsAction = false;
  // hack hack hack testing r_OO dependence
  //SetMode(NEWMODE);
  //double dx=0.03;
  //Array<int,1> activeP(5);
  //activeP(0) = 1;
  //activeP(1) = 3;
  //activeP(2) = 5;
  //activeP(3) = 7;
  //activeP(4) = 9;
  ////cerr << "hi activeP "<< activeP << endl;
  //for(int pIndex=0; pIndex<activeP.size(); pIndex++) {
  //  int ptcl=activeP(pIndex);
  //  dVec R = PathData.Path(0, ptcl);
  //  R(0) -= 3;
  //  PathData.Path.PutInBox(R);
  //  PathData.Path.SetPos(0, ptcl, R);
  //}
  //for(int s=0; s<200; s++) {
  //  for(int pIndex=0; pIndex<activeP.size(); pIndex++) {
  //    int ptcl=activeP(pIndex);
  //    dVec R = PathData.Path(0, ptcl);
  //    R(0) += dx;
  //    PathData.Path.PutInBox(R);
  //    PathData.Path.SetPos(0, ptcl, R);
  //    //cerr << s << " " << ptcl << " " << R << PathData.Path(0,ptcl) << endl;
  //  }
	//  ComputeEnergy(0, 0, activeP, level, false);
  //}
  //exit(0);
	double TotalU = ComputeEnergy(startSlice, endSlice-1, activeParticles, level, IsAction);
	//cerr << " RETURNING " << TotalU << endl;
  return TotalU;
}

// very inefficient but we call a global upate of the induced dipoles no matter what
// don't think this can be avoided!
double CummingsWaterPotentialClass::ComputeEnergy(int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level, bool isAction)
{/*
  for (int counter=0; counter<Path.DoPtcl.size(); counter++)
    Path.DoPtcl(counter)=true;

  double TotalU = 0.0;
	double TotalExp6 = 0.0;
	double TotalCharge = 0.0;
	double kspace = 0.0;
	double TotalDipole = 0.0;
  double TotalReactionField = 0.0;
  double checkDi = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;

  // long range kspace initialization
  if(firstTime){
    vector<int> activeSpec(0);
    int numActiveSpec = 0;
    for (int species=0; species<Path.NumSpecies(); species++) {
		  if(Interacting(species,1)) {
        activeSpec.push_back(species);
        numActiveSpec++;
      }
    }
    int kPts = Path.kVecs.size();
    k2.resize(kPts);
    for(int i=0; i<kPts; i++){
      double kSq = dot(Path.kVecs(i), Path.kVecs(i));
      k2(i) = kSq;
    }
    phi.resize(numActiveSpec, numActiveSpec, kPts);
    for(int si=0; si<numActiveSpec; si++) {
      int spec1 = activeSpec[si];
      for(int sj=0; sj<numActiveSpec; sj++) {
        int spec2 = activeSpec[sj];
        double kalpha = 0.5/(PathData.Species(spec1).chargeSpread + PathData.Species(spec2).chargeSpread);
        for(int i=0; i<kPts; i++){
          phi(spec1, spec2, i) = exp(-k2(i)/(4*kalpha*kalpha));
        }
      }
    }
    cerr << "GCPM Long Range initialized phi of size " << phi.size() << endl;
    firstTime= false;
    //cerr << "activespec " << activeSpecies << endl;

    // initialize Rho_k
    for (int slice=0; slice<=Path.TotalNumSlices; slice+=1) {
      for (int species=0; species<Path.NumSpecies(); species++) {
		    if(Interacting(species,1)){
          Path.CalcRho_ks_Fast(slice, species);
        }
      }
    }
    volfactor = 1.0;
    for(int i=0; i<3; i++)
      volfactor *= 1.0/Path.GetBox()(i);
  }
  // update rho_k
  if (GetMode() == NEWMODE)
  {
    for (int slice=startSlice; slice<=endSlice; slice+=skip) {
      for (int species=0; species<Path.NumSpecies(); species++) {
		    if(Interacting(species,1)){
          Path.CalcRho_ks_Fast(slice, species);
        }
      }
    }
  }
  // kspace summation of gaussian charge interactions
  for (int species1=0; species1<Path.NumSpecies(); species1++) {
	  if(Interacting(species1,1)) {
      for (int slice=startSlice; slice<=endSlice; slice+=skip) {
        // First, do the homologous (same species) terms
        double qi = PathData.Species(species1).Charge;
        //cerr << "homo: spec " << species << " with q " << qi << " Path.NumSpec is " << Path.NumSpecies() << endl;
        for (int ki=0; ki<Path.kVecs.size(); ki++) {
	        double rhok2 = mag2(Path.Rho_k(slice,species1,ki));
          double myKterm = prefactor*volfactor * 2*M_PI * qi*qi * rhok2 * phi(species1, species1, ki)/k2(ki);
          if(isnan(myKterm))
            cerr << "NAN " << species1 << " " << ki << " " << mag2 << " " << myKterm << endl;
	        kspace += myKterm;
	      }

        // Now do the heterologous terms
        for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
          if (Interacting(species2,1)){
            double qi = PathData.Species(species1).Charge;
            double qj = PathData.Species(species2).Charge;
	          for (int ki=0; ki<Path.kVecs.size(); ki++) {
	            double rhorho = 
	              Path.Rho_k(slice, species1, ki).real() *
	              Path.Rho_k(slice, species2, ki).real() + 
	              Path.Rho_k(slice, species1, ki).imag() *
	              Path.Rho_k(slice, species2, ki).imag();
	            double myKterm = prefactor*volfactor * 2*M_PI * qi*qj * rhorho * phi(species1, species2, ki)/k2(ki);
              if(isnan(myKterm))
                cerr << "NAN het " << species1 << " " << species2 << " " << ki << " " << mag2 << " " << myKterm << endl;
	            kspace += myKterm;
              //cerr << "hetero " << ki << " qi " << qi << " qj " << qj << " rhorho " << rhorho << " phi_k " << phi(ki) << " k^-2 " << 1.0/k2(ki) << " total " << myKterm << endl;
            }
          }
        }
      }
    }
  }
  TotalU += kspace;

  double rSeparation;
  // update and sum up contributions from induced dipole
	for (int slice=startSlice; slice<=endSlice; slice+=skip) {
    // this call is where everything is updated
    SolveDipole(slice);
    // this is done molecule-wise
    for (int mol=0; mol<PathData.Mol.NumMol(); mol++) {
      dVec netE = -1*Eq(slice,mol) - 0.5*Ep(slice,mol) * prefactor_Efield + 1./(2*alpha)*p(slice,mol) * prefactor_Efield; // eqn 8
      double piece1_BC=0.0;
      double piece2_BC=0.0;
      double piece3_BC=0.0;
      piece1_BC+=dot(p(slice,mol),Eq(slice,mol));
      piece2_BC+=dot(p(slice,mol),0.5*Ep(slice,mol) * prefactor_Efield);
      piece3_BC+=dot(p(slice,mol),1./(2*alpha)*p(slice,mol) * prefactor_Efield);
      checkDi += dot(p(slice,mol), netE); // eqn 8
      TotalDipole -= 0.5 * dot(p(slice,mol), Eq(slice,mol)); // eqn 9s
      //cout << mol << "p " << p(slice,mol) << endl;
      //cout << mol << "Eq " << Eq(slice,mol) << endl;
      // long-range reaction field contribution
      TotalReactionField -= 0.5 * p_permanent(slice, mol) * RF_q(slice, mol);
      //cout << mol << " " << piece1_BC << " " << piece2_BC << " " << piece3_BC << endl;
    }
  }
  TotalU += TotalDipole;
  TotalU += TotalReactionField;

  // Particle-wise terms (Exp6, and for comparison short-range charges)
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);

		//	Calculate exp-6 interaction
    if (Interacting(species1,0)){
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,0) && Path.DoPtcl(ptcl2)){
	  			for (int slice=startSlice; slice<=endSlice; slice+=skip){
	    			dVec r;
	    			double rmag;
	    			PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
						//	disregard interactions outside spherical cutoff
            //cerr << "exp6 " << ptcl1 << " " << ptcl2 << " r " << rmag << endl;
            if (rmag <= CUTOFF){
							//cerr << "<"<<CUTOFF;
              double rinv = 1.0/rmag;
              double sigR = PathData.Species(species1).Sigma*rinv;
              double sigR6 = sigR*sigR*sigR*sigR*sigR*sigR;
							double shift = 0.0;
              double expPart = exp(gamma * (1 - (rmag/PathData.Species(species1).Sigma)));
              
              // should do proper tail correction here!
							if(with_truncations){
  						  double sigma_over_cutoff = PathData.Species(species1).Sigma/CUTOFF;
  						  double shiftExpPart = exp(gamma * (1 - (CUTOFF/PathData.Species(species1).Sigma)));
	      			  shift = conversion / (1 - 6.0/gamma)
                  * PathData.Species(species1).Epsilon
                  * (6*shiftExpPart/gamma - pow(sigma_over_cutoff,6)); // this is in kcal/mol by default

							}
	      			double exp_6 = conversion / (1 - 6.0/gamma)
                * PathData.Species(species1).Epsilon
                * (6*expPart/gamma - sigR6) - shift; // this is in kcal/mol by default
							TotalExp6 += exp_6;
              //cerr << ptcl1 << " " << ptcl2 << " " << rmag << " " << exp_6 << endl;
	      			TotalU += exp_6;
              rSeparation = rmag;
            }
	  			}
				}
      }
    }

    /// calculating electrostatic coulomb interactions
		if(Interacting(species1,1)){
      // original calculations
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,1) && Path.DoPtcl(ptcl2)){
					if(PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  				for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    				//double Ormag;
	    				//dVec Or;
	    				//PathData.Path.DistDisp(slice,PathData.Mol(ptcl1),PathData.Mol(ptcl2),Ormag,Or);
              //if (Ormag <= CUTOFF){
                double rmag;
                dVec r;
                PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
              if(rmag < CUTOFF) {
                double coulomb;
                // we don't double count so factor of 0.5 is not needed
                //coulomb = 0.5 * conversion * prefactor / rmag
                coulomb = conversion * prefactor / rmag
                  * PathData.Species(species1).pseudoCharge
                  * PathData.Species(species2).pseudoCharge
                  * erf(rmag / sqrt(2 * (PathData.Species(species1).chargeSpread + PathData.Species(species2).chargeSpread)));
                //cerr << "charge " << slice << " " << ptcl1 << " " << ptcl2 << " " << Ormag  << " " << rmag << " " << coulomb << endl;
                //cerr << ptcl1 << " " << ptcl2 << " " << rmag << " " << coulomb << endl;
                TotalCharge += coulomb;
                //TotalU += coulomb;
              }
  	        }
					}
  	    }
  	  }
  	}
  }
  //cerr<<"Bryan check: "<<piece1_BC<<" "<<piece2_BC<<" "<<piece3_BC<<" "<< checkDi <<" "<<TotalDipole<<endl;
  
  if(!isAction) {
    //cout << "Dipole Exp6 Charge_short Charge_kspace_total Total" << endl;
    cout << TotalDipole << " " << TotalExp6 << " " << TotalCharge << " " << kspace << " " << TotalU << endl;
    //cout << TotalDipole << "  " << TotalExp6 << " " << TotalCharge << " " << TotalU << " " << rSeparation << endl;
    //cout << "Erf(1) " << erf(1) << " prefactor " << conversion*prefactor << endl;
    }*/
  return 0;//(TotalU);
}

void CummingsWaterPotentialClass::Read (IOSectionClass &in)
{
  cerr<<"CUMMINGS READ BEGIN "<<endl;
	if(!ReadComplete){
		cerr << "In CummingsWaterPotentialClass::Read" << endl;
		// DEFAULTS
		prefactor = SI * angstrom_to_m
      * elementary_charge * elementary_charge
      * N_Avogadro / kcal_to_joule;
		prefactor_Efield = SI * angstrom_to_m
      * elementary_charge * elementary_charge
      * N_Avogadro / kcal_to_joule;


	  // conversion factor for mdyn*angstrom^-1 --> kcal*mol^-1*angstrom^-2
	  Dyn2kcal = 143.929;

    string units = "angstrom";
    in.ReadVar("Units",units);
    if(units == "bohr"){
      prefactor *= bohr_per_angstrom;
      prefactor_Efield *= bohr_per_angstrom;

      // these are converted into bohr from angstrom
	    // conversion factor for mdyn*bohr^-1 --> kcal*mol^-1*bohr^-2
	    Dyn2kcal = 76.138;
    }
    cerr << "UNITS PREFACTOR " << prefactor << endl;
    cerr << "UNITS k_B " << 1.38e-23*N_Avogadro/kcal_to_joule << endl;

		CUTOFF = Path.GetBox()(0)/2;
    Rcut_RF = CUTOFF;

    conversion = 1.0;
		if(in.ReadVar("Prefactor",conversion)){
      cerr << "Setting conversion factor to " << conversion << ". Default units are kcal/mol with length in " << units << " overall prefactor is now " << prefactor*conversion << endl;
    }
		in.ReadVar("Cutoff",CUTOFF);
		assert(in.ReadVar("sigma",sigma_COM));
		assert(in.ReadVar("alpha",alpha));
		assert(in.ReadVar("gamma",gamma));
		assert(in.ReadVar("dipole",mu_0));
		assert(in.ReadVar("Tolerance",tolerance));
		assert(in.ReadVar("dielectric",dielectric_RF));

    prefactor_RF = prefactor * (dielectric_RF-1) * 2 / ((2*dielectric_RF+1) * Rcut_RF*Rcut_RF*Rcut_RF);
		Interacting.resize(PathData.NumSpecies(),2);
		Interacting = false;
		LJSpecies.resize(0);
		ChargeSpecies.resize(0);
		in.ReadVar("Exp6Species",LJSpecies);
		in.ReadVar("ChargeSpecies",ChargeSpecies);

		for(int s=0; s<LJSpecies.size(); s++){
			Interacting(Path.SpeciesNum(LJSpecies(s)), 0) = true;
			//cerr << "Setting " << LJSpecies(s) << " to interact via LJ!!" << endl;
		}
		for(int s=0; s<ChargeSpecies.size(); s++)
			Interacting(Path.SpeciesNum(ChargeSpecies(s)), 1) = true;
		cerr << "MoleculeInteractions::Read I have loaded Interacting table " << Interacting << endl;

		// make sure something is filled
		bool empty = true;
		for(int i=0; i<Interacting.rows(); i++){
			for(int j=0; j<Interacting.cols(); j++){
				if(Interacting(i,j) != 0){
					empty = (empty && false);
				}
			}
		}
		assert(!empty);
				
	
    int numSlices = PathData.Path.NumTimeSlices();
    int numMol = PathData.Mol.NumMol();
    Eq.resize(numSlices, numMol);
    Ep.resize(numSlices, numMol);
    Ep = 0.0;
    p.resize(numSlices, numMol);
    p_permanent.resize(numSlices, numMol);
    RF_p.resize(numSlices, numMol);
    RF_q.resize(numSlices, numMol);
	  T.resize(numSlices, numMol, numMol);
    for(int s=0; s<numSlices; s++) {
      for(int m=0; m<numMol; m++) {
        //// initialize p... should there be a better initial value?
        //dVec b;
        //b = PathData.Path(s,m+numMol) - PathData.Path(s,m);
        //PathData.Path.PutInBox(b);
        ////b = Normalize(b);
        //p(s,m) = mu_0 * b;
        //cerr << m << " " << PathData.Path(s,m+numMol) << " " << PathData.Path(s,m) << " " << b << " " << p(s,m) << endl;
        ////p(s,m) = ((0.0,0.0,0.0));
        for(int n=0; n<numMol; n++)
          T(s,m,n).resize(3,3);
      }
    }
    firstTime = true;
		ReadComplete = true;
	}
  with_truncations = true;
  in.ReadVar("Truncate",with_truncations);
	cerr<<"CUMMINGS READ END "<<endl;
}
