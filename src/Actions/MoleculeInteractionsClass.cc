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
#include "MoleculeInteractionsClass.h"
#include "../Moves/MoveUtils.h"

// hack some extra output
//ofstream yout("MolIntAct.dat");
//int bCount;

MoleculeInteractionsClass::MoleculeInteractionsClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  //bCount = 0;
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

	// radial cutoffs for ST2 modulation function
  // these are in angstrom
	RL = 2.0160;
	RU = 3.1287;

  // these are initialized in Read() now
	// SPC/F2 intramolecular potential parameters
	//rho = 2.361;
	//D = 0.708;
	//alpha = 108.0*M_PI/180.0;
	//R_OH_0 = 1.0;
	//R_HH_0 = 2*R_OH_0*sin(alpha/2);
	//b = 1.803;
	//c = -1.469;
	//d = 0.776;
	//// conversion factor for mdyn*angstrom^-1 --> kcal*mol^-1*angstrom^-2
	//Dyn2kcal = 143.929;

	ReadComplete = false;
  TIPPIMC = 1;
}

string MoleculeInteractionsClass::GetName(){
	return("MoleculeInteractionsClass");
}

double MoleculeInteractionsClass::SingleAction (int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level){

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
  else if (startSlice == 0 && endSlice == PathData.Path.TotalNumSlices) {
    startSlice -= 1;
  }
  //cerr << "MolInAct over slices " << startSlice+1 << " " << endSlice-1 << endl;
	bool IsAction = true;
	double TotalU = ComputeEnergy(startSlice+1, endSlice-1, activeParticles, level, TruncateAction, IsAction);
	//cerr << " RETURNING " << TotalU*PathData.Path.tau << " tau " << PathData.Path.tau << endl;
  return(TotalU*PathData.Path.tau);
}

double MoleculeInteractionsClass::d_dBeta (int startSlice, int endSlice,  int level)
{
	//cerr << "MoleculeInteractions::d_dBeta__________________ for slices " << startSlice << " to " << endSlice;// << endl;
	assert(ReadComplete);
  Array<int,1> activeParticles(PathData.Path.NumParticles());
  for (int i=0;i<PathData.Path.NumParticles();i++)
    activeParticles(i)=i;
	
	bool IsAction = false;
	double TotalU = ComputeEnergy(startSlice, endSlice-1, activeParticles, level, TruncateEnergy, IsAction);
	//cerr << " RETURNING " << TotalU << endl;
  return TotalU;
}

double MoleculeInteractionsClass::ComputeEnergy(int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level, bool with_truncations, bool isAction){
  double int1, int2, int3, int4;
  int1 = int2 = int3 = int4 = 0.0;

  //cerr << isAction << " MoleculeInteraction computing SPC energy over slices " << startSlice << " to " << endSlice << " activeP " << activeParticles << endl;
	Updated = false;
  for (int counter=0; counter<Path.DoPtcl.size(); counter++)
    Path.DoPtcl(counter)=true;

  double TotalU = 0.0;
	double TotalLJ = 0.0;
	double TotalCharge = 0.0;
	double TotalSpring = 0.0;
	double TotalKinetic = 0.0;
	double TotalHarmonic = 0.0;
	double TotalPair = 0.0;
	double TotalCore = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;

	//cerr << "Hello.  This is MoleculeInteractionsClass::ComputeEnergy.  ActivePtcls are " << activeParticles << endl;
	//cerr << "I'm doing Intramolecular " << IntraMolecular << " and with_truncations " << with_truncations << " and with modulation " << withS << endl;
	//cerr << "CUTOFF is " << CUTOFF << endl;
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
		//cerr << "DoPtcl bools are after killing " << ptcl1 << Path.DoPtcl << endl;
    int species1=Path.ParticleSpeciesNum(ptcl1);

		//	Calculate Lennard-Jones interaction
    if (Interacting(species1,0)){
			//cerr << "LJ for ptcl " << ptcl1 << " of species " << species1 << endl;
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
				//cerr << "Considering LJ between " << ptcl1 << " and " << ptcl2 << " for which DoPtcl is " << Path.DoPtcl(ptcl2) << "...";
    		if (Interacting(species2,0) && Path.DoPtcl(ptcl2)){
					//cerr << " Going to compute between " << ptcl1 << " and " << ptcl2;
	  			for (int slice=startSlice; slice<=endSlice; slice+=skip){
	    			dVec r;
	    			double rmag;
	    			PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
						//cerr << "... rmag is " << rmag;
						//if(ptcl1==PathData.Mol(ptcl1) && ptcl2==PathData.Mol(ptcl2) && !Updated(ptcl1,ptcl2)){
						//	Updated(ptcl1,ptcl2) = Updated(ptcl2,ptcl1) = true;
						//	COMTable(ptcl1,ptcl2) = COMTable(ptcl2,ptcl1) = rmag;
						//	COMVecs(ptcl1,ptcl2) = r;
						//	COMVecs(ptcl2,ptcl1) = -1*r;
						//}
						//	disregard interactions outside spherical cutoff
            if (rmag <= CUTOFF){
							//cerr << "<"<<CUTOFF;
              double rinv = 1.0/rmag;
              double sigR = PtclSigma(species1)*rinv;
              double sigR6 = sigR*sigR*sigR*sigR*sigR*sigR;
							double offset = 0.0;

							if(with_truncations){
  						  double sigma_over_cutoff = PtclSigma(species1)/CUTOFF;
  						  offset = pow(sigma_over_cutoff,12) - pow(sigma_over_cutoff,6);
							}
              //
							//cerr << "; using offset " << offset << endl;
	      			//double lj = conversion*4*PathData.Species(species1).Epsilon*(sigR6*(sigR6-1) - offset); // this is in kcal/mol 
	      			double lj = conversion*4*PtclEpsilon(species1)*(sigR6*(sigR6-1) - offset); // this is in kcal/mol 
							TotalLJ += lj;
	      			TotalU += lj;

            }
						//cerr << "... outside" << endl;
	  			}
				}
				else{
					//cerr << "Skipped" << endl;
				}
      }
    }
		if(Interacting(species1,1)){
      /// calculating coulomb interactions
      /// calculating charge-charge interactions

  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,1) && Path.DoPtcl(ptcl2)){
					//	don't compute intramolecular interactions
					//  unless told to
					if(IntraMolecular || PathData.Mol(ptcl1)!=PathData.Mol(ptcl2)){
	  				for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    				double Ormag;
	    				dVec Or;
	    				PathData.Path.DistDisp(slice,PathData.Mol(ptcl1),PathData.Mol(ptcl2),Ormag,Or);
              //bool useFirstImage = true;
              /// hack: we're going to sum over the first set of images too
              if(isAction && useFirstImage){
                dVec L = Path.GetBox();
                double tempOrmag = Ormag;
						    for (int x=-1; x<=1; x++) {
						      for (int y=-1; y<=1; y++) {
						        for (int z=-1; z<=1; z++) {
                      Or(0) += x*L(0);
                      Or(1) += y*L(1);
                      Or(2) += z*L(2);
                      Ormag = Mag(Or);
							        // implement spherical cutoff
            	        //double Ormag = COMSeparation(slice,ptcl1,ptcl2);
            	        if (Ormag <= CUTOFF){
                        dVec op1 = Path(slice,ptcl1) - Path(slice,PathData.Mol(ptcl1));
                        dVec op2 = Path(slice,ptcl2) - Path(slice,PathData.Mol(ptcl2));
                        Path.PutInBox(op1);
                        Path.PutInBox(op2);
                        dVec r = Or - op1 + op2;
                        double rmag = Mag(r);
							        	double truncate = 0.0;
                      	double ptclCutoff = 1.0;
							        	if(with_truncations){
							        		truncate = 1.0;
							        		//ptclCutoff = CalcCutoff(ptcl1,ptcl2,slice,CUTOFF);
							        		//ptclCutoff = Mag(r - Or + Scale(Or,CUTOFF));
							        		ptclCutoff = Mag(r - Or + Renormalize(Or,CUTOFF));
							        	}
							        	double modulation=1.0;
							        	if(withS)
							        		modulation = S(Ormag);
                        double coulomb;
                      	//coulomb = conversion*prefactor*PathData.Species(species1).pseudoCharge
							        	//								*PathData.Species(species2).pseudoCharge*modulation
                        //                *(1.0/rmag - truncate/ptclCutoff);
                      	coulomb = conversion*prefactor* PtclCharge(species1)
							        									* PtclCharge(species2) * modulation
                                        * (1.0/rmag - truncate/ptclCutoff);
							        	TotalCharge += coulomb;
	      			        	TotalU += coulomb;
                        //cerr << "  " << ptcl1 << " - " << ptcl2 << " rmag " << rmag << " coulomb " << coulomb << " tmpOrmag " << tempOrmag << endl;
							        	//cerr  << TotalU << " added " << coulomb << " from charge-charge interaction between " << ptcl1 << " and " << ptcl2 << endl;
  	                  }
                    }
                  }
						    }
              }
              // using conventional truncation at spherical cutoff
              else{
                //cerr << "CALLING CONVENTIONAL COULOMB EVALUATION" << endl;
							  // implement spherical cutoff
            	  if (Ormag <= CUTOFF){
                  double rmag;
                  dVec r;
	    				    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
							  	double truncate = 0.0;
                	double ptclCutoff = 1.0;
							  	if(with_truncations){
							  		truncate = 1.0;
							  		ptclCutoff = CalcCutoff(ptcl1,ptcl2,slice,CUTOFF);
							  		//ptclCutoff = Mag(r - Or + Renormalize(Or,CUTOFF));
							  	}
							  	double modulation=1.0;
							  	if(withS)
							  		modulation = S(Ormag);
                  double coulomb;
                	coulomb = conversion*prefactor * PtclCharge(species1) * PtclCharge(species2)
							  									* modulation * (1.0/rmag - truncate/ptclCutoff);
                  //if(ptcl2%(ptcl1-1) == 0)
                  //  cerr << ptcl1 << ", " << ptcl2 << " " << rmag << " " << coulomb << " " << truncate/ptclCutoff << " " << conversion*prefactor*PathData.Species(species1).pseudoCharge*PathData.Species(species2).pseudoCharge*modulation << endl;
                  //coulomb = prefactor*PathData.Species(species1).pseudoCharge*PathData.Species(species2).pseudoCharge/rmag;
							  	TotalCharge += coulomb;
	      			  	TotalU += coulomb;
                  //cerr << "  " << ptcl1 << " " << ptcl2 << " " << rmag << " " << coulomb << " " << TotalCharge << endl;
                }
              }
  	        }
					}
  	    }
  	  }
  	}
		if(Interacting(species1,2) && (ptcl1 == PathData.Mol(ptcl1))){
			// quadratic intramolecular potential
			// SPC/F2; see Lobaugh and Voth, JCP 106, 2400 (1997)
			double spring = 0.0;
			vector<int> activeP(0);
			activeP.push_back(ptcl1);
			//cerr << " INTRA added " << ptcl1;
  		for (int index2=1; index2<PathData.Mol.MembersOf(ptcl1).size(); index2++){
				int ptcl2 = PathData.Mol.MembersOf(ptcl1)(index2);
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,2) && Path.DoPtcl(ptcl2)){
					activeP.push_back(ptcl2);
					//cerr << " INTRA added " << ptcl2;
				}
			}
	  	for (int slice=startSlice;slice<=endSlice;slice+=skip){
        assert(activeP.size() == 3);
				dVec r;
				double ROH1, ROH2, RHH;
				PathData.Path.DistDisp(slice,activeP[0],activeP[1],ROH1,r);
				PathData.Path.DistDisp(slice,activeP[0],activeP[2],ROH2,r);
				PathData.Path.DistDisp(slice,activeP[1],activeP[2],RHH,r);

				double term1 = rho*rho*D*((ROH1 - R_OH_0)*(ROH1 - R_OH_0) + (ROH2 - R_OH_0)*(ROH2 - R_OH_0));
				double term2 = 0.5*b*(RHH - R_HH_0)*(RHH - R_HH_0);
				double term3 = c*(ROH1 + ROH2 - 2*R_OH_0)*(RHH - R_HH_0);
				double term4 = d*(ROH1 - R_OH_0)*(ROH2 - R_OH_0);
				//cerr << "slice " << slice << " INTRA: ROH1 " << ROH1 << " ROH2 " << ROH2 << " RHH " << RHH << " term1 " << term1 << " term2 " << term2 << " term3 " << term3 << " term4 " << term4 << endl;
				spring = conversion*Dyn2kcal*(term1 + term2 + term3 + term4);
				//outfile << spring << " " << term1 << " " << term2 << " " << term3 << " " << term4 << endl;
        int1 += conversion*Dyn2kcal*term1;
        int2 += conversion*Dyn2kcal*term2;
        int3 += conversion*Dyn2kcal*term3;
        int4 += conversion*Dyn2kcal*term4;
				TotalSpring += spring;
				TotalU += spring;
			}
		}
		if(Interacting(species1,3)){
			// compute kinetic action: imaginary time spring term
			// Based on KineticClass::SingleAction, execpt that displacements are WRT the molecule COM

			//cerr << "species " << species1 << ", ptcl " << ptcl1 << " computing kinetic...";
		  double TotalK = 0.0;
		  int skip = 1<<level;
		  double levelTau = Path.tau* (1<<level);
		  double lambda = lambdas(species1);
		  if (lambda != 0.0){
		    double FourLambdaTauInv=1.0/(4.0*lambda*levelTau*levelTau);
        if(isAction){
		      for (int slice=startSlice-1; slice<= endSlice;slice+=skip) {
		        dVec vel;
				  	//vel = COMVelocity(slice, slice+skip, ptcl1);
				  	vel = PathData.Path.Velocity(slice, slice+skip, ptcl1);
			
		        double GaussProd = 1.0;
		        for (int dim=0; dim<NDIM; dim++) {
				  		double GaussSum=0.0;
				  		for (int image=-NumImages; image<=NumImages; image++) {
				  		  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
				  		  GaussSum += exp(-dist*dist*FourLambdaTauInv);
				  		}
				  		GaussProd *= GaussSum;
		        }
				  	TotalK -= log(GaussProd);
		      }
        } else {
		      for (int slice=startSlice; slice <= endSlice;slice+=skip) {
		        dVec vel;
				  	//vel = COMVelocity(slice, slice+skip, ptcl1);
				  	vel = PathData.Path.Velocity(slice, slice+skip, ptcl1);
			
		        double GaussProd = 1.0;
		        for (int dim=0; dim<NDIM; dim++) {
				  		double GaussSum=0.0;
				  		for (int image=-NumImages; image<=NumImages; image++) {
				  		  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
				  		  GaussSum += exp(-dist*dist*FourLambdaTauInv);
				  		}
				  		GaussProd *= GaussSum;
		        }
				  	TotalK -= log(GaussProd);
		      }
        }
		  }
			TotalKinetic += TotalK;
			TotalU += TotalK;
			//cerr << TotalKinetic << endl;
		}
		if(Interacting(species1,4)){
			// harmonic intermolecular potential
			// hard-wired parameters for ST2 water dimer
			// could be generalized
			double harmonic = 0.0;
			//double omega = 26; // ps^-1
			// for Rossky ST2 in kcal/mol*s^2 m^-2
			//double m_H2O = 0.043265; // ???? this seems off by a factor of 2
      //double m_H2O = 0.02163; // kcal/mol*ps^2 angstrom^-2
      double m_H2O = 0.0060537; // kcal/mol*ps^2 bohr^-2
			// Lobaugh & Voth in amu
			//double m_H2O = 9.0;
			//double R0 = 2.85; // angstrom
			double R0 = 5.39; // bohr
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,4) && Path.DoPtcl(ptcl2)){
	  			for (int slice=startSlice;slice<=endSlice;slice+=skip){
  					dVec COMr;
  					double COMrmag;
  					PathData.Path.DistDisp(slice, ptcl1, ptcl2, COMrmag, COMr);
						//cerr << "Harmonic: COMrmag is " << COMrmag << " between " << ptcl1 << " and " << ptcl2 << " and R0 is " << R0 << endl;
						harmonic = conversion*0.5*m_H2O*omega*omega*(COMrmag - R0)*(COMrmag - R0);
						TotalHarmonic += harmonic;
						TotalU += harmonic;
					}
				}
			}
		}
    // Species-wise pair potential 
		if(Interacting(species1,5)){
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,5)){
      //cerr << "Pair Potential ptcls " << Path.Species(species1).FirstPtcl << " to " << Path.Species(species1).LastPtcl << " activeP " << activeParticles << " DoPtcl " << Path.DoPtcl << endl;
          if(ptcl1 != ptcl2 && Path.DoPtcl(ptcl2)) {
	  			  for (int slice=startSlice; slice<=endSlice; slice+=skip){
              double rmag;
              dVec r;
  	          PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
              //cerr << ptcl1 << " " << ptcl2 << " " << slice << " " << rmag << endl;
              double pair = 0.0;
              if(rmag < pairCutoff) {
                //cerr << "accessing " << species1 << " " << species2 << " " << rmag << endl;
                pair = (*PairVTable(species1,species2))(rmag);
                //cerr << "  " << pair << endl;
              }

              TotalPair += pair;
              TotalU += pair;
            }
          }
        }
      }
    }

		if(Interacting(species1,6)){
  		for (int ptcl2=Path.Species(species1).FirstPtcl; ptcl2<=Path.Species(species1).LastPtcl; ptcl2++){
        if(ptcl1 != ptcl2) {
	  			for (int slice=startSlice; slice<=endSlice; slice+=skip){
            double rmag;
            dVec r;
  	        PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
            double U;
            if(boundary=="zero") {
              if(rmag > StartCore) {
                U = 0.;
              }
              else {
                //U = A + (StartCore-rmag)/StartCore*10*A;
                U = (StartCore-rmag)/StartCore*10*A;
              }
            }

            else if(boundary=="infinity") {
              if(rmag < StartCore) {
                U = 0.;
              }
              else {
                //U = A + (rmag-StartCore)/StartCore*10*A;
                U = (rmag-StartCore)/StartCore*10*A;
              }
            }
            //cerr << ptcl1 << " " << ptcl2 << " " << slice << " " << boundary << " " << rmag << " " << U << endl;

            TotalCore += U;
            TotalU += U;
          }
        }
      }
    }
	}
  //cerr << "MolINtAct returning energy " << TotalU << " and action " << TotalU*PathData.Path.tau << endl;// << " tau is " << PathData.Path.tau << " " << Path.tau << endl;
  //cerr << "  It consists of LJ " << TotalLJ << " and coulomb " << TotalCharge << endl;
	//cerr << "Empirical potential contributions: " << TotalU << " " << TotalLJ << " " << TotalCharge << " " << TotalSpring << " " << TotalHarmonic << " " << TotalKinetic;// << endl;
	// write additional output file

  //cerr << "MI " << TotalU << " " << TotalLJ << " " << TotalCharge << " " << TotalSpring << " " << TotalHarmonic << " " << TotalKinetic << " " << TotalPair << " " << endl;
	if(!isAction && special){
		outfile << TotalU << " " << TotalLJ << " " << TotalCharge << " " << TotalSpring << " " << TotalHarmonic << " " << TotalKinetic << " " << TotalPair << " " << int1 << " " << int2 << " " << int3 << " " << int4 << endl;
	}
  return (TotalU);
}


double MoleculeInteractionsClass::CalcCutoff(int ptcl1, int ptcl2, int slice, double Rcmag){
  // get oxygen particle ids
  int Optcl1 = PathData.Mol(ptcl1);
  int Optcl2 = PathData.Mol(ptcl2);
  // get vectors of oxygens
  dVec O1 = Path(slice,Optcl1);
  dVec O2 = Path(slice,Optcl2);
  // get vector between oxygens
  dVec Roo;
  double Ormag;
	//if(!Updated(Optcl1,Optcl2)){
  	Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Roo);
	//	Updated(Optcl1,Optcl2) = Updated(Optcl2,Optcl1) = true;
	//	COMTable(Optcl1,Optcl2) = COMTable(Optcl2,Optcl1) = Ormag;
	//	COMVecs(Optcl1,Optcl2) = Roo;
	//	COMVecs(Optcl2,Optcl1) = -1*Roo;
	//} else {
	//	Ormag = COMTable(Optcl1,Optcl2);
	//	Roo = COMVecs(Optcl1,Optcl2);
	//}
  //dVec Rc = Scale(Roo,Rcmag);
  dVec Rc = Renormalize(Roo,Rcmag);
  // get constituent coordinates WRT oxygen COM
  dVec P1 = Path(slice,ptcl1);
  dVec P2 = Path(slice,ptcl2);
  P1 -= O1;
  Path.PutInBox(P1);
  P2 -= O2;
  Path.PutInBox(P2);
  // solve for vector between constituents (if molecule 2 is at cutoff radius)
  dVec R12 = Rc + P2 - P1;
  // obtain cutoff magnitude
  double r12 = Mag(R12);
  //if(ptcl2%(ptcl1-1)==0) {
  //  cerr << "CCO " << ptcl1 << " " << ptcl2 << " " << r12 << " " << Optcl1 << " " << Optcl2 << " " << O1 << " " << O2 << endl;
  //}
  return r12;
}

double MoleculeInteractionsClass::COMSeparation (int slice,int ptcl1,int ptcl2)
{
  // get oxygen particle ids
  int Optcl1 = PathData.Mol(ptcl1);
  int Optcl2 = PathData.Mol(ptcl2);
  dVec Or;
  double Ormag;

	//if(!Updated(Optcl1,Optcl2)){
  	PathData.Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Or);
//		Updated(Optcl1,Optcl2) = Updated(Optcl2,Optcl1) = true; COMTable(Optcl1,Optcl2) = COMTable(Optcl2,Optcl1) = Ormag;
//		COMVecs(Optcl1,Optcl2) = Or;
//		COMVecs(Optcl2,Optcl1) = -1*Or;
//	} else {
//		Ormag = COMTable(Optcl1,Optcl2);
//	}

  return Ormag;
}

double MoleculeInteractionsClass::S(double r)
{
  double mod;
  if(r<RL)
    mod = 0.0;
  else if(r>RU)
    mod = 1.0;
  else{
    double diff1 = r - RL;
    double diff2 = 3*RU - RL - 2*r;
    double diff3 = RU - RL;
    mod = diff1*diff1*diff2/(diff3*diff3*diff3);
  }
  return mod;
}

dVec MoleculeInteractionsClass::COMVelocity(int sliceA, int sliceB, int ptcl){
	dVec p1 = Path(sliceB, ptcl);
	dVec p2 = Path(sliceA, ptcl);
	if(ptcl != PathData.Mol(ptcl)){
		int COMptcl = PathData.Mol(ptcl);
 		p1 -= Path(sliceB, COMptcl);
    Path.PutInBox(p1);
 		p2 -= Path(sliceA, COMptcl);
    Path.PutInBox(p2);
	}
 	dVec vel = p1 - p2;
  //PathData.Path.PutInBox(vel);
	return vel;
}

void MoleculeInteractionsClass::Read (IOSectionClass &in)
{
	if(!ReadComplete){
		cerr << "In MoleculeInteractionsClass::Read" << endl;
		// DEFAULTS
		prefactor = SI*angstrom_to_m*elementary_charge*
								elementary_charge*N_Avogadro/kcal_to_joule;

    // SPC/F2 intramolecular potential parameters
    // Lobaugh and Voth, JCP 106 2400 (1997)
    // these are the stated parameters from the paper in units of angstrom
	  rho = 2.361;
	  D = 0.708;
	  alpha = 108.0*M_PI/180.0;
	  R_OH_0 = 1.0;
	  R_HH_0 = 2*R_OH_0*sin(alpha/2);
	  b = 1.803;
	  c = -1.469;
	  d = 0.776;
	  // conversion factor for mdyn*angstrom^-1 --> kcal*mol^-1*angstrom^-2
	  Dyn2kcal = 143.929;

    string units = "angstrom";
    in.ReadVar("Units",units);
    if(units == "bohr"){
      prefactor *= bohr_per_angstrom;

      // these are converted into bohr from angstrom
	    rho = 1.24897;
	    D = 1.338;
	    //R_OH_0 = 1.8904;
      // optimized value
	    R_OH_0 = 1.8696;
	    //R_HH_0 = 2*R_OH_0*sin(alpha/2);
      // optimized value
	    R_HH_0 = 2.98;
	    b = 0.9538;
	    c = -0.7771;
	    d = 0.4105;
	    // conversion factor for mdyn*bohr^-1 --> kcal*mol^-1*bohr^-2
	    Dyn2kcal = 76.138;
    }

		CUTOFF = Path.GetBox()(0)/2;
		IntraMolecular = false;
		withS = false;
		TruncateAction = true;
		TruncateEnergy = false;
    useFirstImage = true;

    conversion = 1.0;
		if(in.ReadVar("Prefactor",conversion)){
      cerr << "Setting conversion factor to " << conversion << ". Default units are kcal/mol with length in " << units << endl;
      cerr << "Absolute energy prefactor will be prodcut of " << prefactor <<"*" << conversion <<" " << prefactor*conversion << endl;
    }
		in.ReadVar("Cutoff",CUTOFF);
		in.ReadVar("Modulated",withS);
		in.ReadVar("Intramolecular",IntraMolecular);
		in.ReadVar("TruncateAction",TruncateAction);
		in.ReadVar("TruncateEnergy",TruncateEnergy);
		in.ReadVar("UseImage",useFirstImage);

		Interacting.resize(PathData.NumSpecies(),7);
		Interacting = false;
		LJSpecies.resize(0);
		ChargeSpecies.resize(0);
		SpringSpecies.resize(0);
		KineticSpecies.resize(0);
		QuadSpecies.resize(0);
		PairSpecies.resize(0);
		CoreSpecies.resize(0);
		in.ReadVar("LJSpecies",LJSpecies);
		in.ReadVar("ChargeSpecies",ChargeSpecies);
		in.ReadVar("IntraMolecularSpecies",SpringSpecies);
		in.ReadVar("KineticActionSpecies",KineticSpecies);
		in.ReadVar("QuadraticSpecies",QuadSpecies);
		in.ReadVar("PairActionSpecies",PairSpecies);
		in.ReadVar("CoreSpecies",CoreSpecies);

		cerr << "Read LJSpec " << LJSpecies << " and chargeSpec " << ChargeSpecies << " etc..." << endl;
    PtclEpsilon.resize(PathData.NumSpecies());
    PtclSigma.resize(PathData.NumSpecies());
    for(int s=0; s<PathData.NumSpecies(); s++) {
      PtclEpsilon(s) = PathData.Species(s).Epsilon;
      PtclSigma(s) = PathData.Species(s).Sigma;
    }
    cerr << "Original Species LJ params:" << PtclEpsilon << PtclSigma << endl;
    if(LJSpecies.size() > 0) {
      Array<double,1> setEps, setSig;
      assert(in.ReadVar("ActiveEpsilon",setEps));
      assert(in.ReadVar("ActiveSigma",setSig));
      int epsIndex = 0;
		  for(int s=0; s<LJSpecies.size(); s++) {
			  int myS = Path.SpeciesNum(LJSpecies(s));
        PtclEpsilon(myS) = setEps(epsIndex);
        PtclSigma(myS) = setSig(epsIndex);
        epsIndex++;
      }
    }
    cerr << "Assigned LJ Params for each species:" << PtclEpsilon << PtclSigma << endl;
    PtclCharge.resize(PathData.NumSpecies());
    for(int s=0; s<PathData.NumSpecies(); s++)
      PtclCharge(s) = PathData.Species(s).pseudoCharge;
    cerr << "Original Species charges:" << PtclCharge << endl;
    if(ChargeSpecies.size() > 0) {
      Array<double,1> setCharge;
      assert(in.ReadVar("ActiveCharges",setCharge));
      int chargeIndex = 0;
		  for(int s=0; s<ChargeSpecies.size(); s++) {
			  int myS = Path.SpeciesNum(ChargeSpecies(s));
        PtclCharge(myS) = setCharge(chargeIndex);
        chargeIndex++;
      }
    }
    cerr << "Assigned charges for each species:" << PtclCharge << endl;

		if(KineticSpecies.size() > 0){
			assert(in.ReadVar("Lambdas",lambdas));
			cerr << "Read lambda for each species: " << lambdas << endl;
		}

    if(PairSpecies.size() > 0) {
      //if(PairSpecies.size() > 1) {
      //  cerr << "ERROR MOLECULE PAIR ACTIONS ONLY IMPLEMENTED FOR SINGLE PAIR INTERACTION!" << endl;
      //  assert(0);
      //}
      int numP;
      double start;
      assert(in.ReadVar("NumGridPoints",numP));
      assert(in.ReadVar("StartGrid",start));
      assert(start == 0.0);
      assert(in.ReadVar("EndGrid",pairCutoff));
      grid = new LinearGrid(start, pairCutoff, numP);

      PairVTable.resize(PairSpecies.size());
      for(int s1=0; s1<PairSpecies.size(); s1++) {
        for(int s2=s1; s2<PairSpecies.size(); s2++) {
          string mySpline = "Spline_" + PairSpecies(s1) + "_" + PairSpecies(s2);
          cerr << "Looking for Pair action gridpoint data called " << mySpline << endl;
          Array<double,1> values;
          assert(in.ReadVar(mySpline.c_str(),values));
          assert(values.size() == numP);
          PairVTable(s1,s2) = new CubicSpline(grid, values);
          if(s1 != s2)
            PairVTable(s2,s1) = PairVTable(s1,s2);
        }
      }
      cerr << "PairVTable contains " << PairVTable << endl;
    }

    if(CoreSpecies.size() > 0) {
      if(PairSpecies.size() > 1) {
        cerr << "ERROR CORE INTERACTION ONLY IMPLEMENTED FOR SINGLE PAIR INTERACTION!" << endl;
        assert(0);
      }
      assert(in.ReadVar("Onset",StartCore));
      assert(in.ReadVar("Amplitude",A));
      bool valid=false;
      while(!valid) {
        assert(in.ReadVar("Boundary",boundary));
        if(boundary == "zero") {
          cerr << "Core set between zero and " << StartCore << endl;
          valid = true;
          LO = A;
          HI = 0.;
        } else if (boundary == "infinity") {
          cerr << "Core set between " << StartCore << " and infinity" << endl;
          valid = true;
          LO = 0.;
          HI = A;
        } else {
          cerr << "Boundary " << boundary << " not valid (options are 'zero' or 'infinity')" << endl;
        }
      }
    }
    if(QuadSpecies.size() > 0) {
      omega = 26;
      in.ReadVar("omega",omega);
      cerr << "HARMONIC POT omega " << omega << endl;
    }

		for(int s=0; s<LJSpecies.size(); s++){
			Interacting(Path.SpeciesNum(LJSpecies(s)), 0) = true;
			//cerr << "Setting " << LJSpecies(s) << " to interact via LJ!!" << endl;
		}
		for(int s=0; s<ChargeSpecies.size(); s++)
			Interacting(Path.SpeciesNum(ChargeSpecies(s)), 1) = true;
		for(int s=0; s<SpringSpecies.size(); s++)
			Interacting(Path.SpeciesNum(SpringSpecies(s)), 2) = true;
		for(int s=0; s<KineticSpecies.size(); s++)
			Interacting(Path.SpeciesNum(KineticSpecies(s)), 3) = true;
		for(int s=0; s<QuadSpecies.size(); s++)
			Interacting(Path.SpeciesNum(QuadSpecies(s)), 4) = true;
		for(int s=0; s<PairSpecies.size(); s++)
			Interacting(Path.SpeciesNum(PairSpecies(s)), 5) = true;
		for(int s=0; s<CoreSpecies.size(); s++)
			Interacting(Path.SpeciesNum(CoreSpecies(s)), 6) = true;
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
				
	
		Updated.resize(PathData.Mol.NumMol());
	  COMTable.resize(PathData.Mol.NumMol());
		COMVecs.resize(PathData.Mol.NumMol());

		special = false;
		in.ReadVar("ExtraOutput",special);
		if(special){
			string filename;
			assert(in.ReadVar("File",filename));
			outfile.open(filename.c_str());
			outfile << "# Total LJ Coulomb Intramolecular Quadratic Kinetic Pair Int1 Int2 Int3 Int4" << endl;
		}

    if(PathData.Mol.checkNumParticles != PathData.Path.NumParticles()) {
      cerr << "ERROR: SIZE OF MOLREF DOES NOT MATCH TOTAL NUMBER OF PARTICLES " << PathData.Mol.checkNumParticles << " " << PathData.Path.NumParticles() << endl;
      exit(1);
    }

		ReadComplete = true;
	}
}

dVec MoleculeInteractionsClass::Force(int slice, int ptcl)
{
  Array<int,1> activeP(PathData.Path.NumParticles());
  for (int p=0; p<PathData.Path.NumParticles(); p++)
    activeP(p) = p;
  return(Force(slice, ptcl, activeP));
}

dVec MoleculeInteractionsClass::Force(int slice, int ptcl, int ptcl2)
{
  Array<int,1> activeP(1);
  activeP(0) = ptcl2;
  return(Force(slice, ptcl, activeP));
}

dVec MoleculeInteractionsClass::Force(int slice, int ptcl, Array<int,1>& activePtcls)
{
  if(withS)
    cerr << "WARNING MOLECULEINTERACTIONSCLASS::FORCE NOT BUILT TO INCLUDE ST2 MODULATION" << endl;

  for (int counter=0;counter<Path.DoPtcl.size();counter++)
    Path.DoPtcl(counter)=true;
  Path.DoPtcl(ptcl) = false;
  int species1=Path.ParticleSpeciesNum(ptcl);
  //if(Interacting(species1,2))
  //  cerr << "WARNING MOLECULEINTERACTIONSCLSS::FORCE NOT COMPUTED FOR INTRAMOLECULAR TERM" << endl;

  dVec F; // the force
  F[0] = 0;
  F[1] = 0;
  F[2] = 0;

  int numChangedPtcls = 1;
  //int speciesO=Path.SpeciesNum("O");
  //int speciesp=Path.SpeciesNum("p");
  //int speciese=Path.SpeciesNum("e");

  if (Interacting(species1,0)){
    for (int ptcl2Index=0; ptcl2Index<activePtcls.size(); ptcl2Index++){
      int ptcl2 = activePtcls(ptcl2Index);
      int species2=Path.ParticleSpeciesNum(ptcl2);
      if (Interacting(species2,0) && Path.DoPtcl(ptcl2)
          && PathData.Mol(ptcl)!=PathData.Mol(ptcl2)){
        //cerr << "LJ " << ptcl << " " << ptcl2 << endl;
        dVec r, unitr;
        double rmag;
        PathData.Path.DistDisp(slice, ptcl, ptcl2, rmag, r);
        unitr[0] = r[0]/rmag;
        unitr[1] = r[1]/rmag;
        unitr[2] = r[2]/rmag;
        double rinv = 1.0/rmag;
        //double sigR = PathData.Species(species1).Sigma*rinv;
        double sigR = PtclSigma(species1)*rinv;
        double sigR6 = sigR*sigR*sigR*sigR*sigR*sigR;
        // this is in kcal/mol 
        //double ljmag = conversion*24*PathData.Species(species1).Epsilon*rinv*sigR6*(sigR6-1);
        double ljmag = conversion*24*PtclEpsilon(species1)*rinv*sigR6*(sigR6-1);
        F[0] += ljmag*unitr[0];
        F[1] += ljmag*unitr[1];
        F[2] += ljmag*unitr[2];
      }
    }
  }
  if(Interacting(species1,1)){
    /// calculating charge-charge interactions

    for (int ptcl2Index=0; ptcl2Index<activePtcls.size(); ptcl2Index++){
      int ptcl2 = activePtcls(ptcl2Index);
      int species2=Path.ParticleSpeciesNum(ptcl2);
      if (Interacting(species2,1) && Path.DoPtcl(ptcl2)){
        //cerr << "QQ " << ptcl << " " << ptcl2 << endl;
        //	don't compute intramolecular interactions
        //  unless told to
        if(IntraMolecular || PathData.Mol(ptcl)!=PathData.Mol(ptcl2)){
          double Ormag;
          dVec Or;
          PathData.Path.DistDisp(slice,PathData.Mol(ptcl),PathData.Mol(ptcl2),Ormag,Or);
          double rmag;
          dVec r, unitr;
          PathData.Path.DistDisp(slice,ptcl,ptcl2,rmag,r);
          unitr[0] = r[0]/rmag;
          unitr[1] = r[1]/rmag;
          unitr[2] = r[2]/rmag;
          double coulomb;
          coulomb = conversion * prefactor * PtclCharge(species1) * PtclCharge(species2) / (rmag*rmag);
          //coulomb = conversion*prefactor*PathData.Species(species1).pseudoCharge
          //  *PathData.Species(species2).pseudoCharge/(rmag*rmag);
          F[0] += coulomb*unitr[0];
          F[1] += coulomb*unitr[1];
          F[2] += coulomb*unitr[2];
        }
      }
    }
  }
  if(Interacting(species1,4)){
    // force from harmonic potential
  	double harmonic = 0.0;
  	//double omega = 26; // ps^-1
    double m_H2O = 0.0060537; // kcal/mol*ps^2 bohr^-2
  	double R0 = 5.39; // bohr
  	for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
  		int species2=Path.ParticleSpeciesNum(ptcl2);
  		if (Interacting(species2,4) && Path.DoPtcl(ptcl2)){
        if(PathData.Mol(ptcl)!=PathData.Mol(ptcl2)){
  		    dVec COMr;
  			  double COMrmag;
  			  PathData.Path.DistDisp(slice, ptcl, ptcl2, COMrmag, COMr);
          dVec unitr;
          unitr[0] = COMr[0]/COMrmag;
          unitr[1] = COMr[1]/COMrmag;
          unitr[2] = COMr[2]/COMrmag;
  			  harmonic = -1 * conversion * m_H2O *omega*omega * (COMrmag - R0);
  			  F[0] += harmonic * unitr[0];
  			  F[1] += harmonic * unitr[1];
  			  F[2] += harmonic * unitr[2];
        }
  		}
  	}
  }
  if(Interacting(species1,5)){
  	for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
  		int species2=Path.ParticleSpeciesNum(ptcl2);
  		if (Interacting(species2,5)){
    //cerr << "Pair Potential ptcls " << Path.Species(species1).FirstPtcl << " to " << Path.Species(species1).LastPtcl << " activeP " << activeParticles << " DoPtcl " << Path.DoPtcl << endl;
        if(ptcl != ptcl2 && Path.DoPtcl(ptcl2)) {
          double rmag;
          dVec r;
          PathData.Path.DistDisp(slice, ptcl, ptcl2, rmag, r);
          //cerr << ptcl << " " << ptcl2 << " " << slice << " " << rmag << endl;
          double pair = 0.0;
          if(rmag < pairCutoff) {
            //cerr << "accessing " << species1 << " " << species2 << " " << rmag << endl;
            pair = (*PairVTable(species1,species2)).Deriv(rmag);
            //cerr << "  " << pair << endl;
          }
          dVec unitr;
          unitr[0] = r[0]/rmag;
          unitr[1] = r[1]/rmag;
          unitr[2] = r[2]/rmag;
          F[0] += pair*unitr[0];
          F[1] += pair*unitr[1];
          F[2] += pair*unitr[2];
        }
      }
    }
  }
        //ofstream harm("harmonicForce.dat");
        //double x = 2.0;
        //double dx = 6.0/500;
        //for (int i=0; i<500; i++) {
        //  x += dx;
  			//  harmonic = -1 * conversion * m_H2O *omega*omega * (x - R0);
        //  harm << x << " " << harmonic << endl;
        //}
        //harm.close();
        //cerr << "Done" << endl;
        //exit(1);

  //cerr << "Force " << ptcl << " " << activePtcls << " " << F << endl;
  return F;
}

void MoleculeInteractionsClass::SetNumImages (int num)
{
	NumImages = num;
}
