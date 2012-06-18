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

#ifndef WATER_CLASS_H
#define WATER_CLASS_H
// #include "../PathDataClass.h"
#include "ActionBase.h"
#define CO 18.0

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;
#include <fstream>

#include <iomanip>
#include "Timer.h"
#define angs2bohr 1.889725989


class dVecp
{
public:
  vector<double> vec;
  dVecp()
  {
    vec.resize(NDIM);
    for (int i=0;i<NDIM;i++){
       vec[i]=0.0;
    }
    
  }
  
};

class HelpClass
{
 public:
void ZeroVecp(vector<dVecp> &myVec)
{
  for (int i=0;i<myVec.size();i++){
    for (int dim=0;dim<NDIM;dim++){
      myVec[i].vec[dim]=0.0;
    }

  }

}

void ZeroVecp(vector<double> &myVec)
{
  for (int i=0;i<myVec.size();i++){
      myVec[i]=0.0;
  }

}



void ZeroVecp(vector<complex<double> > &myVec)
{
  for (int i=0;i<myVec.size();i++){
      myVec[i]=0.0;
  }

}


void PrintVec(vector<dVecp> &myVec)
{
  cerr<<setprecision (9);
  for (int i=0;i<myVec.size();i++)
    cerr<<myVec[i].vec[0]<<" "<<myVec[i].vec[1]<<" "<<myVec[i].vec[2]<<endl;
  
}
double dot(dVecp &a, dVecp &b)
{
  double ans=0.0;
  for (int dim=0;dim<NDIM;dim++)
    ans+=a.vec[dim]*b.vec[dim];
  return ans;

}

};

// oxygen -> 0 
// hydrogen -> 1


class ParamsClass
{
public:
  double beta;
  double dipfac; 
  double qO;
  double qH;

  double cOO;
  double cOH;
  
  double bOO;
  double bOH;
  double alpha_O;
  double raggio;


  double D1_OO;
  double D1_OH;
  double D1_HH;

  double D2_OO;
  double D2_OH;
  double D2_HH;

  double gamma1_OO;
  double gamma1_OH;
  double gamma1_HH;

  double gamma2_OO;
  double gamma2_OH;
  double gamma2_HH;


  double r1_OO;
  double r1_OH;
  double r1_HH;



  double r2_OO;
  double r2_OH;
  double r2_HH;



  dVecp Box;
  dVecp BoxInv;
  void Init()
  {
    Box.vec[0]= 12.431071665375*angs2bohr;
    Box.vec[1]= 12.431071665375*angs2bohr;
    Box.vec[2]= 12.431071665375*angs2bohr;
    beta=0.75;
    for (int dim=0;dim<NDIM;dim++)
      BoxInv.vec[dim]=1.0/Box.vec[dim];
      
    raggio=2.9655460304201733; 

    qO=-1.1995807E+00;
      //    qO=-1.1995809;
    qH=0.59979035;

    cOO=6.8676628;
    //    cOO=6.86766;
    //    cOH=-2.04454;
    cOH=-2.0445413;

    //    bOO=2.24850;
    bOO=2.2485039;
    
    //    bOH=3.92473;
    bOH=3.9247332;

    //    alpha_O=4.08675;
    alpha_O=4.0867573;
    dipfac=1.01977602741587826E-002;

    D1_OO=-1.26367e-4;
    D1_OH=3.77059e-5;
    D1_HH=8.71664e-1;

    D2_OO=3.08029e-4;
    D2_OH=8.67779e-6;
    D2_HH=1.93115e-7;

    

    gamma1_OO=13.09317;
    gamma1_OH=15.20544;
    gamma1_HH=13.13611;

    gamma2_OO=13.96521;
    gamma2_OH=12.38136;
    gamma2_HH=16.13997;

    r1_OO=7.473237;
    r1_OH=3.08451;
    r1_HH=0.38486;

    r2_OO=7.03818;
    r2_OH=5.63316;
    r2_HH=8.12431;



  }

};

class SystemClass
{
public:
  HelpClass h;
  vector<dVecp> r;
  vector<int> atom;
  vector<double> k_factor;
  vector<dVecp> kPoints;
  vector<double> q;
  
  int NumAtoms;
  int NumMolecules;
  ParamsClass params;

  void DistDisp (dVecp &r1, dVecp &r2,
	       double &dist, dVecp &disp)
{
  for (int dim=0;dim<NDIM;dim++)
    disp.vec[dim] =r1.vec[dim]-r2.vec[dim];
  dVecp n;
  n.vec[0] = nearbyint(disp.vec[0]*params.BoxInv.vec[0]);
  n.vec[1] = nearbyint(disp.vec[1]*params.BoxInv.vec[1]);
  n.vec[2] = nearbyint(disp.vec[2]*params.BoxInv.vec[2]);
  disp.vec[0] -= n.vec[0]*params.Box.vec[0];
  disp.vec[1] -= n.vec[1]*params.Box.vec[1];
  disp.vec[2] -= n.vec[2]*params.Box.vec[2];
  dist =sqrt(disp.vec[0]*disp.vec[0]+
	     disp.vec[1]*disp.vec[1]+disp.vec[2]*disp.vec[2]);
}
  
  
  void ReadPositions(int numMolecules)
  {
    cerr<<"Reading positions"<<endl;
    ifstream infile;
    infile.open("wat.pos");
    assert(infile);
    string garbage;
    infile>>garbage;
    r.resize(numMolecules*3);
    atom.resize(numMolecules*3);
    for (int i=0;i<r.size();i++){
      for (int dim=0;dim<NDIM;dim++){
	infile>>r[i].vec[dim];
	r[i].vec[dim]=r[i].vec[dim]*angs2bohr;
      }
      infile>>atom[i];
      infile>>atom[i];
      atom[i]=atom[i]-1;
      if (atom[i]==0)
	q.push_back(params.qO);
      else 
	q.push_back(params.qH);
    }
    infile.close();
    cerr<<"Done Reading positions"<<endl;
  }

  void ReadKPoints()
  {
  ifstream infile;
  infile.open("kPoints.txt");
  assert(infile);
  while (!infile.eof()){
    dVecp k;
    infile>>k.vec[0];
    if (!infile.eof()){
      
      infile>>k.vec[1];
      infile>>k.vec[2];
      double kf;
      infile>>kf;
      k_factor.push_back(kf);
      kPoints.push_back(k);
    }
  }
  //  cerr<<"The number of k points is "<<kPoints.size()<<endl;
  infile.close();
  }
  
  
  
  
};

class WaterEnergy{
public:
  SystemClass system;
  HelpClass h;
  void Init(int numMolecules)
  {
    system.ReadKPoints();
    system.params.Init();
    system.params.Init();
    system.NumMolecules=numMolecules;
    system.NumAtoms=3*system.NumMolecules;
    efield.resize(system.NumAtoms);
    efield_old.resize(system.NumAtoms);
    h.ZeroVecp(efield);
    h.ZeroVecp(efield_old);
  }

  vector<dVecp> efield;
  vector<dVecp> efield_old;




double converged(SystemClass &system,
	       vector<dVecp> &efield,
	       vector<dVecp> &efield_old)
{
  double diff=0.0;
  for (int i=0;i<efield.size();i++){
    if (system.atom[i]==0){
      dVecp sub;
      for (int dim=0;dim<NDIM;dim++)
	sub.vec[dim]=efield[i].vec[dim]-efield_old[i].vec[dim];
      diff=diff+h.dot(sub,sub)*system.params.alpha_O*system.params.alpha_O/efield.size();
    }

  }
  diff =sqrt(diff);
  //  cerr<<"My diff is "<<diff<<endl;
  //  return (diff<5e-4);
  return diff;
}


double ShortRangeEnergy(SystemClass &system)
{
  double energy=0.0;
  for (int i=0;i<system.r.size();i++){
    dVecp ii;
    for (int j=i+1;j<system.r.size();j++){
      //      if (i>j){
      //	for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	  {
	    //	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	    {
	      //	    for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	      {
	      dVecp r12;
	      double dist;
	      system.DistDisp (system.r[i], system.r[j],dist,r12);
	      //	      r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	      //	      r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	      //	      r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	      double dist2=h.dot(r12,r12);
	      dist=sqrt(dist2);
	      double D1=0.0; double D2=0.0; double gamma1=0.0; double gamma2=0.0;
	      double r1=0.0; double r2=0.0;
	      if (system.atom[i]==0 && system.atom[j]==0){
		D1=system.params.D1_OO;
		D2=system.params.D2_OO;
		gamma1=system.params.gamma1_OO;
		gamma2=system.params.gamma2_OO;
		r1=system.params.r1_OO;
		r2=system.params.r2_OO;
	      }
	      else if (system.atom[i]==1 && system.atom[j]==1){
		D1=system.params.D1_HH;
		D2=system.params.D2_HH;
		gamma1=system.params.gamma1_HH;
		gamma2=system.params.gamma2_HH;
		r1=system.params.r1_HH;
		r2=system.params.r2_HH;


	      }
	      else {
		D1=system.params.D1_OH;
		D2=system.params.D2_OH;
		gamma1=system.params.gamma1_OH;
		gamma2=system.params.gamma2_OH;
		r1=system.params.r1_OH;
		r2=system.params.r2_OH;


	      }
	      if (dist<CO){
		energy+=D1*(exp(gamma1*(1-(dist/r1)))-2*exp(gamma1/2.0*(1-dist/r1)))+
		  D2*(exp(gamma2*(1-dist/r2))-2*exp((gamma2/2.0)*(1-dist/r2)));
		if (dist<0.3)
		  energy+=99999999;
	      }
	    }
	  }
	  //	}
      }
    }
  }
  return energy;
}
		     
double ewald(SystemClass  &system,
	   vector<dVecp> &force)
{



  double eself=0.0;
  double recipEnergy=0.0;
  for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
    eself=eself-system.q[ptcl]*system.q[ptcl]/(sqrt(2*M_PI)*system.params.raggio);
  }
  

  vector<complex<double> > phi(system.NumAtoms);
  for (int ki=0;ki<system.kPoints.size();ki++){
    complex<double> sumPhi=0.0;
    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){


      double sdotprod=h.dot(system.kPoints[ki],system.r[ptcl]);
      complex<double> e_iks(cos(sdotprod),sin(sdotprod));
      phi[ptcl]=system.q[ptcl]*e_iks;
      sumPhi+=phi[ptcl];
    }
    recipEnergy+=2.0*system.k_factor[ki]*(sumPhi.real()*sumPhi.real()+sumPhi.imag()*sumPhi.imag());
    
    for (int dim=0;dim<NDIM;dim++){
      double twofg=2.0*system.k_factor[ki]*system.kPoints[ki].vec[dim];
      for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
	complex<double> dotprod=phi[ptcl]*conj(sumPhi);
	force[ptcl].vec[dim]+=2.0*dotprod.imag()*twofg; //2.0 for mysterious reasons
      }
    }
  }
  //reciprocal done. Now let's try real

  double realEnergy=0.0;
  for (int i=0;i<system.r.size();i++){
    dVecp ii;
    for (int j=0;j<system.r.size();j++){
      if (i!=j){
	//	for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	  {
	    //	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	    {
	      //	    for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	      {
	       dVecp r12;
	       double dist;
	       system.DistDisp (system.r[i], system.r[j],dist,r12);
	       //	       r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	       //	       r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	       //	       r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	       double dist2=h.dot(r12,r12);
	       dist=sqrt(dist2);
	       if (dist<CO){
		 double qi=system.q[i];
		 double qj=system.q[j];
		 double gamjir2=1.0/dist2;
		 double rckj=system.params.raggio*sqrt(2.0);
		 double jkcr=1.0/rckj;
		 double factor1=qi*qj*erfc(dist*jkcr)/dist;
		 double oiggar=1.0/system.params.raggio;
		 double factor2=qi*qj*exp(-dist2*jkcr*jkcr)*
		   gamjir2*oiggar;
		 for (int dim=0;dim<NDIM;dim++){
		   double dforce=(factor1*gamjir2+sqrt(2.0/M_PI)*factor2)*r12.vec[dim];
		   force[i].vec[dim]+=dforce;
		   
		 }
		 realEnergy+=factor1*0.5;
	       }
	    }
	  }
	}
      }
    }
  }
  //  cerr<<"EWALD ENERGIES: "<<realEnergy<<" "<<recipEnergy<<" "<<eself<<endl;
  return realEnergy+recipEnergy+eself;
}



void Dipole_sr(SystemClass &system, vector<dVecp> &dip_sr)
{
  h.ZeroVecp(dip_sr);

  for (int i=0;i<system.r.size();i++){
    if (system.atom[i]==0){

      double fij=0.0;
      double c=0.0;
      double b=0.0;
      double q=0.0;

      dVecp ii;
      for (int j=0;j<system.r.size();j++){
 	if (i!=j){
	  //	  for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	    {
	      //	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	      {
		//	  for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	    {
	  
	  double dist;
	  dVecp r12;
	  system.DistDisp (system.r[i], system.r[j],dist,r12);
	  //	  r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	  //	  r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	  //	  r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	  
	  double dist2=h.dot(r12,r12);
	  dist=sqrt(dist2);
	  if (dist<CO){
	    
	    q=system.q[j];
	    if (system.atom[i]==0 && system.atom[j]==0){
	      c=system.params.cOO;
	      b=system.params.bOO;
	    }
	    else {
	      c=system.params.cOH;
	      b=system.params.bOH;
	    }
	    double sum=0.0;
	    double k_factorial=1;
	    for (int k=0;k<=4;k++){
	      sum+=pow(b*dist,k)/k_factorial;
	      k_factorial=k_factorial*(k+1);
	    }
	    fij=c*exp(-b*dist)*sum;
	    //	    cerr<<"rijA: "<<system.params.alpha_O<<" "<<q<<" "<<" "<<fij<<" "<<dist<<" "<<r12.vec[0]<<" "<<r12.vec[1]<<" "<<r12.vec[2]<<endl;
	    for (int dim=0;dim<NDIM;dim++){
	      dip_sr[i].vec[dim]+=system.params.alpha_O*q*fij*
		(r12.vec[dim])/(dist*dist*dist);
		//		(r12.vec[dim]+ii.vec[dim]*system.params.Box.vec[dim])/(dist*dist*dist);
	    }
	    //	    cerr<<dip_sr[i].vec[0]<<" "<<dip_sr[i].vec[1]<<" "<<dip_sr[i].vec[2]<<endl;
	  }
	  }
	  }
	  }
 	}  

       }
     }
   }
}


double Dipole_sr_energy(SystemClass &system, vector<dVecp> &dip_sr)
{
  double energy=0.0;
  for (int i=0;i<system.r.size();i++){
    //    if (system.atom[i]==0){

      double fij=0.0;
      double c=0.0;
      double b=0.0;
      double q=0.0;

      dVecp ii;
      for (int j=0;j<system.r.size();j++){
 	if (i!=j){
	  //	  for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	  {
	    //	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	    {
	      //	  for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	    {
	  
	  double dist;
	  dVecp r12;
	  system.DistDisp (system.r[i], system.r[j],dist,r12);
	  //	  r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	  //	  r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	  //	  r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	  
	  double dist2=h.dot(r12,r12);
	  dist=sqrt(dist2);
	  double dist3=dist*dist2;
	  if (dist<CO){
	    
	    q=system.q[j];
	    if (system.atom[i]==0 && system.atom[j]==0){
	      c=system.params.cOO;
	      b=system.params.bOO;
	    }
	    else {
	      c=system.params.cOH;
	      b=system.params.bOH;
	    }
	    double sum=0.0;
	    double k_factorial=1;

	    for (int k=0;k<=4;k++){
	      sum+=pow(b*dist,k)/k_factorial;
	      k_factorial=k_factorial*(k+1);
	    }
	    fij=c*exp(-b*dist)*sum;
	    //	    dgijdrij*=exp(-b*dist);

	    
	    //	    cerr<<"dip_sr_i "<<dip_sr[i].vec[0]<<" "<<dip_sr[i].vec[1]<<" "<<dip_sr[i].vec[2]<<endl;
	    double pri=h.dot(dip_sr[i],r12);
	    double prj=h.dot(dip_sr[j],r12);
	    double denergy=(system.q[i]*prj-system.q[j]*pri)*fij*(1.0/(dist3));
	    //	    cerr<<"Denergy is "<<0.5*denergy<<" "<<prj<<" "<<pri<<" "<<fij<<" "<<(1.0/dist3)<<" "<<dist<<endl;

	    energy+=0.5*denergy;
	    //	    double const1=gij*(dist2/dist);
	    //	    double const2=(system.q[j]*pri-system.q[i]*prj)*
	    //	      (dist*dgijdrij-3.0*gij)*(dist3*dist2);
	    
	    

// 	    for (int dim=0;dim<NDIM;dim++){
// 	      dip_sr[i].vec[dim]+=system.params.alpha_O*q*fij*
// 		(r12.vec[dim])/(dist*dist*dist);
// 	    }


	  }
	  }
	  }
	  }
	  // 	}  

       }
     }
   }
  //  cerr<<"My dipole sr energy is "<<energy<<endl;
  return energy;
}


double ComputeEnergy()
{



  TimerClass Timer("AllTimer");


  double oiggar=1.0/system.params.raggio;
  double oiggar3=oiggar*oiggar*oiggar;
  double oiggar5=oiggar3*oiggar*oiggar;
  double rckj=system.params.raggio*sqrt(2.0);     
  double jkcr=1.0/rckj;
  double myConst=1.0/(3.0*(pow(system.params.raggio,3))*sqrt(2*M_PI));


  Timer.Start();
  vector<dVecp> force(system.NumAtoms);
  h.ZeroVecp(force);
  double ewald_energy=ewald(system,force);
  Timer.Stop();
  //  cerr<<"Ewald takes "<<Timer.Time()<<endl;
  Timer.Clear();
  //  cerr<<"EWALD FORCE: "<<endl;
  //  PrintVec(force);
  //  cerr<<"EWALD FORCE DONE: "<<endl;
  
  Timer.Start();
  vector<dVecp> dip_sr;
  dip_sr.resize(system.NumAtoms);
  Dipole_sr(system,dip_sr);
  Timer.Stop();
  //  cerr<<"Dipole SR takes "<<Timer.Time()<<endl;
  Timer.Clear();

  //  cerr<<"DIPOLE SR"<<endl;
  //  PrintVec(dip_sr);
  //  cerr<<"DIPOLE SR DONE"<<endl;


  vector<dVecp> eftot;
  vector<dVecp> dip;
  //  HACK efield.resize(system.NumAtoms);
  //  HACK efield_old.resize(system.NumAtoms);
  eftot.resize(system.NumAtoms);
  dip.resize(system.NumAtoms);



  //  h.ZeroVecp(eftot);
  //Maybe you want to set efield to some old stuff
  //HACK  h.ZeroVecp(efield);
  Timer.Start();
  double diff=1.0;
  double diff_old=0.0;
  int step=-1;
  //  for (int step=0;step<5;step++){
  do {
    step++;
    h.ZeroVecp(dip);

    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
      if (system.atom[ptcl]==0){
        for (int dim=0;dim<NDIM;dim++){
	  if (step==0)
	    eftot[ptcl].vec[dim]=efield[ptcl].vec[dim]
	      +force[ptcl].vec[dim]/system.params.qO;
	  else 
	    eftot[ptcl].vec[dim]=efield[ptcl].vec[dim]*system.params.beta + 
	      (1-system.params.beta)*efield_old[ptcl].vec[dim]+
	      force[ptcl].vec[dim]/system.params.qO;;
	  dip[ptcl].vec[dim]=eftot[ptcl].vec[dim]*system.params.alpha_O+dip_sr[ptcl].vec[dim];
	}
      }
    }
    
    //    cerr<<"DIPOLE aa"<<endl;
    //    PrintVec(dip);
    //    cerr<<"DIPOLE done aa"<<endl;


    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++)
      for (int dim=0;dim<NDIM;dim++)
	efield[ptcl].vec[dim]=dip[ptcl].vec[dim]*system.params.dipfac;
  
    //    cerr<<"Efield before recpolewald"<<endl;
    //    PrintVec(efield);
    //    cerr<<"Done"<<endl;

    

// // ///recpolewald_short.f
    vector<complex<double > > rhok_dipole;
   
    rhok_dipole.resize(system.kPoints.size());
    h.ZeroVecp(rhok_dipole);
	       
    for (int ptcl=0;ptcl<system.r.size();ptcl++){
      for (int ki=0;ki<system.kPoints.size();ki++){

        double kdotr=h.dot(system.kPoints[ki],system.r[ptcl]);
        double dipdotr=h.dot(system.kPoints[ki],dip[ptcl]);
	complex<double> e_ikr(cos(kdotr),sin(kdotr));
        rhok_dipole[ki]+=e_ikr*dipdotr;
      }
    }

    for (int ki=0;ki<system.kPoints.size();ki++){
      for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
	double kdotr=h.dot(system.kPoints[ki],system.r[ptcl]);
        complex<double> e_mikr(cos(kdotr),-1*sin(kdotr));
        for (int dim=0;dim<NDIM;dim++)
	  efield[ptcl].vec[dim]=efield[ptcl].vec[dim]-
	    4*system.kPoints[ki].vec[dim]*system.k_factor[ki]* (e_mikr*rhok_dipole[ki]).real();
      }
    }


   
// //    //realpolewald_short


    for (int i=0;i<system.r.size();i++){
      for (int j=0;j<system.r.size();j++){
      dVecp ii;
      if (i!=j){
	//	for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	  {
	    // 	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	    {
	      //	    for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	      {

	      dVecp r12;
	      double dist;
	      system.DistDisp (system.r[i], system.r[j],dist,r12);
	      //	      r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	      //	      r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	      //	      r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	      double dist2=h.dot(r12,r12);
	      dist=sqrt(dist2);
	      if (dist<CO){

 	     double rijmag=dist;
 	     double rijmag2=dist2;
 	     double gamjir2=1.0/rijmag2;
 	     double gamjir3=gamjir2/rijmag;

 	     double factor2=sqrt(2.0/M_PI)*exp(-dist2*jkcr*jkcr)*gamjir2;
 	     double factor1=gamjir3*erfc(rijmag*jkcr)+oiggar*factor2;
 	     double factor3=3.0*factor1*gamjir2+factor2*oiggar3;
	     
 	     double prj=h.dot(dip[j],r12);
 	     for (int dim=0;dim<NDIM;dim++)
	       efield[i].vec[dim]=efield[i].vec[dim]-dip[j].vec[dim]*factor1+prj*(r12.vec[dim])*factor3;
	      }
	    }
 	  }
	}
      }
      }
    }
    //  cerr<<"This is the efield"<<endl;
    //  PrintVec(efield);
    //  cerr<<"This is the efield done"<<endl;
    diff_old=diff;
    diff=converged(system,efield,efield_old);
    //    cerr<<abs(diff_old-diff)<<endl;
    //    cerr<<"I am converged: "<<(abs(diff_old-diff)<5e-4)<<" "<<(1==1)<<" "<<abs(diff_old-diff)<<" "<<diff<<" "<<diff_old<<endl;
    //    cerr<<"I am converged: "<<" "<<abs(diff_old-diff)<<" "<<diff<<" "<<diff_old<<endl;

    for (int i=0;i<efield.size();i++){
      for (int dim=0;dim<NDIM;dim++)
	efield_old[i].vec[dim]=efield[i].vec[dim];
    }


  } while (abs(diff-diff_old)>5e-4);
  Timer.Stop();
  cerr<<"iterations take "<<Timer.Time()<<" "<<step<<" "<<diff<<" "<<diff_old<<endl;
  Timer.Clear();

//    ///now compute the energy
  Timer.Start();
  double e_ind=0.0;
  double e_dd=0.0;
    double e_dd_tmp=0.0;
    double t_e_dd=0.0;
  double e_qd=0.0;
  double t_e_qd=0.0;
  for (int i=0;i<system.r.size();i++){
    for (int j=0;j<system.r.size();j++){
      dVecp ii;
        if (i!=j){
	  // 	 for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++)
	   {
	     // 	   for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++)
	     {
	       // 	     for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++)
	       {

 	       dVecp r12;
 	       double dist;
 	       system.DistDisp (system.r[i], system.r[j],dist,r12);
	       // 	       r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0];
	       // 	       r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1];
	       // 	       r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2];
	       double dist2=h.dot(r12,r12);
	       dist=sqrt(dist2);
 	       if (dist<CO){
 		 double prj=0.0;
		 double pri=0.0;
 		 double pp=0.0;

 		 double rijmag2 = dist2;
 		 double rijmag = dist;
 		 double gamjir2=1.0/rijmag2;
 		 double gamjir3=gamjir2/rijmag;
 		 double factor2=sqrt(2.0/M_PI)*exp(-dist2*jkcr*jkcr)*gamjir2;
 		 double factor1=gamjir3*erfc(rijmag*jkcr)+oiggar*factor2;
 		 double factor3=3.0*factor1*gamjir2+factor2*oiggar3;
		 
 		 if (system.atom[i]==0){
		   pri=h.dot(dip[i],r12);
 		 }
 		 if (system.atom[j]==0){
		   prj=h.dot(dip[j],r12);
 		   e_qd+=system.q[i]*prj*factor1;		   
 		 }
 		 if (system.atom[i]==0 && system.atom[j]==0){
		   pp=h.dot(dip[i],dip[j]);
 		   e_dd+=0.5*(pp*factor1-pri*prj*factor3);
		 }
 	       }
 	     }
 	   }
 	 }
        }
    }
    
    double pp=h.dot(dip[i],dip[i]);
    
    e_dd_tmp=e_dd_tmp-pp*myConst;
  
    //    cerr<<"post: "<<pp<<" "<<myConst<<" "<<e_dd_tmp<<endl;
    if (system.atom[i]==0)
      e_ind+=0.5*pp/system.params.alpha_O;

  }
  e_dd=e_dd+e_dd_tmp;
  
  vector<complex<double > > rhok_dipole;
  
  rhok_dipole.resize(system.kPoints.size());
  
  

  h.ZeroVecp(rhok_dipole);    
  for (int ki=0;ki<system.kPoints.size();ki++){
    complex<double> qjeSum=0.0;
    for (int ptcl=0;ptcl<system.r.size();ptcl++){
      
      
      double kdotr=h.dot(system.kPoints[ki],system.r[ptcl]);
      double dipdotr=h.dot(system.kPoints[ki],dip[ptcl]);
      complex<double> e_ikr(cos(kdotr),sin(kdotr));
      rhok_dipole[ki]+=e_ikr*dipdotr;
      qjeSum+=system.q[ptcl]*conj(e_ikr);

      //	cerr<<"the dip is "<<rhok_dipole[ki]<<" "<<dip[ptcl].vec[0]<<" "<<dip[ptcl].vec[1]<<" "<<dip[ptcl].vec[2]<<endl;
    }
    t_e_qd+=-4.0*system.k_factor[ki]*  (qjeSum*rhok_dipole[ki]).imag();
    t_e_dd=t_e_dd+2.0*system.k_factor[ki]*
      (rhok_dipole[ki].real()*rhok_dipole[ki].real()+rhok_dipole[ki].imag()*rhok_dipole[ki].imag());
    //    cerr<<"pkjesum: "<<rhok_dipole[ki]<<" "<<t_e_dd<<endl;
  }
  //  cerr<<"It is "<<t_e_qd<<endl;
  e_dd+=t_e_dd;
  e_qd+=t_e_qd;
  Timer.Stop();
  //  cerr<<"energy takes "<<Timer.Time()<<endl;
  Timer.Clear();

  double dipole_sr_energy=Dipole_sr_energy(system, dip);
  //  cerr<<dipole_sr_energy<<endl;
  //  cerr<<e_dd<<" "<<e_qd<<" "<<e_ind<<" "<<endl;  
  Timer.Start();
  double sr_energy=ShortRangeEnergy(system);
  Timer.Stop();
  //  cerr<<"Short range time take "<<Timer.Time();
  Timer.Clear();
  //  cerr<<"short range energy: "<<sr_energy<<endl;
  //  cerr<<setprecision(9);
  double total=sr_energy+ewald_energy+e_dd+e_qd+e_ind+dipole_sr_energy;
  
  //  cerr<<"TOTAL: "<<total<<endl;
  return total;

}
};


/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class WaterClass : public ActionBaseClass
{
protected:
  int TotalTime;
  /// These are the coefficients used for the low-variance estimator
  /// for the gradient
public:
  WaterEnergy we;

  void Read (IOSectionClass &in)
  {
    we.Init(64);
    we.system.ReadPositions(64);

  }

  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
/*   { */
/*     //need to actually do something here to put the stuff from the path into the system */
/*     double en=we.ComputeEnergy(); */
/*     cerr<<en; */
/*     return PathData.Path.tau*en; */
/*     //    cerr<<we.ComputeEnergy()<<endl; */

/*   } */
   double d_dBeta (int slice1, int slice2, int level) 
  {
    double en=we.ComputeEnergy();
    return en;
  }

  string GetName()
  {
    return "water";
  }
 WaterClass(PathDataClass &pathData) :   ActionBaseClass (pathData)
  {
  }
	     
};

#endif
