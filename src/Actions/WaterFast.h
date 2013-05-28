
#ifndef WATER_H
#define WATER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;
#include <fstream>
#define NDIM 3
#include <iomanip>
#include "Timer.h"
#define angs2bohr 1.889725989
#include "dVec.h"
#include "SystemClass.h"
#include "ActionBase.h"
#include "Timer.h"


// oxygen -> 0 
// hydrogen -> 1

class WaterEnergy{
public:

 double converged(SystemClass &system,
	       vector<dVecp> &efield,
		 vector<dVecp> &efield_old);

 double ShortRangeEnergy(SystemClass &system,DistanceClass &d);

 double ShortRangeEnergy_slow(SystemClass &system,DistanceClass &d);

 void Dipole_sr(SystemClass &system, vector<dVecp> &dip_sr,DistanceClass &d);

 double Dipole_sr_energy(SystemClass &system, vector<dVecp> &dip_sr,DistanceClass &d);




  bool Fast;
  SystemClass system;
  DistanceClass d;
  RhoClass rho;
  double eself;
  vector<complex<double> > phi;
  vector<complex<double > > rhok_dipole;
  
  vector<dVecp> disp;
  vector<double> dist;
  vector<dVecp> dip_sr;
  vector<dVecp> force;
  vector<dVecp> efield;
  vector<dVecp> efield_old;
  vector<dVecp> eftot;
  vector<dVecp> dip;
  

  




  void Init(int numMolecules)
  {
    Fast=true;
    TimerClass Timer("Init");
    Timer.Start();
    system.ReadKPoints();
    system.params.Init();
    //    system.params.Init();
    system.NumMolecules=numMolecules;
    system.NumAtoms=3*system.NumMolecules;
    system.ReadPositions(numMolecules);
    eself=0.0;
    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){ 
      eself=eself-system.q[ptcl]*system.q[ptcl]/(sqrt(2*M_PI)*system.params.raggio); 
    } 
    phi.resize(system.NumAtoms);


    rhok_dipole.resize(system.kPoints.size());
    d.Init(3*numMolecules);
    rho.Init(system);

    //    d.ComputeDistDisp(system);
    //    rho.Compute(system);

    force.resize(system.NumAtoms);
    dip_sr.resize(system.NumAtoms);
    efield.resize(system.NumAtoms);
    efield_old.resize(system.NumAtoms);
    eftot.resize(system.NumAtoms);
    dip.resize(system.NumAtoms);
    Timer.Stop();
    //    cerr<<"Init Takes "<<Timer.Time()<<endl;
  }


double ewald(SystemClass  &system,
	   vector<dVecp> &force)
{
  
  double recipEnergy=0.0;
  for (int ki=0;ki<system.kPoints.size();ki++){
    complex<double> sumPhi=0.0;
    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
      //      double sdotprod=dot(system.kPoints[ki],system.r[ptcl]);
      //      complex<double> e_iks(cos(sdotprod),sin(sdotprod));
      complex<double> e_iks;
      e_iks=rho.Get_eikr(ptcl,ki);
      phi[ptcl]=system.q[ptcl]*e_iks;
      sumPhi+=phi[ptcl];
    }
    recipEnergy+=2.0*system.k_factor[ki]*
      (sumPhi.real()*sumPhi.real()+sumPhi.imag()*sumPhi.imag());

      for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
	complex<double> dotprod=phi[ptcl]*conj(sumPhi);
	for (int dim=0;dim<NDIM;dim++){
	  double twofg=2.0*system.k_factor[ki]*system.kPoints[ki].vec[dim];
	  force[ptcl].vec[dim]+=2.0*dotprod.imag()*twofg; //2.0 for mysterious reasons
	}
      }
  }
  //reciprocal done. Now let's try real


  double rckj=system.params.raggio*sqrt(2.0);
  double jkcr=1.0/rckj;
  double oiggar=1.0/system.params.raggio;

  double realEnergy=0.0;
  for (int i=0;i<system.r.size();i++){
    double qi=system.q[i];
    for (int j=i+1;j<system.r.size();j++){
      dVecp r12;
      double dist;
      double dist2;
      d.GetDistDisp(i,j,dist,dist2,r12);
      double qj=system.q[j];
      double gamjir2=1.0/dist2;
      double factor1=qi*qj*erfc(dist*jkcr)/dist;
      double factor2=qi*qj*exp(-dist2*jkcr*jkcr)*
	gamjir2*oiggar;
      for (int dim=0;dim<NDIM;dim++){
	double dforce=(factor1*gamjir2+sqrt(2.0/M_PI)*factor2)*r12.vec[dim];
	force[i].vec[dim]+=dforce;
	force[j].vec[dim]-=dforce;
      }
      realEnergy+=factor1; // *0.5;
    }
  }
  //A  cerr<<"EWALD ENERGIES: "<<realEnergy<<" "<<recipEnergy<<" "<<eself<<endl;
  return realEnergy+recipEnergy+eself;
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
  ZeroVec(force);
  double ewald_energy=ewald(system,force);
  Timer.Stop();
  //  cerr<<"Ewald takes "<<Timer.Time()<<endl;
  Timer.Clear();


  
  Timer.Start();
  Dipole_sr(system,dip_sr,d);
  Timer.Stop();
  //  cerr<<"Dipole SR takes "<<Timer.Time()<<endl;
  Timer.Clear();
  //  cerr<<"DIPOLE SR"<<endl;
  //  PrintVec(dip_sr);
  //  cerr<<"DIPOLE SR DONE"<<endl;


  Timer.Start();
  ZeroVec(efield);
  double diff=1.0;
  double diff_old=0.0;
  int step=-1;
  int tries=0;
  do {
    step++;
    ZeroVec(dip);
    
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
    //    PrintVec(dip);
    for (int ptcl=0;ptcl<system.NumAtoms;ptcl++)
      for (int dim=0;dim<NDIM;dim++)
	efield[ptcl].vec[dim]=dip[ptcl].vec[dim]*system.params.dipfac;
    //    PrintVec(efield);


// // ///recpolewald_short.f
    ZeroVec(rhok_dipole);
	       
    for (int ptcl=0;ptcl<system.r.size();ptcl++){
      for (int ki=0;ki<system.kPoints.size();ki++){


        double dipdotr=dot(system.kPoints[ki],dip[ptcl]);

	//        double kdotr=dot(system.kPoints[ki],system.r[ptcl]);
	//	complex<double> e_ikr(cos(kdotr),sin(kdotr));
	//	double kdotr;
	complex<double> e_ikr=rho.Get_eikr(ptcl,ki);
	//	rho.Get(ptcl,ki,kdotr,e_ikr);
        rhok_dipole[ki]+=e_ikr*dipdotr;
      }
    }

    for (int ki=0;ki<system.kPoints.size();ki++){
      for (int ptcl=0;ptcl<system.NumAtoms;ptcl++){
	//	double kdotr=dot(system.kPoints[ki],system.r[ptcl]);
	//        complex<double> e_mikr(cos(kdotr),-1*sin(kdotr));
	complex<double> e_mikr=conj(rho.Get_eikr(ptcl,ki));
        for (int dim=0;dim<NDIM;dim++)
	  efield[ptcl].vec[dim]=efield[ptcl].vec[dim]-
	    4*system.kPoints[ki].vec[dim]*
	    system.k_factor[ki]* (e_mikr*rhok_dipole[ki]).real();
      }
    }


   
// //    //realpolewald_short
    for (int i=0;i<system.r.size();i++){
      for (int j=0;j<system.r.size();j++){
      if (i!=j){
	      dVecp r12;
	      double dist;
	      double dist2;
	      d.GetDistDisp(i,j,dist,dist2,r12);

 	     double gamjir2=1.0/dist2;
 	     double gamjir3=gamjir2/dist;

 	     double factor2=sqrt(2.0/M_PI)*exp(-dist2*jkcr*jkcr)*gamjir2;
 	     double factor1=gamjir3*erfc(dist*jkcr)+oiggar*factor2;
 	     double factor3=3.0*factor1*gamjir2+factor2*oiggar3;
	     
 	     double prj=dot(dip[j],r12);
 	     for (int dim=0;dim<NDIM;dim++)
	       efield[i].vec[dim]=efield[i].vec[dim]-dip[j].vec[dim]*factor1+prj*(r12.vec[dim])*factor3;
	      }
	    }
    }
    //  PrintVec(efield);
    diff_old=diff;
    diff=converged(system,efield,efield_old);

    for (int i=0;i<efield.size();i++){
      for (int dim=0;dim<NDIM;dim++)
	efield_old[i].vec[dim]=efield[i].vec[dim];
    }

    tries++;
  } while (abs(diff-diff_old)>5e-4 && tries<10);
  Timer.Stop();
  //  cerr<<"iterations take "<<Timer.Time()<<" "<<tries<<endl;
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
      if (i!=j){
	dVecp r12;
	double dist;
	double dist2;
	d.GetDistDisp(i,j,dist,dist2,r12);
	//	system.DistDisp (system.r[i], system.r[j],dist,r12);
	//	double dist2=dot(r12,r12);
	
	double prj=0.0;
	double pri=0.0;
	double pp=0.0;
	
	//	double rijmag2 = dist2;
	//	double rijmag = dist;
	double gamjir2=1.0/dist2;
	double gamjir3=gamjir2/dist;
	double factor2=sqrt(2.0/M_PI)*exp(-dist2*jkcr*jkcr)*gamjir2;
	double factor1=gamjir3*erfc(dist*jkcr)+oiggar*factor2;
	double factor3=3.0*factor1*gamjir2+factor2*oiggar3;
	
	if (system.atom[i]==0){
	  pri=dot(dip[i],r12);
	}
	if (system.atom[j]==0){
	  prj=dot(dip[j],r12);
	  e_qd+=system.q[i]*prj*factor1;		   
	}
	if (system.atom[i]==0 && system.atom[j]==0){
	  pp=dot(dip[i],dip[j]);
	  e_dd+=0.5*(pp*factor1-pri*prj*factor3);
	}
      }
    }
  
    double pp=dot(dip[i],dip[i]);
    e_dd_tmp=e_dd_tmp-pp*myConst;

    if (system.atom[i]==0)
      e_ind+=0.5*pp/system.params.alpha_O;
  }
  e_dd=e_dd+e_dd_tmp;
  
  
  ZeroVec(rhok_dipole);    
  for (int ki=0;ki<system.kPoints.size();ki++){
    complex<double> qjeSum=0.0;
    for (int ptcl=0;ptcl<system.r.size();ptcl++){
      

      double dipdotr=dot(system.kPoints[ki],dip[ptcl]);
      complex<double> e_ikr=rho.Get_eikr(ptcl,ki);
      //      double kdotr=dot(system.kPoints[ki],system.r[ptcl]);
      //      complex<double> e_ikr(cos(kdotr),sin(kdotr));
      rhok_dipole[ki]+=e_ikr*dipdotr;
      qjeSum+=system.q[ptcl]*conj(e_ikr);

    }
    t_e_qd+=-4.0*system.k_factor[ki]*  (qjeSum*rhok_dipole[ki]).imag();
    t_e_dd=t_e_dd+2.0*system.k_factor[ki]*
      (rhok_dipole[ki].real()*rhok_dipole[ki].real()+rhok_dipole[ki].imag()*rhok_dipole[ki].imag());

  }

  e_dd+=t_e_dd;
  e_qd+=t_e_qd;
  Timer.Stop();
  //  cerr<<"energy takes "<<Timer.Time()<<endl;
  Timer.Clear();

  Timer.Start();
  double dipole_sr_energy=Dipole_sr_energy(system, dip,d);
  double sr_energy;
  if (Fast)
    sr_energy=ShortRangeEnergy(system,d); 
  else 
    sr_energy=ShortRangeEnergy_slow(system,d); 
  Timer.Stop();
  //  cerr<<"Short range time take "<<Timer.Time();
  Timer.Clear();
  double total=sr_energy+ewald_energy+e_dd+e_qd+e_ind+dipole_sr_energy;
  if (tries==10)
    total=999;
  //  cerr<<"TOTAL: "<<total<<endl;
  return total;

}
};


class WaterClass :public ActionBaseClass
{
  TimerClass t;
 public:
 WaterClass(PathDataClass &pathData) :   ActionBaseClass (pathData),  t("FullTime")
  {

  }
  
  void SetBox(dVecp &Box)
  {
    we.system.params.Box=Box;
    for (int dim=0;dim<NDIM;dim++)
      we.system.params.BoxInv.vec[dim]=1.0/Box.vec[dim];
  }
  void SetBox(double L)
  {
    dVecp Box;
    for (int dim=0;dim<Box.size();dim++)
      Box.vec[dim]=L;

    we.system.params.Box=Box;
    for (int dim=0;dim<NDIM;dim++)
      we.system.params.BoxInv.vec[dim]=1.0/Box.vec[dim];
  }

  WaterEnergy we;
  void Init()
  {

    we.Init(54);
    //    we.system.ReadPositions(54);

  }
  double SingleAction(int, int, const blitz::Array<int, 1>&, int);
  void Read (IOSectionClass &in);

  double d_dBeta(int, int, int)
  {

  } 
  string GetName()
  {
    return "PolarizedWater";
  }



};
#endif
