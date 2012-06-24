#include "WaterFast.h"

double WaterEnergy::converged(SystemClass &system,
	       vector<dVecp> &efield,
	       vector<dVecp> &efield_old)
{
  double diff=0.0;
  for (int i=0;i<efield.size();i++){
    if (system.atom[i]==0){
      dVecp sub;
      for (int dim=0;dim<NDIM;dim++)
	sub.vec[dim]=efield[i].vec[dim]-efield_old[i].vec[dim];
      diff=diff+dot(sub,sub)*system.params.alpha_O*system.params.alpha_O/efield.size();
    }

  }
  diff =sqrt(diff);
  //  cerr<<"My diff is "<<diff<<endl;
  //  return (diff<5e-4);
  return diff;
}


double WaterEnergy::ShortRangeEnergy(SystemClass &system,DistanceClass &d)
{

  double energy=0.0;
  for (int i=0;i<system.r.size();i++){
    dVecp ii;
    for (int j=i+1;j<system.r.size();j++){
      //if (i>j){
      //	for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++){
      //	  for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++){
      //	    for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++){
	      dVecp r12;
	      double dist;
	      double dist2;
	      d.GetDistDisp(i,j,dist,dist2,r12);

/* 	      system.DistDisp (system.r[i], system.r[j],dist,r12); */
/* 	      //	      r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0]; */
/* 	      //	      r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1]; */
/* 	      //	      r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2]; */
/* 	      dist2=dot(r12,r12); */
/* 	      dist=sqrt(dist2); */


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
	      //	      if (dist<18.0){
		energy+=D1*(exp(gamma1*(1-(dist/r1)))-2*exp(gamma1/2.0*(1-dist/r1)))+
		  D2*(exp(gamma2*(1-dist/r2))-2*exp((gamma2/2.0)*(1-dist/r2)));
		//	    }
		//	  }
		//	}
	  //	}
		//      }
    }
  }
  return energy;
}
		     

double WaterEnergy::ShortRangeEnergy_slow(SystemClass &system,DistanceClass &d)
{

  double energy=0.0;
  for (int i=0;i<system.r.size();i++){
    dVecp ii;
    for (int j=i+1;j<system.r.size();j++){
      //if (i>j){
      for (ii.vec[0]=-1;ii.vec[0]<=1;ii.vec[0]++){
	for (ii.vec[1]=-1;ii.vec[1]<=1;ii.vec[1]++){
	  for (ii.vec[2]=-1;ii.vec[2]<=1;ii.vec[2]++){
	      dVecp r12;
	      double dist;
	      double dist2;
	      //	      d.GetDistDisp(i,j,dist,dist2,r12);

 	      system.DistDisp (system.r[i], system.r[j],dist,r12); 
	      r12.vec[0]=r12.vec[0]+ii.vec[0]*system.params.Box.vec[0]; 
	      r12.vec[1]=r12.vec[1]+ii.vec[1]*system.params.Box.vec[1]; 
	      r12.vec[2]=r12.vec[2]+ii.vec[2]*system.params.Box.vec[2]; 
 	      dist2=dot(r12,r12); 
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
	      if (dist<18.0){
		energy+=D1*(exp(gamma1*(1-(dist/r1)))-2*exp(gamma1/2.0*(1-dist/r1)))+
		  D2*(exp(gamma2*(1-dist/r2))-2*exp((gamma2/2.0)*(1-dist/r2)));
	      }
	  }
	}
      }
    }
  }
  return energy;
}
		     



void WaterEnergy::Dipole_sr(SystemClass &system, vector<dVecp> &dip_sr,DistanceClass &d)
{
  ZeroVec(dip_sr);

  for (int i=0;i<system.r.size();i++){
    if (system.atom[i]==0){
      double fij=0.0;
      double c=0.0;
      double b=0.0;
      double q=0.0;
      for (int j=0;j<system.r.size();j++){
 	if (i!=j){
	  double dist;
	  dVecp r12;
	  //	  system.DistDisp (system.r[i], system.r[j],dist,r12);
	  double dist2; //=dot(r12,r12);
	  d.GetDistDisp(i,j,dist,dist2,r12);
	  q=system.q[j];
	  if (system.atom[j]==0){
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
	  for (int dim=0;dim<NDIM;dim++){
	    dip_sr[i].vec[dim]+=system.params.alpha_O*q*fij*
	      (r12.vec[dim])/(dist*dist*dist);
	  }
	  
 	}  
	
      }
    }
  }
}


double WaterEnergy::Dipole_sr_energy(SystemClass &system, vector<dVecp> &dip_sr,DistanceClass &d)
{
  double energy=0.0;
  for (int i=0;i<system.r.size();i++){
      double fij=0.0;
      double c=0.0;
      double b=0.0;
      double q=0.0;
      for (int j=0;j<system.r.size();j++){
 	if (i!=j){
	  double dist;
	  dVecp r12;
	  double dist2;
	  d.GetDistDisp(i,j,dist,dist2,r12);
	  //	  system.DistDisp (system.r[i], system.r[j],dist,r12);
	  //	  double dist2=dot(r12,r12);
	  double dist3=dist*dist2;
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
	    
	    double pri=dot(dip_sr[i],r12);
	    double prj=dot(dip_sr[j],r12);
	    double denergy=(system.q[i]*prj-system.q[j]*pri)*fij*(1.0/(dist3));
	    energy+=0.5*denergy;
	  }
      }
   }
  //  cerr<<"My dipole sr energy is "<<energy<<endl;
  return energy;
}

