#ifndef SYSTEM_CLASS_H
#define SYSTEM_CLASS_H

#include "ParamsClass.h"
#include "dVec.h"

class SystemClass
{
public:
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

  void DistDisp (dVecp &r1, dVecp &r2,
		 double &dist, double &dist2,dVecp &disp)
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
  dist2 =disp.vec[0]*disp.vec[0]+
	     disp.vec[1]*disp.vec[1]+disp.vec[2]*disp.vec[2];
  dist=sqrt(dist2);
}


  void ReadPositions_old(int numMolecules)
  {
    ifstream infile;
    
    infile.open("oxygenPositions2");
    string garbage;
    r.resize(numMolecules);
    atom.resize(numMolecules);
    for (int i=0;i<r.size();i++){
      for (int dim=0;dim<NDIM;dim++){
	infile>>r[i].vec[dim];
	//	r[i].vec[dim]=r[i].vec[dim]*angs2bohr;
      }
      atom[i]=0;
      if (atom[i]==0)
	q.push_back(params.qO);
      else 
	q.push_back(params.qH);
    }
    infile.close();
  }
  void ScalePositions(double s)
  {
    for (int i=0;i<r.size();i++){
      for (int dim=0;dim<NDIM;dim++){
	r[i].vec[dim]*=s;
      }

    }
  }


  void ReadPositions(int numMolecules)
  {
    ifstream infile;
    infile.open("wat.pos");

    string garbage;
    infile>>garbage;
    r.resize(numMolecules*3);
    atom.resize(numMolecules*3);
    for (int i=0;i<r.size();i++){
      for (int dim=0;dim<NDIM;dim++){
	infile>>r[i].vec[dim];
	//	r[i].vec[dim]*=0.75*0.8; //0.9;
	//	r[i].vec[dim]=r[i].vec[dim]*angs2bohr;
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
  }

  void ReadKPoints()
  {
  ifstream infile;
  infile.open("kPoints.txt");
  while (!infile.eof()){
    dVecp k;
    infile>>k.vec[0];
    k.vec[0]=k.vec[0]; // /0.75;
    if (!infile.eof()){
      
      infile>>k.vec[1];
      infile>>k.vec[2];
      k.vec[1]=k.vec[1]; // /0.75;
      k.vec[2]=k.vec[2]; // /0.75;

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



class DistanceClass
{
 public:
  vector<double> dist2;
  vector<double> dist;
  vector<dVecp> disp;
  int NumParticles;
  void Init(int t_NumParticles)
  {
    NumParticles=t_NumParticles;
    dist.resize(NumParticles*NumParticles);
    dist2.resize(NumParticles*NumParticles);
    disp.resize(NumParticles*NumParticles);
  }
  void ComputeDistDisp(SystemClass &system)
  {
    int pos=-1;
    for (int i=0;i<NumParticles;i++){
      for (int j=0;j<NumParticles;j++){
	double pos=i*NumParticles+j;
	system.DistDisp(system.r[i],system.r[j],dist[pos],dist2[pos],disp[pos]);
      }
    }
  }
  inline double GetDist(int i, int j)
  {
    return dist[i*NumParticles+j];
  }

  inline dVecp  GetDisp(int i, int j)
  {
    return disp[i*NumParticles+j];
  }
  inline void GetDistDisp(int i,int j,double &t_dist,double &t_dist2,dVecp &t_disp)
  {
    int pos=i*NumParticles+j;
    t_dist=dist[pos];
    t_disp=disp[pos];
    t_dist2=dist2[pos];
  }
  
};

class RhoClass
{
 public:
  vector<double> kdotr;
  vector<complex<double> > eikdotr;
  int Num_kpoints;
  int NumParticles;
  void Init(SystemClass &system)
  {
    Num_kpoints=system.kPoints.size();
    NumParticles=system.NumAtoms;
    kdotr.resize(system.kPoints.size()*system.NumAtoms);
    eikdotr.resize(system.kPoints.size()*system.NumAtoms);
  }
  void Compute(SystemClass &system)
  {
    for (int ptcl=0;ptcl<NumParticles;ptcl++){
      for (int ki=0;ki<Num_kpoints;ki++){
	int pos=ptcl*Num_kpoints+ki;
	kdotr[pos]=dot(system.kPoints[ki],system.r[ptcl]);
	complex<double> e_ikr(cos(kdotr[pos]),sin(kdotr[pos]));
	eikdotr[pos]=e_ikr;
      }
    }
    
  }
  complex<double> Get_eikr(int ptcl, int ki)
  {
    int pos=ptcl*Num_kpoints+ki;
    return eikdotr[pos];
  }
  double Get_kdotr(int ptcl, int ki)
  {
    int pos=ptcl*Num_kpoints+ki;
    return kdotr[pos];
  }

  void Get(int ptcl, int ki,double &t_kdotr, complex<double> &t_eikr)
  {
    int pos=ptcl*Num_kpoints+ki;
    t_kdotr=kdotr[pos];
    t_eikr=eikdotr[pos];
  }

};
#endif
