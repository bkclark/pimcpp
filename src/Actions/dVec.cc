#include "dVec.h"

void ZeroVec(vector<dVecp> &myVec)
{
  for (int i=0;i<myVec.size();i++){
    for (int dim=0;dim<NDIM;dim++){
      myVec[i].vec[dim]=0.0;
    }

  }

}

void ZeroVec(vector<double> &myVec)
{
  for (int i=0;i<myVec.size();i++){
      myVec[i]=0.0;
  }

}



void ZeroVec(vector<complex<double> > &myVec)
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
