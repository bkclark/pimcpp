#ifndef DVECP_H
#define DVECP_H

#include <vector>
#include "../Blitz.h"
using namespace std;

class dVecp
{
public:
  int size()
  {
    return vec.size();
  }

  vector<double> vec;
  dVecp()
  {
    vec.resize(NDIM);
    for (int i=0;i<NDIM;i++){
       vec[i]=0.0;
    }
    
  }

  friend dVecp operator*(double const& d, dVecp const& v)
    {
      dVecp out;
      for (int i=0;i<out.size();i++){
	out.vec[i]=d*v.vec[i];
      }
      return out;
}

  friend dVecp operator+(dVecp const& v1, dVecp const& v2)
    {
      dVecp out;
      for (int i=0;i<out.size();i++){
	out.vec[i]=v1.vec[i]+v2.vec[i];
      }
      return out;
}

  friend std::ostream& operator<< (std::ostream& output, dVecp& v)
    {
      for (int i=0;i<v.size();i++){
	output<<v.vec[i]<<" ";
      }
      return output;
    }
};



void ZeroVec(vector<dVecp> &myVec);


void ZeroVec(vector<double> &myVec);



void ZeroVec(vector<complex<double> > &myVec);


void PrintVec(vector<dVecp> &myVec);

double dot(dVecp &a, dVecp &b);


#endif

