#include "../Blitz.h"
#include <fstream>
#include <vector>

class HistogramClass
{
 public:
  vector<double> histogram;
  int NumPoints;
  double totalPoints;
  double delta;
  double startVal;
  double endVal;
  void Init(int t_NumPoints, double t_startVal, double t_endVal)
    {
      NumPoints=t_NumPoints;
      startVal=t_startVal;
      endVal=t_endVal;
      delta=(endVal-startVal)/NumPoints;
      histogram.resize(NumPoints);
      totalPoints=0.0;
      for (int i=0;i<histogram.size();i++) histogram[i]=0.0;
      //    histogram=0.0;
    }
  void Clear()
  {
    for (int i=0;i<histogram.size();i++) histogram[i]=0.0;
    totalPoints=0.0;
  }
  HistogramClass()
    {

    }
  HistogramClass(int t_NumPoints, double t_startVal, double t_endVal) : 
    NumPoints(t_NumPoints), startVal(t_startVal), endVal(t_endVal)
    {
      delta=(endVal-startVal)/NumPoints;
      histogram.resize(NumPoints);
      totalPoints=0.0;
    }
  
  void add(double loc,double val)
    {
      int n_loc=(int)(trunc((loc-startVal)/delta));
      if (n_loc>=0 && n_loc<histogram.size())
	histogram[n_loc]+=val;
      totalPoints++;
    }
  void print()
    {
      for (int i=0;i<NumPoints;i++){
	cerr<<startVal+i*delta<<" "<<histogram[i]/totalPoints<<endl;
      }
    }

  void print(ofstream &outFile)
    {
      for (int i=0;i<NumPoints;i++){
	outFile<<startVal+i*delta<<" "<<histogram[i]/totalPoints<<endl;
      }
    }

};

