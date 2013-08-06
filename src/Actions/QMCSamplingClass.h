#ifndef QMC_SAMPLING_CLASS_H
#define QMC_SAMPLING_CLASS_H

//#include <QMCApp/QMCInterface.h>
#include "ActionBase.h"
#include <vector>
#include <fstream>

#ifdef USE_QMC
class CEIMCActionClass: public ActionBaseClass
{
	ofstream out;
  std::vector<double>* qmcData;
  int EnergyIndex0, EnergyIndex1, EnergyDiffIndex;
	Array<string, 1> ptclSet0, ptclSet1;
  
public:
	double dt;
	int walkers, chains, steps, blocks;
	bool correlated;
	string QMCMethod;
	double SingleAction(int slice1,int slice2,
		      const Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
  void Read (IOSectionClass &in);
  
  CEIMCActionClass(PathDataClass &pathData);
};
#endif

class IonIonActionClass: public ActionBaseClass
{
  double prefactor;
  double cutoff;
public:
  double SingleAction(int slice1,int slice2,
		      const Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
  void Read (IOSectionClass &in);
  
  IonIonActionClass(PathDataClass &pathData);
};

class QMCSamplingClass: public ActionBaseClass
{
public:
  int manager;
  int myProc;
  int workerSize;
  int workerList[];
  double SingleAction(int slice1,int slice2,const 
		      Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
  void Read (IOSectionClass &in);

  QMCSamplingClass(PathDataClass &pathData);
};

//double QuickAvg(std::vector<double>* values);
void QuickAvg(std::vector<double>* values, double& avg, double& var);

#endif
