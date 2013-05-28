#ifndef READ_FROM_FILE_ACTION_CLASS_H
#define READ_FROM_FILE_ACTION_CLASS_H

#include "ActionBase.h"

class ReadFromFileActionClass: public ActionBaseClass
{
  ifstream infile;
  int qCount;
  bool safe_mode;
  double conversion;
	public:

	double SingleAction(int slice1,int slice2,
		      const Array<int,1> &activeParticles,int level);
	double ComputeEnergy(int slice1,int slice2,
					const Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
  void Read (IOSectionClass &in);
  
  ReadFromFileActionClass(PathDataClass &pathData);
};

#endif
