#include "../PathDataClass.h"
#include "ReadAction.h"

std::string
ReadFromFileActionClass::GetName()
{
  return "ReadFromFileAction";
}
	
ReadFromFileActionClass::ReadFromFileActionClass(PathDataClass &pathData) : ActionBaseClass (pathData){
  qCount = 0;
}
	
double ReadFromFileActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  if(safe_mode) {
    assert(slice1==slice2);
    assert(slice1==0);
  }
  qCount++;
	double Utotal = 0.0;
  Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
	return Utotal*PathData.Path.tau;
}

double ReadFromFileActionClass::ComputeEnergy(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  double Utotal = 0.0;
  //cerr << "pwscf Compute Energy " << slice1 << " " << slice2 << endl;
  infile >> Utotal; 
  return (conversion * Utotal);
}

double ReadFromFileActionClass::d_dBeta (int slice1, int slice2, int level){
	//int slice = 0;
  cerr << "ReadFromFileAction dBeta returning null" << endl;
  double Utotal = -999;
	return Utotal;
}

void ReadFromFileActionClass::Read (IOSectionClass &in){
	cerr << "ReadFromFileAction Read" << endl;
  string filename;
  assert(in.ReadVar("File",filename));
  infile.open(filename.c_str());
  cerr << "ReadFRomFileAction Read finished." << endl;
  safe_mode = true;
  in.ReadVar("Safe",safe_mode);
  conversion = 1.0;
  in.ReadVar("Conversion",conversion);
}
