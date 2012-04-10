 
#include "../PathDataClass.h"
#include "DiagonalActionClass.h"
#include "ctime"

string 
DiagonalActionClass::GetName()
{
  return "DiagonalActionClass";
}


double 
DiagonalActionClass::d_dBeta (int slice1, int slice2, int level)
{
  assert(1==2);
}


double 
DiagonalActionClass::SingleAction (int slice1, int slice2,
			       const Array<int,1> &changedParticles,
			       int level)
{
  cerr<<"Calling diagonal action"<<endl;
  double TotalU=0.0;
  PathClass &Path = PathData.Path;
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;

  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
      if (Path.DoPtcl(ptcl2)){
	int species2=Path.ParticleSpeciesNum(ptcl2);
	PairActionFitClass &PA = *(PairMatrix(species1, species2));
	for (int slice=slice1;slice<=slice2;slice+=skip){
	  dVec r;
	  double rmag;
	  PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
	  double factor=1.0;
	  if (slice==slice1 || slice==slice2)
	    factor=0.5;
	  //	  double U=factor*PA.U(rmag,0,0, level);
	  double U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	  //	  double U = factor*PA.Udiag(rmag, level);
	  //	  double U = factor*((DavidPAClass*)(&PA)).Udiag_exact(rmag,level);
	  // Subtract off long-range part from short-range action
	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	    U -=  (PA.Ulong(level)(rmag));
	  TotalU+=U;
	}
      }
    }
  }
//  cerr<<"Num total U is "<<TotalU<<endl;
////  cerr<<"Worm short range action is "<<TotalU<<endl;
//  cerr<<"I'm out of the short range action"<<endl;
  return (TotalU);
}
void 
DiagonalActionClass::Read(IOSectionClass &in)
{

}

DiagonalActionClass::DiagonalActionClass(PathDataClass &pathData,
					 Array<PairActionFitClass* ,2> 
					 &pairMatrix) : 

  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{

}

