#include "WaterFast.h"
#include "../PathDataClass.h"

double WaterClass::SingleAction (int slice1, int slice2, 
				 const Array<int,1> &activeParticles, int level)
{
  double total=0.0;
  //need to actually do something here to put the stuff from the path into the system
  for (int slice=slice1;slice<=slice2;slice++){
    double factor = ( (slice==slice1) || (slice==slice2)) ? 0.5 : 1.0;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      if (ptcl<54){
	for (int dim=0;dim<NDIM;dim++)
	  we.system.r[3*ptcl].vec[dim]=PathData.Path(slice,ptcl)[dim];
      }
      else if (ptcl<54*2){
	for (int dim=0;dim<NDIM;dim++)
	  we.system.r[3*(ptcl-54)+1].vec[dim]=PathData.Path(slice,ptcl)[dim];
      }
      else {
	for (int dim=0;dim<NDIM;dim++)
	  we.system.r[3*(ptcl-54*2)+2].vec[dim]=PathData.Path(slice,ptcl)[dim];
      }

    }
    we.d.ComputeDistDisp(we.system);
    we.rho.Compute(we.system);

    double en=we.ComputeEnergy();

    //    cerr<<"En is "<<en<<endl;;
    total+=en*factor;
  }
  
  return PathData.Path.tau*total; // *3.1577504e5;
  //    cerr<<we.ComputeEnergy()<<endl;
    
}

void WaterClass::Read(IOSectionClass &in)
{
    we.Init(54);
    we.system.ReadPositions(54);
    SetBox(PathData.Path.GetBox()[0]);
    we.Fast=true;

}
