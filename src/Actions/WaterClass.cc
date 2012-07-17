#include "WaterFast.h"
#include "../PathDataClass.h"

double WaterClass::SingleAction (int slice1, int slice2, 
				 const Array<int,1> &activeParticles, int level)
{
  t.Start();
  //  cerr<<"Calling for "<<slice1<<" and "<<slice2<<endl;
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

  for (int i=0;i<we.system.r.size();i++){
      for (int j=0;j<we.system.r.size();j++){
 	if (i!=j){
	  double dist;
	  dVecp r12;
	  double dist2;
	  we.d.GetDistDisp(i,j,dist,dist2,r12);
	  if ((we.system.atom[i]!=we.system.atom[j])  &&
	      dist<0.2)
	    return PathData.Path.tau*999;
	}
      }
  }
  

    we.rho.Compute(we.system);

    double en=we.ComputeEnergy();

    //    cerr<<"En is "<<en<<endl;;
    total+=en*factor;
  }
    t.Stop();
    //    cerr<<"Total time was "<<t.Time()<<endl;
    t.Clear();
    //    cerr<<"The value is "<<total<<" "<<PathData.Path.tau*total<<endl;;
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
