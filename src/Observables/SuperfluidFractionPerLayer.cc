#include "SuperfluidFractionPerLayer.h"


void 
SuperfluidFractionPerLayerClass::Read(IOSectionClass& IO)
{
  WindingNumberClass::Read(IO);
  assert(SpeciesList.size()==1);
  ParticleAlreadyCounted.resize(PathData.Path.Permutation.size());
  int numGridPoints=200;
  grid.Init(-PathData.Path.GetBox()[2]/2.0,PathData.Path.GetBox()[2]/2.0,numGridPoints);
  Histogram.resize(numGridPoints);
}

//Gets the winding of the loop that ptcl is in for dimension dim.
//Also updates the ParticleAlreadyCounted array with a flag
//that indicates which particles have been touched
double 
SuperfluidFractionPerLayerClass::GetWindingNumber(int ptcl,int dim,
						  int &cycleLength)
{
  int currentPtcl=ptcl;
  double dist=0.0;
  cycleLength=0;
  do{
    cycleLength=cycleLength+1;
    ParticleAlreadyCounted(ptcl)=true;
    for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++){
      dVec disp=PathData.Path.Velocity(slice,slice+1,currentPtcl);
      dist=dist+disp[dim];
    }
    currentPtcl= PathData.Path.Permutation(currentPtcl);
  }
  while (currentPtcl!=ptcl);
  return round(dist/PathData.Path.GetBox()[dim]);
  
}

void SuperfluidFractionPerLayerClass::BinZWindingNumber(double ptclWinding,int ptcl,
							  int dim,int cycleLength)
{
  double windingPerSlice=ptclWinding/(PathData.Path.TotalNumSlices*cycleLength);
  int currentPtcl=ptcl;
  do{
    for (int slice=0;slice<PathData.Path.TotalNumSlices;slice++){
      dVec disp=PathData.Path(slice,currentPtcl);
      PathData.Path.PutInBox(disp);
      if (disp[2]<grid.End && disp[2]>grid.Start) {
	int index=grid.ReverseMap(disp[2]);
	Histogram(index)=Histogram(index)+windingPerSlice;
      }
    }
    currentPtcl= PathData.Path.Permutation(currentPtcl);
  }
  while (currentPtcl!=ptcl);



}

void SuperfluidFractionPerLayerClass::Accumulate()
{
  SamplesInBlock++;
  // Move the join to the end so we don't have to worry about
  // permutations 
  int cycleLength;
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  for (int dim=0;dim<NDIM;dim++){
    ParticleAlreadyCounted=false;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      if (!ParticleAlreadyCounted(ptcl)){
        double ptclWinding=GetWindingNumber(ptcl,dim,cycleLength);
	if (ptclWinding!=0)
	  BinZWindingNumber(ptclWinding,ptcl,dim,cycleLength);
      }
    }
  }


}


void
SuperfluidFractionPerLayerClass::WriteBlock()
{
  int species=SpeciesList(0);
  double beta=PathData.Path.tau*PathData.Path.TotalNumSlices;
  int numParticles=
    PathData.Path.Species(species).LastPtcl-PathData.Path.Species(species).FirstPtcl;
  double factor=1.0/((double)SamplesInBlock*(2*PathData.Path.Species(species).lambda*beta*numParticles));
  

  // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type",string("Vector"));
    }
    Histogram *=factor;
    HistogramVar.Write(Histogram);
  }
  Histogram=0.0;
  SamplesInBlock=0;
}

