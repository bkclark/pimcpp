#ifndef SUPERFLUID_FRACTION_PER_LAYER_H
#define SUPERFLUID_FRACTION_PER_LAYER_H

#include "ObservableBase.h"
#include "WindingNumber.h"



class SuperfluidFractionPerLayerClass : public WindingNumberClass
{
 protected:
  ObservableVecDouble1 HistogramVar;
  Array<bool,1> ParticleAlreadyCounted;
  void BinZWindingNumber(double ptclWinding,int ptcl,
			 int dim,int cycleLength);
  double GetWindingNumber(int ptcl,int dim,int &cycleLength);

  void Accumulate();
  Array<double,1> Histogram;
 public:
  LinearGrid grid;
  void Read(IOSectionClass& IO);
  void WriteBlock();
  SuperfluidFractionPerLayerClass(PathDataClass &myPathData,IOSectionClass &ioSection) :
    WindingNumberClass(myPathData,ioSection),
    HistogramVar("SuperfluidFractionPerLayer", IOSection, myPathData.Path.Communicator)
  {
    W2Sum = 0.0;
    SamplesInBlock = 0;
  }

};



#endif
