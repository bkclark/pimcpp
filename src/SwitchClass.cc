#include "PathDataClass.h"
#include "SwitchClass.h"
#include "time.h"

void SwitchClass::DoEvent()
{
  TimesCalled++;
  if (PathData.Path.Equilibrate)
    PathData.Path.Equilibrate = 0;
  else
    PathData.Path.Equilibrate = 1;
}

void SwitchClass::Read(IOSectionClass &in)
{
  //do nothing for now
}
