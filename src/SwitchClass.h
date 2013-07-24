#ifndef SWITCH_CLASS_H
#define SWITCH_CLASS_H

#include "EventClass.h"

/// This is the generic parent class for all switches, including "real switches"
/// which actually move particles and "pseudo switches", which just shift around
/// data, but don't move anything physical.
class SwitchClass : public EventClass
{
protected:

public:
  //virtual void Read(IOSectionClass &input)=0;
  void DoEvent();
  void Read(IOSectionClass &input);

  /// SwitchClass constructor. Sets reference to the PathData object
  SwitchClass(PathDataClass &pathData, IOSectionClass &out) :
    EventClass(pathData, out)
    {
      // do nothing
    }
};

#endif
