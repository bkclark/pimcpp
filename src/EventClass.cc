#include "PathDataClass.h"
#include "EventClass.h"

EventClass::EventClass(PathDataClass &pathData, IOSectionClass& out) : 
  PathData(pathData), IOSection(out), TimeSpent(0.0), TimesCalled(0),
  Path(pathData.Path)
{
    // do nothing else for now
  //		cerr << "EventClass construct" << endl;
}
