#include "ParseCommand.h"

CommandLineParserClass::CommandLineParserClass(list<ParamClass> &argList)
{
  list<ParamClass>::iterator iter;
  for (iter=argList.begin(); iter!=argList.end(); iter++) 
    ArgMap[iter->GetName()] = (*iter);
  for (iter=argList.begin(); iter!=argList.end(); iter++)
    if (iter->GetShortName() != "")
      ShortMap[iter->GetShortName()] = (*iter);
}

bool
CommandLineParserClass::Parse(int argc, char **argv)
{
  for (int i=1; i<argc; i++) {
    if ((argv[i][0]=='-') && (argv[i][1]=='-')) {
      string name = &(argv[i][2]);
      map<string, ParamClass>::iterator iter;
      iter = ArgMap.find(name);
      if (iter != ArgMap.end()) {
	ArgMap[name].Found = true;
	if (ArgMap[name].NeedsArg) {
	  if ((i+1) < argc) 
	    ArgMap[name].SetArg(argv[++i]);
	  else
	    return false;
	}
      }
      else {
	cerr << "Unrecognized argument """ << name << """" << endl;
	return false;
      }
    }
    else if (argv[i][0] == '-') {
      string name = &(argv[i][1]);
      map<string, ParamClass>::iterator iter;
      iter = ShortMap.find(name);
      if (iter != ArgMap.end()) {
	ShortMap[name].Found = true;
	if (ShortMap[name].NeedsArg) {
	  if ((i+1) < argc) 
	    ShortMap[name].SetArg(argv[++i]);
	  else
	    return false;
	}
      }
      else {
	cerr << "Unrecognized short argument """ << name << """" << endl;
	return false;
      } 
    }
    else 
      Files.push_back (argv[i]);
  }
  return true;
}
