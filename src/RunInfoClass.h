/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef RUN_INFO_H
#define RUN_INFO_H

#include <Common/IO/IO.h>
#include <unistd.h>
#include <cstdio>
#include <sys/types.h>
#include <pwd.h>
#include <time.h>

#ifndef COMMONVERSION
#define COMMONVERSION "unknown"
#endif

#ifndef VERSION
#define VERSION "unknown"
#endif


using namespace IO;

class RunInfoClass
{
public:
  string ProgramName;
  string Version;
  string CommonVersion;
  string UserName;
  string RunTime;
  string BuildDate;
  string BuildTime;
  string HostName;
  inline void Write(IOSectionClass &outSec)
  {
    outSec.WriteVar("ProgramName", ProgramName);
    outSec.WriteVar("Version", Version);
    outSec.WriteVar("CommonVersion", CommonVersion);
    outSec.WriteVar("UserName", UserName);
    outSec.WriteVar("RunTime", RunTime);
    outSec.WriteVar("BuildDate", BuildDate);
    outSec.WriteVar("BuildTime", BuildTime);
    outSec.WriteVar("HostName", HostName);
  }
  inline void Read(IOSectionClass &inSec)
  {
    assert (inSec.ReadVar("ProgramName", ProgramName));
    assert (inSec.ReadVar("Version", Version));
    if (!inSec.ReadVar ("CommonVersion", CommonVersion))
      CommonVersion = "unknown";
    assert (inSec.ReadVar("UserName", UserName));
    assert (inSec.ReadVar("RunTime", RunTime));
    assert (inSec.ReadVar("BuildDate", BuildDate));
    assert (inSec.ReadVar("BuildTime", BuildTime));
    assert (inSec.ReadVar("HostName", HostName));
  }
  inline RunInfoClass()
  {
    //struct passwd* pwInfo = getpwuid(getuid());
    UserName = "ethan"; //pwInfo->pw_name;
    BuildDate = __DATE__;
    BuildTime = __TIME__;
    char hostname[300];
    gethostname(hostname, 300);
    HostName = hostname;
    time_t seconds = time(NULL);
    RunTime = ctime(&seconds);
    Version = VERSION;
    CommonVersion = COMMONVERSION;
  }
};


#endif
