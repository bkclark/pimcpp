/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
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
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "PIMCClass.h"
#include "MirroredClass.h"
#include "ParseCommand.h"

/// \mainpage An overview of PIMC++ 
///
/// \section intro_sec Introduction
///
/// Understanding PIMC++ can be broken down into understanding four
/// main pieces -- 
/// the moves, the observables, and the actions and how they rest on
/// top of the representation of the path.  
///
///  \section path_sec  The Path
/// There is a file called PathClass.cc that contains all the
/// information relevant to the path. PathClass stores information
/// concerning --the location of all the particles at all places in
/// imaginary time
/// --the species information
/// --the box information
///
/// \subsection species_subsection The Species
/// Every particle is assigned a specific species.  Each species has
/// certain information associated with it. This includes
///
/// \section moves_sec The Moves
/// There are two ways that a move can be built in PIMC++
///

main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  string version = VERSION;
  cout << "pimc++ v. " << version << endl;

  list<ParamClass> argList;
  argList.push_back (ParamClass("verbose", "v", false));
  CommandLineParserClass parser(argList);
  parser.Parse (argc, argv);
  if (parser.NumFiles() != 1)
    cout << "Usage:\n" << "  pimc++ [-v] myfile.in\n";
  else {
    if (parser.Found ("verbose"))
      IO::SetVerbose(true);
    PIMCClass PIMC;
    string inputName = parser.GetFile(0);
    IOSectionClass in;
    if (!in.OpenFile(inputName)) {
      cerr << "Could not open " << inputName << " for reading.  Exitting.\n";
      exit(-1);
    }
    if (PIMC.Read(in)) {
      PIMC.Run();
    } else {
      PIMC.Dummy();
    }
  }

  COMM::Finalize();
}
