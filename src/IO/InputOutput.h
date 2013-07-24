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

#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include "InputOutputBase.h"
#include "InputOutputHDF5.h"
#include "InputOutputASCII.h"


#include <stack>


/// In the file name format name.extn, returns the extension.
/// Actually returns everything after the trailing.
inline string Extension (string fileName)
{
  string extn;
  stack<char> bwExtn;
  int pos = fileName.length()-1;
  while ((pos >= 0) && fileName[pos]!='.') {
    bwExtn.push(fileName[pos]);
    pos--;
  }
  
  if (fileName[pos] == '.') 
    while (!bwExtn.empty()) {
      extn += bwExtn.top();
      bwExtn.pop();
    }
  else
    extn = "";
  return (extn);
}


/// This function takes a filename, determines it extension, creates a
/// new IOTreeASCIIClass or IOTreeHDF5Class based on the
/// extension, and calls OpenFile on the new object.
/// Extensions:  
/// .h5:            HDF5
/// .xml:           XML
/// .anything_else  ASCII
inline IOTreeClass *ReadTree (string fileName, 
			      string myName,
			      IOTreeClass *parent)
{
  IOTreeClass *newTree;
  string extn = Extension (fileName);
  if (extn == "h5")
    newTree = new IOTreeHDF5Class;
  else
    newTree = new IOTreeASCIIClass;
  
  newTree->FileName = fileName;
  bool success = newTree->OpenFile (fileName, myName, parent);
  cerr << "success = " << success << endl;
  if (success)
    return (newTree);
  else{
    delete newTree;
    return (NULL);
  }
}


inline IOTreeClass *NewTree (string fileName,
			     string myName,
			     IOTreeClass *parent)
{
  IOTreeClass *newTree;
  string extn = Extension (fileName);
  if (extn == "h5")
    newTree = new IOTreeHDF5Class;
  //  else if (extn == "xml")
  //    newTree = newIOTreeXMLClass;
  else
    newTree = new IOTreeASCIIClass;

  bool success = newTree->NewFile (fileName, myName, parent);
  if (success)
    return (newTree);
  else{
    delete newTree;
    return (NULL);
  }
}


///  Wrapper class for IOTreeClass that gives a nearly identical
///  interface as the OutputSectionClass.
class IOSectionClass
{
private:
  IOTreeClass *CurrentSection;
public:

  /// Opens the file reference by fileName and reads the contents into
  /// the tree in CurrentSection.  Creates a new object based on the
  /// extnesion of the filename.  For ".h5", it creates an
  /// IOTreeHDF5Class.  For ".xml" it creaes an IOTreeXMLClass.
  /// After creating the object, it calls the objects virtual OpenFile
  /// function, reading the contents of the file into the tree.
  bool OpenFile (string fileName);
  string GetName(){ return CurrentSection->Name;}
  inline string GetFileName();
  string GetVarName(int num){ return GetVarPtr(num)->Name;}
  /// Creates a file at the top level, choosing the appropriate type
  /// based on the file extension.
  bool NewFile (string fileName);

  /// Calls CurrentSections close file and then deletes the
  /// CurrentSection.  
  void CloseFile ();

  /// Flush all buffers to disk for safety
  inline void FlushFile();

  /// Opens the num'th section with the given name.  The default
  /// value for num is 0.
  inline bool OpenSection (string name, int num=0);

  /// Opens the num'th section below CurrentSection.
  inline bool OpenSection (int num);

  /// This mounts a file in the current tree under CurrentSection at
  /// the end of CurrentsSection's SectionList.  It does not change
  /// what CurrentSection points to, ie. it does not descend to the
  /// newly-opened section.
  inline bool IncludeSection (string name, string fileName);

  /// Creates a new section of the same type as currentSection under
  /// currentSection.  Pushes the new section to the end of the
  /// section list.
  inline void NewSection (string name)
  {
    CurrentSection = CurrentSection->NewSection(name);
  }

  /// This function creates a new file of the appropriate type as
  /// determined by the extension of fileName and mounts it at the end
  /// of the list under CurrentSection.  Returns false if the file
  /// couldn't be created.
  inline bool NewSection (string name, string fileName);

  /// Closes the current section.  That is, CurrentSection becomes
  /// CurrentSection's parent.
  inline void CloseSection ();

  /// Template function which reads a variable in the present section
  /// into the passed-by-reference T variable.
  template<class T>
  bool ReadVar(string name, T &var)
  {  return (CurrentSection->ReadVar(name, var)); }

  template<class T>
  bool ReadVar(string name, T &var, T Default)
  { 
    bool success = ReadVar(name, var);
    if (!success)
      var = Default;
    return (success);
  }

  /// Writes a variable under the current section.
  /* template <class T>
  void WriteVar(string name, T val)
  { CurrentSection->WriteVar(name, val); }*/

  inline void WriteVar (string name, double val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<double,1> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<double,2> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<double,3> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<double,4> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, int val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<int,1> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<int,2> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<int,3> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<int,4> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, bool val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<bool,1> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<bool,2> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<bool,3> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<bool,4> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, const char *val)
  { CurrentSection->WriteVar(name, (string)val);}
  inline void WriteVar (string name, string val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<string,1> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<string,2> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<string,3> &val)
  { CurrentSection->WriteVar(name, val); }
  inline void WriteVar (string name, blitz::Array<string,4> &val)
  { CurrentSection->WriteVar(name, val); }
  
  template<class T>
  bool AppendVar(string name, T val)
  { return CurrentSection->AppendVar(name, val); }

  inline VarClass *GetVarPtr(string name)
  {
    return (CurrentSection->GetVarPtr(name));
  }

  inline VarClass *GetVarPtr(int num)
  {
    return (CurrentSection->GetVarPtr(num));
  }



  /// Returns the number of subsections within the present section
  /// which have the name name.  If called without a name, it returns
  /// the total number of sections.
  inline int CountSections(string name="")
  { return (CurrentSection->CountSections(name)); }
  inline int CountVars()
  {return (CurrentSection->CountVars());}
  /// Calls CurrentSections virtual PrintTree() function.  This is for
  /// debugging purposes.  It spits out a hierarchy of the sections
  /// and variable names.
  void PrintTree()
  { CurrentSection->PrintTree(); }

  IOSectionClass() 
  {
    CurrentSection=NULL;
  }
	
};



inline bool IOSectionClass::OpenFile (string fileName)
{
  CurrentSection = ReadTree (fileName, "Root", NULL);
  cerr << "After ReadTree.\n";
  if (CurrentSection == NULL)
    return (false);
  else
    return (true);
}

inline bool IOSectionClass::NewFile (string fileName)
{
  CurrentSection = NewTree (fileName, "Root", NULL);
  if (CurrentSection == NULL)
    return (false);
  else
    return (true);
}


inline void IOSectionClass::CloseFile ()
{
  while (CurrentSection->Parent != NULL)
    CloseSection();
  CurrentSection->FlushFile();
  CurrentSection->CloseFile();
  delete (CurrentSection);
}


inline void IOSectionClass::FlushFile()
{
  IOTreeClass *tree = CurrentSection;
  while (tree->Parent != NULL)
    tree = tree->Parent;
  tree->FlushFile();
}


inline bool IOSectionClass::OpenSection (string name, int num)
{
  IOTreeClass *newSection;
  bool success;
  success = CurrentSection->FindSection(name, newSection, num);
  if (success)
    CurrentSection=newSection;
  return success;
}


inline bool IOSectionClass::OpenSection (int num)
{
  IOTreeClass *newSection;
  list<IOTreeClass*>::iterator Iter=CurrentSection->SectionList.begin();
  int i = 0;
  while ((i<num) && 
	 (Iter != CurrentSection->SectionList.end())){
    i++;
    Iter++;
  }
  if (i<num)
    return false;
  else {
    CurrentSection = *Iter;
    return true;
  }
}


inline bool IOSectionClass::IncludeSection (string name, string fileName)
{
  IOTreeClass *newSection;
  newSection = ReadTree (fileName, name, CurrentSection);
  if (newSection == NULL)
    return false;
  else {
    CurrentSection->IncludeSection(newSection);
    return true;
  }
}

///Don't think this pushes to back of list like it should nor does newfile
inline bool IOSectionClass::NewSection (string name, string fileName)
{
  IOTreeClass *newSection;
  newSection = NewTree (fileName, name, CurrentSection);
  if (newSection == NULL)
    return false;
  else {
    CurrentSection->IncludeSection(newSection);
    CurrentSection = newSection;
    return true;
  }
}

inline void IOSectionClass::CloseSection ()
{
  //cerr << "Closing Section " << CurrentSection->Name << endl;
  assert (CurrentSection->Parent != NULL);
  CurrentSection = CurrentSection->Parent;
}


inline string
IOSectionClass::GetFileName()
{
  IOTreeClass *tree = CurrentSection;
  while (tree->Parent != NULL)
    tree = tree->Parent;
  return tree->FileName;
}


#endif
