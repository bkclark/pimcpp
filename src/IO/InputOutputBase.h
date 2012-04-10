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

#ifndef INPUT_OUTPUT_BASE_H
#define INPUT_OUTPUT_BASE_H

#include <string>
#include <list>
#include <stack>
#include <blitz/array.h>
//#include "../Blitz.h"
#include <fstream>

using namespace std;
//using namespace blitz;
/// This typedef enumerates the atomic datatype that are supported in
/// inputfiles.   Note that 1,2, and 3-dimensional arrays of these
/// types are also supported
typedef enum{NOT_ATOMIC,INT_TYPE,DOUBLE_TYPE,STRING_TYPE,BOOL_TYPE} AtomicType; typedef double scalar;

/// This is the abstract base class for the storage of variables in
/// the IOTreeClass lists.   It contains the name of the variable,
/// it's atomic type (see above), and it's dimensionality.
/// Specializations of this class will also have some reference to the
/// data in the variable itself.
class VarClass 
{  
public:
  string Name;
  AtomicType Type;
  int Dim;
  

  virtual bool ReadInto (double &val)              = 0;
  virtual bool ReadInto (blitz::Array<double,1> &val)     = 0;
  virtual bool ReadInto (blitz::Array<double,2> &val)     = 0;
  virtual bool ReadInto (blitz::Array<double,3> &val)     = 0;
  virtual bool ReadInto (blitz::Array<double,4> &val)     = 0;
  virtual bool ReadInto (int &val)                 = 0;
  virtual bool ReadInto (blitz::Array<int,1> &val)        = 0;
  virtual bool ReadInto (blitz::Array<int,2> &val)        = 0;
  virtual bool ReadInto (blitz::Array<int,3> &val)        = 0;
  virtual bool ReadInto (blitz::Array<int,4> &val)        = 0;
  virtual bool ReadInto (string &val)              = 0;
  virtual bool ReadInto (blitz::Array<string,1> &val)     = 0;
  virtual bool ReadInto (blitz::Array<string,2> &val)     = 0;
  virtual bool ReadInto (blitz::Array<string,3> &val)     = 0; 
  virtual bool ReadInto (blitz::Array<string,4> &val)     = 0; 
  virtual bool ReadInto (bool &val)                = 0;
  virtual bool ReadInto (blitz::Array<bool,1> &val)       = 0;
  virtual bool ReadInto (blitz::Array<bool,2> &val)       = 0;
  virtual bool ReadInto (blitz::Array<bool,3> &val)       = 0;
  virtual bool ReadInto (blitz::Array<bool,4> &val)       = 0;

  virtual bool Append (double val) = 0;
  virtual bool Append (blitz::Array<double,1> &val) = 0;
  virtual bool Append (blitz::Array<double,2> &val) = 0;  
  virtual bool Append (blitz::Array<double,3> &val) = 0;
  virtual bool Append (int val) = 0;
  virtual bool Append (blitz::Array<int,1> &val) = 0;
  virtual bool Append (blitz::Array<int,2> &val) = 0;
  virtual bool Append (blitz::Array<int,3> &val) = 0;
  virtual bool Append (string val) = 0;
  virtual bool Append (blitz::Array<string,1> &val) = 0;
  virtual bool Append (blitz::Array<string,2> &val) = 0;
  virtual bool Append (blitz::Array<string,3> &val) = 0;
  virtual bool Append (bool val) = 0;
  virtual bool Append (blitz::Array<bool,1> &val) = 0;
  virtual bool Append (blitz::Array<bool,2> &val) = 0;
  virtual bool Append (blitz::Array<bool,3> &val) = 0;
  virtual ~VarClass() { /* do nothing */ }
};




/// This class stores a tree of input file sections.  Each member of
/// the tree contains a list of tree nodes below it and a list of
/// variables contained in the present node.
class IOTreeClass
{
protected:
  // USE ME!  I'm not being used yet.
  bool IsModified;
public:
  list<VarClass*> VarList;
  list<IOTreeClass*> SectionList;

  inline void MarkModified();
  /// This is used to ensure proper ordering of sections in the HDF
  /// version in which there is no guarantee that the sections will
  /// come out of the file in the same order you put them in.
  int MyNumber, CurrSecNum;
  virtual void PrintTree()=0;
  virtual void PrintTree(int numIndent)=0;

  IOTreeClass* Parent;
  /// This is the empty string unless I'm the root node of some file. 
  string FileName;
  string Name;
  inline void InsertSection (IOTreeClass *newSec);
  inline bool FindSection (string name, IOTreeClass * &sectionPtr, 
			   int num=0);
  inline int CountSections(string name);
  inline int CountVars();

  template<class T>
  bool ReadVar(string name, T &var)
  {
    bool readVarSuccess;
    list<VarClass*>::iterator varIter=VarList.begin();
    while ((varIter!=VarList.end()) && ((*varIter)->Name!=name)) 
      varIter++;

    bool found = varIter != VarList.end();
    if (found){
      readVarSuccess=(*varIter)->ReadInto(var);
    }
    else if (Parent!=NULL){
      readVarSuccess=Parent->ReadVar(name,var);
    }
    else {
      //cerr<<"Couldn't find variable "<<name<<endl;
      return false;
    }  
    return readVarSuccess;	 
  }

  inline VarClass* GetVarPtr(string name)
  {
    list<VarClass *>::iterator iter = VarList.begin();
    while ((iter != VarList.end()) && ((*iter)->Name != name)){
      iter++;
    }
    if (iter == VarList.end())
      return NULL;
    else {
      MarkModified();  // If we're getting the pointer, we're probably 
                       // gonna modify the variable.
      return *iter;
    }
  }
  inline VarClass* GetVarPtr(int num)
  {
    list<VarClass *>::iterator iter = VarList.begin();
    int i=0;
    while ((iter != VarList.end()) && (i != num)){
      iter++;
      i++;
    }
    if (iter == VarList.end())
      return NULL;
    else {
      MarkModified();  // If we're getting the pointer, we're probably 
                       // gonna modify the variable.
      return *iter;
    }
  }




  /// Write me!
  virtual IOTreeClass* NewSection(string name)=0;
  
  virtual bool OpenFile (string fileName, 
			 string mySectionName, 
			 IOTreeClass *parent) = 0;
  virtual bool NewFile (string fileName,
			string mySectionName,
			IOTreeClass *parent) = 0;
  /// Inserts a new Include directive in the present section.
  virtual void IncludeSection (IOTreeClass *) = 0;
  virtual void CloseFile() = 0;
  virtual void FlushFile() = 0;
  virtual void WriteVar(string name, double val)=0;
  virtual void WriteVar(string name, blitz::Array<double,1> &val)=0;
  virtual void WriteVar(string name, blitz::Array<double,2> &val)=0;
  virtual void WriteVar(string name, blitz::Array<double,3> &val)=0;
  virtual void WriteVar(string name, blitz::Array<double,4> &val)=0;

  virtual void WriteVar(string name, int val)=0;
  virtual void WriteVar(string name, blitz::Array<int,1> &val)=0;
  virtual void WriteVar(string name, blitz::Array<int,2> &val)=0;
  virtual void WriteVar(string name, blitz::Array<int,3> &val)=0;
  virtual void WriteVar(string name, blitz::Array<int,4> &val)=0;

  virtual void WriteVar(string name, string val)=0;
  virtual void WriteVar(string name, blitz::Array<string,1> &val)=0;
  virtual void WriteVar(string name, blitz::Array<string,2> &val)=0;
  virtual void WriteVar(string name, blitz::Array<string,3> &val)=0;
  virtual void WriteVar(string name, blitz::Array<string,4> &val)=0;

  virtual void WriteVar(string name, bool val)=0;
  virtual void WriteVar(string name, blitz::Array<bool,1> &val)=0;
  virtual void WriteVar(string name, blitz::Array<bool,2> &val)=0;
  virtual void WriteVar(string name, blitz::Array<bool,3> &val)=0;
  virtual void WriteVar(string name, blitz::Array<bool,4> &val)=0;

  /// Append a value to a variable of dimension of 1 higher than val.
  /// i.e. Add a double to an blitz::Array<double,1> or add blitz::Array<double,1>
  /// to an blitz::Array<double,2>, etc.
  template<class T>
  inline bool AppendVar(string name, T val);

  inline IOTreeClass() : Parent(NULL)
  { FileName="";}
};

void IOTreeClass::MarkModified()
{
  if (FileName != "")
    IsModified = true;
  else if (Parent!=NULL)
    Parent->MarkModified();
}


template<class T>
inline bool IOTreeClass::AppendVar(string name, T val)
{ 
  VarClass *var = GetVarPtr(name);
  if (var == NULL)
    return false;
  MarkModified();
  return var->Append(val);
}



/// Returns the number of subsections with the given name within the
/// present section.
inline int IOTreeClass::CountSections(string name)
{
  list<IOTreeClass*>::iterator sectionIter;
  sectionIter=SectionList.begin();
  int numSections=0;
  while (sectionIter!=SectionList.end()){
    if ((name==(*sectionIter)->Name) || (name == "")){
      numSections++;
    }
    sectionIter++;
  }
  return numSections;
}

inline int IOTreeClass::CountVars()
{
  list<VarClass*>::iterator varIter;
  varIter=VarList.begin();
  int numVars=0;
  while (varIter!=VarList.end()){
    numVars++;
    varIter++;
  }
  return numVars;
}





/// FindSection locates a subsection with the given name within the
/// section in contains and returns it in the pointer, sectionPtr,
/// which is passed a reference.  Returns true if the section is
/// found.  The final parameter, which default value "true",
/// optionally resets the section iterator to the beginning of the
/// section.  Thus, one may control whether or not order is
/// significant.  
inline bool IOTreeClass::FindSection (string name, 
				      IOTreeClass* &sectionPtr,
				      int num)
{
  
  list<IOTreeClass*>::iterator Iter=SectionList.begin();
  int counter=0;
  while(counter<=num && Iter!=SectionList.end()){
    if ((*Iter)->Name==name){
      counter++;
    }
    if (counter<=num)
      Iter++;
  }
  bool found = Iter != SectionList.end(); 
  if (found){
    sectionPtr = *Iter;
  }
  return (found); 
}



inline void IOTreeClass::InsertSection(IOTreeClass *newSec)
{
  list<IOTreeClass *>::iterator iter;
  
  if (SectionList.empty())
    SectionList.push_back(newSec);
  else
    {
      iter = SectionList.begin();
      while ((iter != SectionList.end()) && 
	     ((*iter)->MyNumber < newSec->MyNumber))
	iter++;
      if (iter!=SectionList.end())
	SectionList.insert(iter, newSec);
      else
	SectionList.push_back(newSec);
    }
}





class OutputSectionClass
{
public:
  virtual bool OpenFile(string fileName)                 = 0;
  virtual void CloseFile()                               = 0;
  virtual void OpenSection(string name)                  = 0;
  virtual void CloseSection()                            = 0;
  virtual void WriteVar(string name, double T)           = 0;
  virtual void WriteVar(string name, blitz::Array<double,1> &v) = 0;
  virtual void WriteVar(string name, blitz::Array<double,2> &v) = 0;
  virtual void WriteVar(string name, blitz::Array<double,3> &v) = 0;
  virtual void WriteVar(string name, int T)              = 0;
  virtual void WriteVar(string name, blitz::Array<int,1> &v)    = 0;
  virtual void WriteVar(string name, blitz::Array<int,2> &v)    = 0;
  virtual void WriteVar(string name, blitz::Array<int,3> &v)    = 0;
  virtual void WriteVar(string name, string str)         = 0;
  virtual void WriteVar(string name, blitz::Array<string,1> &v) = 0;
  virtual void WriteVar(string name, blitz::Array<string,2> &v) = 0;
  virtual void WriteVar(string name, blitz::Array<string,3> &v) = 0;

};



#endif
