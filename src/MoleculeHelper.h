#ifndef MOLECULE_HELPER_H
#define MOLECULE_HELPER_H

#include "Common.h"
#include "IO/IO.h"


using namespace IO;

class MoleculeManagerClass
{
  int numMolTypes;
  int totalNumMol;
  Array<string,1> names;
  Array<int,1> num;
  Array<Array<int,1>, 1> Members;
  Array<Array<int,1>, 1> ListByType;
  Array<int,1> MolRef;
  Array<string,1> MolLabel;

  public:

  int checkNumParticles;
  void Read(IOSectionClass& in);
  void Init();
  Array<int,1>& MembersOf(int mol);
  void MembersOf(Array<int,1>& members, int mol);
  int SizeOf(int mol);
  int NumMol(int type);
  int NumMol(string typeLabel);
  int NumMol();
  string NameOf(int mol);
  Array<int,1> MolOfType(int type);
  Array<int,1> MolOfType(string typeLabel);
  void MolOfType(Array<int,1>& list, string typeLabel);
  void MolOfType(Array<int,1>& list, int type);
  int Index(string label);
  inline int operator() (int ptcl);

  MoleculeManagerClass();
};

inline int MoleculeManagerClass::operator()(int ptcl)
{
  return MolRef(ptcl);
}

#endif
