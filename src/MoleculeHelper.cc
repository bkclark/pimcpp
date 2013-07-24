#include "MoleculeHelper.h"

MoleculeManagerClass::MoleculeManagerClass()
{

}

void MoleculeManagerClass::Read(IOSectionClass& in)
{
  cerr << "Reading Molecule information" << endl;
  string tempName;
  if(in.ReadVar("Name",tempName)){
    numMolTypes = 1;
    names.resize(1);
    names(0) = tempName;
    assert(in.ReadVar("NumOfMolecules",totalNumMol));
    num.resize(1);
    num(0) = totalNumMol;
  }
  else{
    assert(in.ReadVar("NumMolTypes",numMolTypes));
    assert(in.ReadVar("Names",names));
    assert(names.size() == numMolTypes);
    assert(in.ReadVar("NumOfMolecules",num));
    assert(num.size() == numMolTypes);
    totalNumMol = 0;
    for(int n=0; n<numMolTypes; n++){
      totalNumMol += num(n);
    }
  }


  cerr << "Going to read array of molecule IDs for each particle" << endl;
  // initialize a bunch of stuff
  int max = 0;
  int count = 0;
  MolLabel.resize(totalNumMol);
  Members.resize(totalNumMol);
  for(int m=0; m<totalNumMol; m++)
    Members(m).resize(0);
  ListByType.resize(numMolTypes);
  cerr << "SIZE OF LIST IS " << ListByType.size() << endl;
  int total = 0;
  for(int t=0; t<numMolTypes; t++){
    ListByType(t).resize(num(t));
    for(int x=0; x<ListByType(t).size(); x++)
      ListByType(t)(x) = total + x;
    total += num(t);
    cerr << "  ARRAY " << t << " is " << ListByType(t) << endl;
  }

  assert(in.ReadVar("MoleculeIDs",MolRef));
  checkNumParticles = MolRef.size();
  cerr << "Read in array of size " << MolRef.size() << endl;
  for(int p=0; p<MolRef.size(); p++){
    cerr << p << " belongs to ";
    cerr << MolRef(p) << endl;
    // get max mol id for a sanity check
    if(MolRef(p) > max)
      max = MolRef(p);

    // initialize array with molecule labels
    int molCount = num(0);
    int index = 0;
    while(MolRef(p) >= molCount){
      cerr << "Huh? incrementing index" << endl;
      index++;
      molCount += num(index);
    }
    //if(count == num(index)){
    //  cerr << "count " << count << " = " << num(index) << endl;
    //  count = 0;
    //  index++;
    //}
    //count ++;
    //cerr << "Problem is...";
    //cerr << " name " << names(index);
    //cerr << " index " << index;
    //cerr << " size " << ListByType(index).size();
    MolLabel(MolRef(p)) = names(index);
    //cerr << " MolLabel assigned " << MolLabel(MolRef(p));
    //int typeSize = ListByType(index).size();
    //ListByType(index).resizeAndPreserve(typeSize+1);
    //ListByType(index)(typeSize) = MolRef(p);
    //cerr << "added molecule " << p << " to list: " << ListByType(index) << endl;
    //cerr << "  added label" << endl;

    // initialize arrays with ptcl ids indexed by molecule
    int myMol = MolRef(p);
    int s = Members(myMol).size();
    Members(myMol).resizeAndPreserve(s+1);
    Members(myMol)(s) = p;
    cerr << "  added member " << p << " to mol " << myMol << " of size " << s+1 << endl;
  }
  assert(max == (totalNumMol-1));
}

void MoleculeManagerClass::Init()
{
  // not sure what to have here...
}

void MoleculeManagerClass::MembersOf(Array<int,1>& members, int mol)
{
  assert(mol < totalNumMol);
  members.resize(Members(mol).size());
  for(int s=0; s<Members(mol).size(); s++)
    members(s) = Members(mol)(s);
}

Array<int,1>& MoleculeManagerClass::MembersOf(int mol)
{
  assert(mol < totalNumMol);
  return Members(mol);
}

int MoleculeManagerClass::SizeOf(int mol)
{
  assert(mol < totalNumMol);
  return Members(mol).size();
}

int MoleculeManagerClass::NumMol()
{
  return totalNumMol;
}

int MoleculeManagerClass::NumMol(int type)
{
  return num(type);
}

int MoleculeManagerClass::NumMol(string typeLabel)
{
  int found = Index(typeLabel);
  return num(found);
}

string MoleculeManagerClass::NameOf(int mol)
{
  assert(mol < totalNumMol);
  return MolLabel(mol);
}

void MoleculeManagerClass::MolOfType(Array<int, 1>& list, int type)
{
  cerr << "Returning MolOfType " << ListByType(type) << endl;
  list.resize(ListByType(type).size());
  for (int s=0; s<ListByType(type).size(); s++)
    list(s) = ListByType(type)(s);
}

void MoleculeManagerClass::MolOfType(Array<int, 1>& list, string typeLabel) 
{
  int found = Index(typeLabel);
  cerr << "MolOfType found MOlID " << found << " for mol " << typeLabel << endl;
  list.resize(ListByType(found).size());
  for (int s=0; s<ListByType(found).size(); s++)
    list(s) = ListByType(found)(s);
}

Array<int,1> MoleculeManagerClass::MolOfType(int type)
{
  cerr << "Returning MolOfType " << ListByType(type) << endl;
  return ListByType(type);
}

Array<int,1> MoleculeManagerClass::MolOfType(string typeLabel)
{
  int found = Index(typeLabel);
  cerr << "MolOfType found MOlID " << found << " for mol " << typeLabel << endl;
  return ListByType(found);
}

int MoleculeManagerClass::Index(string label)
{
  for(int found=0; found<names.size(); found++){
    if(label == names(found)){
      return found;
    }
  }
  cerr << "MoleculeManager ERROR: molecule label " << label << " not found" << endl;
  assert(0);
  return -1;
}
