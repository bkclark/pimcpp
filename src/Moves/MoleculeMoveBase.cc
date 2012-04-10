#include "MoleculeMoveBase.h"
#include "MoveUtils.h"

// default constructor loads MolMembers array of arrays
MolMoveClass::MolMoveClass(PathDataClass& myPathData, IO::IOSectionClass outSection):
  LocalStageClass (myPathData,outSection)
{
	cerr << "MolMoveClass constructor" << endl;
	cerr << "	numMol is " << numMol << endl;
	numAccepted = 0;
	numMoves = 0;

	mode = SINGLE;

  startIndex = 0;
  numActionsToRead = 1;
 
	/* 
  cerr << "In MolMoveClass Constructor.  MolMembers loaded: " << endl;
  for(int i = 0; i < MolMembers.size() ; i++){
    cerr << i << ":  ";
    for(int j = 0; j < MolMembers(i).size(); j++) cerr << MolMembers(i)(j) << ", ";
    cerr << endl;
  }*/
}

//MolMoveClass::MolMoveClass(PathDataClass& myPathData, IO::IOSectionClass outSection, int numToRead, int start):
//  LocalStageClass (myPathData,outSection)
//{
//	numAccepted = 0;
//	numMoves = 0;
//
//	mode = SINGLE;
//
//  startIndex = start;
//  numActionsToRead = numToRead;
//}

void MolMoveClass::Read (IOSectionClass &in){
	/// Read in the molecule to move by string or int
	cerr << "MolMoveClass::Read" << endl;
	string setMolecule;
	int indexMolecule;
	if(in.ReadVar("Molecule",setMolecule)){
		cerr << "got Molecule " << setMolecule << endl;
		molIndex = PathData.Mol.Index(setMolecule);
		numMol = PathData.Mol.NumMol(molIndex);
	}
	else if (in.ReadVar("Molecule",indexMolecule)){
    molIndex = indexMolecule;
		numMol = PathData.Mol.NumMol(molIndex);
	}
	else{
		molIndex = 0;
		numMol = PathData.Mol.NumMol(molIndex);
    cerr << "BY DEFAULT, ";
	}
	cerr << "Selected molecule index " << molIndex << "(" << PathData.Mol.NameOf(molIndex) << ") with " << numMol << endl;
	
/// Read in update mode: GLOBAL, SEQUENTIAL, SINGLE	
	string setMode;
	assert(in.ReadVar("Mode", setMode));
	if(setMode == "GLOBAL"){
		mode = GLOBAL;
    //MoveList = PathData.Mol.MolOfType(molIndex);
    PathData.Mol.MolOfType(MoveList, molIndex);
    cerr << "GLOBAL move init MoveList is " << MoveList << endl;
		//MoveList.resize(numMol);
		//for(int m=0; m<numMol; m++)
		//	MoveList(m) = m + PathData.Path.offset[molIndex];
		//cerr << "Set for GLOBAL particle updates." << endl;
	}
	else if(setMode == "SEQUENTIAL"){
		mode = SEQUENTIAL;
		MoveList.resize(1);
		MoveList(0) = 0;
		cerr << "Set for SEQUENTIAL particle updates." << endl;
	}
	else if(setMode == "MULTIPLE"){
		mode = MULTIPLE;
    assert(in.ReadVar("Radius",rc));
		cerr << "Set for MULTIPLE particle updates." << endl;
	}
	else{
		cerr << "Set for SINGLE particle updates." << endl;
		MoveList.resize(1);
	}
		
  Array<string,1> ActionList;
  assert (in.ReadVar ("Actions", ActionList));
  //cerr << "  Looking for " << numActionsToRead << " actions starting at index " << startIndex << endl;
  //assert ((numActionsToRead + startIndex) <= ActionList.size());
  cerr << "Read in actions " << ActionList << endl;
	for(int a=0; a<ActionList.size(); a++){
		string setAction = ActionList(a);
		//string setAction = ActionList(a+startIndex);
    // hack: needs to be integrated into new scheme!
    //if(setAction == "ShortRange") {
    //  cerr << "TEMP HACK: ADDING ACTIONS.SHORTRANGE TO MOLECULEMOVEBASE.ACTIONS" << endl;
    //  Actions.push_back(&PathData.Actions.ShortRange);
    //}
    if(setAction == "DavidLongRange") {
      cerr << "TEMP HACK: ADDING ACTIONS.DAVIDLONGRANGE TO MOLECULEMOVEBASE.ACTIONS" << endl;
      Actions.push_back(&PathData.Actions.DavidLongRange);
    }
    else {
      ActionBaseClass* newAction = PathData.Actions.GetAction(setAction);
      Actions.push_back(newAction);
      cerr << a+1 << " of " << ActionList.size() << ": Added action with label " << setAction << " and address " << newAction << endl;
    }
	}
}

void MolMoveClass::LoadActions(list<ActionBaseClass*> actions)
{
  cerr << "Wiping out actions loading with " << actions.size();
  Actions = actions;
  cerr << " " << Actions.size() << "... ok, done" << endl;
}

dVec MolMoveClass::GetCOM(int slice, int mol){
  return PathData.Path(slice,PathData.Mol.MembersOf(mol)(0));
}

// translate members of a molecule by a specified vector translate
dVec MolMoveClass::TranslateMol(int slice, Array<int,1>& activePtcls, dVec translate){
  for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
    dVec newP = PathData.Path(slice,activePtcls(ptcl)) + translate;
    PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
  }
	return translate;
}

// Generates a translation vector and moves all
// particles in the molecule
dVec MolMoveClass::TranslateMol(int slice, Array<int,1>& activePtcls, double epsilon){
  dVec translate;
  for(int i = 0; i<3; i++) translate[i] = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
  //cin >> translate(0);
  //cin >> translate(1);
  //cin >> translate(2);
  //cout << "trans_x " << translate(0) <<endl;
  //cout << "trans_y " << translate(1) <<endl;
  //cout << "trans_z " << translate(2) <<endl;
  //cerr << "  translate by " << translate << endl;
  for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
    dVec newP = PathData.Path(slice,activePtcls(ptcl)) + translate;
    //cerr << "  " << ptcl << " " << activePtcls(ptcl) << " " << newP;// << endl;
    //cerr << newP << endl;
    PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
  }
	///////////////////// hackity
	//cerr << "TranslateMol:: I'm just doing this as a favor.  Fixe me right away!!!!!!" << endl;
	//ofstream dimerOut("dimer.in");
	//cerr << "mol is " << activePtcls(0) << " with members " << activePtcls << endl;
	//for(int ptcl=0; ptcl<activePtcls.size(); ptcl++){
	//	dimerOut << PathData.Path(slice,activePtcls(ptcl)) << endl;
	//}
	//dimerOut << endl;
	//dVec myTranslate(2.85, 0.0, 0.0);
  //for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
  //  dVec newP = PathData.Path(slice,activePtcls(ptcl)) + myTranslate;
	//	dimerOut << newP << endl;
  //}
	//////////////////// end hackity
	return translate;
}

// Generates a translation vector and moves all
// particles in the molecule over a range of slices INCLUSIVE
dVec MolMoveClass::TranslateMolAll(Array<int,1>& activePtcls, double epsilon){
  dVec translate;
  for(int i = 0; i<3; i++) translate[i] = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
  int numS = PathData.Path.TotalNumSlices;
  for(int slice=0; slice<=numS; slice++) {
    for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
      dVec newP = PathData.Path(slice,activePtcls(ptcl)) + translate;
      PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
    }
  }
	return translate;
}

void MolMoveClass::TranslatePtcl(int slice, int ptcl, double Sigma){
    dVec translate;
    PathData.Path.Random.CommonGaussianVec (Sigma, translate);

    dVec newP = PathData.Path(slice,ptcl) + translate;
    PathData.Path.SetPos(slice, ptcl, newP);
}

// This will orchestrate the rotation of molecule mol
// about a specified axis, axis, by an angle theta
void MolMoveClass::RotateMol(int slice, Array<int,1>& activePtcls, dVec& axis, double theta){
  //cerr << "RotateMol slice " << slice << " mol " << activePtcls << " about " << axis << " theta " << theta << endl;
	int mol = activePtcls(0);
  axis = Normalize(axis);
  // find COM vector
  dVec O = GetCOM(slice, mol);
	// nonsense to rotate the ptcl at the origin; "O", so loop starts from 1 not 0
  for(int ptcl = 1; ptcl<activePtcls.size(); ptcl++){
    dVec P = PathData.Path(slice, activePtcls(ptcl)) - O;
    //cerr << "  " << ptcl << " " << activePtcls(ptcl) << " old " << P;
    PathData.Path.PutInBox(P);
    //cerr << " putinbox " << P;
    if(Mag(P) > 1e-7) {
      P = ArbitraryRotate(axis, P, theta) + O;
      PathData.Path.SetPos(slice,activePtcls(ptcl),P);
    }
    //cerr << " new " << P << endl;
  }
}

// This will orchestrate and execute an arbitrary rotation
// of a molecule using ArbitraryRotate
// i.e. it generates its own axis of rotation
void MolMoveClass::RotateMol(int slice, Array<int,1>& activePtcls, double theta){
  // generate a unit vector on the unit sphere!
  dVec axis;
  double x = 2 * (PathData.Path.Random.Local()-0.5);
  double A = sqrt(1 - x*x);
  double y = 2 * A * (PathData.Path.Random.Local()-0.5);
  if((PathData.Path.Random.Local()-0.5)<0)
    y *= -1;
  double z = sqrt(1 - x*x - y*y);
  if((PathData.Path.Random.Local()-0.5)<0)
    z *= -1;
  //cout << x << " " << y << " " << z << endl;
  assert(abs(x*x + y*y + z*z - 1) < 1e-5);
  axis(0) = x; axis(1) = y; axis(2) = z;
  RotateMol(slice, activePtcls, axis, theta);
}

// Rotation of molecule about x-, y-, or z-axis, randomly chosen
void MolMoveClass::RotateMolXYZ(int slice, Array<int,1>& activePtcls, double theta){
  int x,y;
  int z = (int)floor(3*PathData.Path.Random.Local());
  if (z == 0){
    x = 1;
    y = 2;
  }
  else if (z == 1){
    x = 2;
    y = 0;
  }
  else if (z == 2){
    x = 0;
    y = 1;
  }
	int mol = activePtcls(0);
  // find COM vector
  dVec O = GetCOM(slice, mol);
	// nonsense to rotate the ptcl at the origin; "O", so loop starts from 1 not 0
  for(int ptcl = 1; ptcl<activePtcls.size(); ptcl++){
    dVec P = PathData.Path(slice, activePtcls(ptcl)) - O;
    PathData.Path.PutInBox(P);
    dVec newP = RotateXYZ(P, x, y, z, theta) + O;
    PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
  }
}

// Rotation of molecule about x-, y-, or z-axis, randomly chosen
// rotates all slices of the molecule about a common COM
void MolMoveClass::RotateMolXYZAll(Array<int,1>& activePtcls, double theta){
  int x,y;
  int z = (int)floor(3*PathData.Path.Random.Local());
  if (z == 0){
    x = 1;
    y = 2;
  }
  else if (z == 1){
    x = 2;
    y = 0;
  }
  else if (z == 2){
    x = 0;
    y = 1;
  }
  int numS = PathData.Path.TotalNumSlices;
	int mol = activePtcls(0);
  // find COM vector
  int norm = 1;
  dVec O_0 = PathData.Path(0, mol);
  dVec COM = O_0;
  for(int slice=1; slice<numS; slice++) {
    dVec O = GetCOM(slice, mol) - O_0;
    PathData.Path.PutInBox(O);
    COM += O + O_0;
    norm += 1;
    //cerr << slice << " " << O << " " << COM << endl;
  }
  COM *= 1./norm;
  for(int slice=0; slice<=numS; slice++) {
    for(int ptcl = 0; ptcl<activePtcls.size(); ptcl++){
      dVec P = PathData.Path(slice, activePtcls(ptcl)) - COM;
      PathData.Path.PutInBox(P);
      dVec newP = RotateXYZ(P, x, y, z, theta) + COM;
      PathData.Path.PutInBox(newP);
      PathData.Path.SetPos(slice,activePtcls(ptcl),newP);
    }
  }
}

void MolMoveClass::MoveDimerSeparation(int slice, Array<int,1> mol1, Array<int,1> mol2, double Sigma){
  dVec translate = PathData.Path(slice,mol1(0)) - PathData.Path(slice,mol2(0));
  PathData.Path.PutInBox(translate);
  //cerr << "D vector between O is " << translate << endl;
  double mag = PathData.Path.Random.LocalGaussian(Sigma);
  //cerr << "D random # is " << mag << endl;
  translate = Renormalize(translate, mag);
  //translate = Scale(translate, mag);
  //cerr << "D rescaled vector is " << translate << endl;
  
  assert(mol1.size() == mol2.size());

  for(int i=0; i<mol1.size(); i++){
    //cerr << "D ptcl " << mol1(i) << " move from " << PathData.Path(slice,mol1(i));
    dVec newR = PathData.Path(slice,mol1(i)) + translate;
    PathData.Path.SetPos(slice, mol1(i), newR);
    //cerr << " to " << PathData.Path(slice,mol1(i)) << endl;
    //cerr << "D ptcl " << mol2(i) << " move from " << PathData.Path(slice,mol2(i));
    newR = PathData.Path(slice,mol2(i)) - translate;
    PathData.Path.SetPos(slice, mol2(i), newR);
    //cerr << " to " << PathData.Path(slice,mol2(i)) << endl;
  }
}

void MolMoveClass::StressBond(int slice, int ptcl, int mol, double s){
  //cerr << "Bond stretch for " << ptcl << " and COM " << mol << endl;
  dVec O = PathData.Path(slice,mol);
  dVec v0 = PathData.Path(slice,ptcl) - O;
  PathData.Path.PutInBox(v0);
  //cerr << "  Initial coords are " << O << " and " << v0 << " wrt O " << endl;
  double scale = 1 + s;
  //dVec v = Scale(v0,scale);
  dVec v = v0*scale;
  //cerr << "  scaled " << v0 << " by " << scale << " to " << v << endl;
  PathData.Path.SetPos(slice, ptcl, O+v);
  //cerr << "  Final coords are " << O << " and " << v << " wrt O " << endl;
}

void MolMoveClass::StressAngle(int slice, int ptcl, dVec axis, double theta){
  int mol = PathData.Mol(ptcl);
  dVec O = PathData.Path(slice,mol);
  dVec v = PathData.Path(slice,ptcl) - O;
  PathData.Path.PutInBox(v);
  v = ArbitraryRotate(axis, v, theta);
  PathData.Path.SetPos(slice, ptcl, O+v);
}

void MolMoveClass::Accept(){
	numAccepted++;
}

void MolMoveClass::Advance(){
	MoveList(0) = MoveList(0) + 1;
	if(MoveList(0) >= PathData.Mol.NumMol())
		MoveList(0) = 0;
}

// ArbitraryRotate takes a unit vector, axis, about
// which to rotate and coord, a coordinate WRT
// the origin of axis, and an angle phi by which
// to rotate.  Then, it returns the rotated coordinate

// Is this really the simplest algorithm to do a rotation
// about an arbitrary axis????
dVec ArbitraryRotate(dVec axis,dVec coord, double phi){
  dVec coord_perp;
  dVec coord_aligned;
  Strip(axis,coord,coord_aligned,coord_perp);
  double perp_mag = Mag(coord_perp);
  coord =  Normalize(coord_perp);
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double gamma = coord(1)/coord(2) - axis(1)/axis(2);
  double eta = coord(0)/coord(2) - axis(0)/axis(2);
  double x = (cosphi/(gamma*coord(2)) - sinphi*axis(2)/coord(0))/(coord(1)/coord(0) + eta/gamma);
  double y = coord(1)*x/coord(0) + sinphi*axis(2)/coord(0);
  double z = -axis(0)*x/axis(2) - axis(1)*y/axis(2);
  dVec newcoord;
  newcoord(0) = x;
  newcoord(1) = y;
  newcoord(2) = z;
  return (Renormalize(newcoord,perp_mag) + coord_aligned);
  //return (Scale(newcoord,perp_mag) + coord_aligned);
}

// Rotate a coodinate about the x-, y-, or z-axis
// indices u1, u2, u3 should be a permutation
// of [0,1,2] to specify the axis about which
// to rotate
dVec RotateXYZ(dVec coord,int u1,int u2,int u3,double theta)
{
  double x0 = coord[u1];
  double y0 = coord[u2];
  double costheta = cos(theta);
  double sintheta = sin(theta);
  dVec new_coord;
  new_coord[u1] = x0*costheta - y0*sintheta;
  new_coord[u2] = y0*costheta + x0*sintheta;
  new_coord[u3] = coord[u3];
  return new_coord;
}
