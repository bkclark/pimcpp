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

#include "ClusterMove.h"
#include "MoveUtils.h"

int firsttotal = 0;
int total = 0;
int clustmax = 0;

void LocalFlip::AssignPtcl(int mol,Array<int,1>& activeParticles){
  for (int i = 0;i<5;i++)
    activeParticles(i) = mol + PathData.Mol.NumMol()*i;
}

double LocalFlip::MolPairAction(int slice,int m,int n){
  // This function generates the activeParticles array for the two specified molecules and feeds it to the action function.
  Array<int,1> activeParticles(10);
  // load values for molecules m and n
  for (int i = 0;i<5;i++){
    activeParticles(i) = m + PathData.Mol.NumMol()*i;
    activeParticles(i+5) = n + PathData.Mol.NumMol()*i;
  }
  //return PathData.Actions.MoleculeInteractions.Action(slice, slice, activeParticles, 0);
  return MolAction->Action(slice, slice, activeParticles, 0);
}

// Original: replacing w/ coordinate inversion
/*
void LocalFlip::RotateMol(int slice,int mol,dVec axis,double phi){
  // generate list of particles to rotate
  Array<int,1> activeParticles(5);
  AssignPtcl(mol,activeParticles);
  for(int i=0;i<5;i++){
    dVec old_coord=PathData.Path(slice,activeParticles(i));
    old_coord -= PathData.Path(slice,activeParticles(0));
    dVec new_coord = ArbitraryRotate(axis,old_coord,phi);
    new_coord += PathData.Path(slice,activeParticles(0));
    PathData.Path.SetPos(slice,activeParticles(i),new_coord);
  }
}*/

// Trying coordinate reflection like for GlobalFlip, but have to subtract O vector I guess.
void LocalFlip::RotateMol(int slice,int mol,dVec Q){
  dVec aligned;
  dVec perp;
  Array<int,1> activeParticles(5);
  AssignPtcl(mol,activeParticles);
  int Optcl = activeParticles(0);
//cerr << "Rotating Molecule " << Optcl<< " about " << Q << endl;
  dVec O = PathData.Path(slice,Optcl);
  for(int i=1;i<5;i++){
    dVec r;
    double rmag;
    PathData.Path.DistDisp(slice,Optcl,activeParticles(i),rmag,r);
    Strip(Q,r,aligned,perp);
    dVec q = aligned - perp;
//cerr << "  ptcl " << activeParticles(i) << ": before " << r << " and after " << q << endl;
    PathData.Path.SetPos(slice,activeParticles(i),(q+O));
  }
}

//ArbitraryRotate
// rotates a coordinate, coord, a specified angle, phi, about a specified axis, axis.  axis and coord are specified as Cartesian coordinates with respect to a common origin; axis should be a unit vector
dVec LocalFlip::ArbitraryRotate(dVec axis,dVec coord, double phi){
  dVec coord_perp;
  dVec coord_aligned;
  Strip(axis,coord,coord_aligned,coord_perp);
  double perp_mag = Mag(coord_perp);
//cerr << "stripped vector " << coord << " off " << axis << " to give aligned " << coord_aligned << " and perp " << coord_perp << endl;
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
//cerr << "then new coord_perp is " << newcoord << endl;
  return (Scale(newcoord,perp_mag) + coord_aligned);
}

void LocalFlip::Strip(dVec R, dVec u,dVec &aligned, dVec &perp){
  //double phi = PathData.Actions.TIP5PWater.GetAngle(R,u);
  //aligned = PathData.Actions.TIP5PWater.Scale(R,PathData.Actions.TIP5PWater.Mag(u)*cos(phi));
  aligned = Scale(R,dotprod(R,u,1.0));
  perp = u - aligned;
}

void LocalFlip::IntegrityCheck(int slice, Array<int,1> activeParticles){
  // constants and tolerances for comparison
  double Amargin = 0.001;
  double Lmargin = 0.001;
  double TetAngle = 1.9106;
  double Hlength = 1.0;
  double Elength = 0.8;
  bool pass = true;

  // get vectors WRT origin (oxygen) to determine molecule's geometry
  Array<dVec,1> V;
  V.resize(activeParticles.size());
  for (int i = 0; i < activeParticles.size(); i++){
    V(i) = PathData.Path(slice,activeParticles(i));
    if (i > 0){
      V(i) -= V(0);
    }
  }

  // test geometry
  // check bond lengths
  for (int i = 1; i < activeParticles.size(); i++){
    double length = Mag(V(i));
    double complength;
    if (i == 1 || i == 2)
      complength = Elength;
    else if (i == 3 || i == 4)
      complength = Hlength;
    if (abs(complength - length) > Lmargin){
      pass = false;
      cerr << activeParticles(0) << ": Integrity Check Failed." << endl;
      cerr << "  Bond length for ptcl " << activeParticles(i) << " " << length << " differs from " << complength << endl;
    }
  }

  // check angular separation 
  for (int i = 1; i < activeParticles.size(); i++){
    for (int j = i+1; j < activeParticles.size(); j++){
      double theta = GetAngle(V(i),V(j));
      if (abs(TetAngle - theta) > Amargin){
        pass = false;
        cerr << activeParticles(0) << ": Integrity Check Failed." << endl;
        cerr << "  Angle between  ptcls " << activeParticles(i) << "  and " << activeParticles(j) << " " << theta << " differs from " << TetAngle << endl;
      }
    }
  }
}

void LocalFlip::Read(IOSectionClass &moveInput) {
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="LocalFlip");
  assert(moveInput.ReadVar("name",Name));
  string actionLabel;
  assert(moveInput.ReadVar("ActionName",actionLabel));
  ActionBaseClass* newAction = PathData.Actions.GetAction(actionLabel);
  cerr << "MOLECULEBIAS ACCEPTED ACTION LABELELD " << actionLabel << " AS ACTION TYPE MOLECULEINTERACTIONSCLASS; MAKE SURE THIS IS VALID" << endl;
  assert(typeid(newAction) == typeid(MoleculeInteractionsClass*));
  MolAction = static_cast<MoleculeInteractionsClass*> (newAction);
}

void LocalFlip::MakeMove()
{
  int numMol = PathData.Mol.NumMol();
  double theta = M_PI;
  Array<int,1> ActiveParticles(5);
  Array<bool,1> AlreadyMoved(numMol);
  AlreadyMoved = false;
  Array<bool,1> InCluster(numMol);
  InCluster = false;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numMol);
  Array<int,1> Stack(numMol);
  Stack = -1;
  int now = 0;
  int next = 0;
  // Generate a random unit vector: axis of rotation
  dVec A;
  A(0) = PathData.Path.Random.Local()-0.5;
  A(1) = PathData.Path.Random.Local()-0.5;
  A(2) = PathData.Path.Random.Local()-0.5;
  A = Normalize(A); 

  // For now make this a classical move
  int slice = 0;

  // Rotate the first molecule and add to stack
  AssignPtcl(choosemol,ActiveParticles);
  RotateMol(slice,choosemol,A);
  PathData.AcceptMove(slice,slice,ActiveParticles);
  Stack(next) = choosemol;
  next++;
  AlreadyMoved(choosemol) = true;
  InCluster(choosemol) = true;

  // Proceed through stack (indexed by now) building cluster
  while (now<Stack.size() && Stack(now)!=-1){
    for (int m = 0;m < numMol; m++){
      if (!(AlreadyMoved(m))){
        AssignPtcl(m,ActiveParticles);
        double oldAction = MolPairAction(slice,Stack(now),m);
        RotateMol(slice,m,A);
        double newAction = MolPairAction(slice,Stack(now),m);
        if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
          PathData.AcceptMove(slice,slice,ActiveParticles);
          Stack(next) = m;
          next++; 
          AlreadyMoved(m) = true;
          InCluster(m) = true;
        }
        else{
          PathData.RejectMove(slice,slice,ActiveParticles);
        }
        if (TimesCalled%10 == 0){
          IntegrityCheck(slice, ActiveParticles);
        }
      }
    }
    now++;
  }//

  // Collect some data
  int firstclustsize = 0;
  for (int n = 0;n < InCluster.size();n++){
    if(InCluster(n))
      firstclustsize++;
  }
cerr << "First Moved a cluster of size " << firstclustsize << endl;
  if (firstclustsize > clustmax)
    clustmax = firstclustsize;
  firsttotal += firstclustsize;
/*///////////////////////////////////////////////

  // For explicit detailed balance, attempt the same move again
  // Re-initialize data storage
  AlreadyMoved = false;
  InCluster = false;
  Stack = -1;
  now = 0;
  next = 0;

  AssignPtcl(choosemol,ActiveParticles);
  RotateMol(slice,choosemol,A);
  PathData.AcceptMove(slice,slice,ActiveParticles);
  Stack(next) = choosemol;
  next++;
  AlreadyMoved(choosemol) = true;
  InCluster(choosemol) = true;

  while (now<Stack.size() && Stack(now)!=-1){
    for (int m = 0;m < numMol; m++){
      if (!(AlreadyMoved(m))){
        AssignPtcl(m,ActiveParticles);
        double oldAction = MolPairAction(slice,Stack(now),m);
        RotateMol(slice,m,A);
        double newAction = MolPairAction(slice,Stack(now),m);
        if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
          PathData.AcceptMove(slice,slice,ActiveParticles);
          Stack(next) = m;
          next++; 
          AlreadyMoved(m) = true;
          InCluster(m) = true;
        }
        else{
          PathData.RejectMove(slice,slice,ActiveParticles);
        }
        if (TimesCalled%10 == 0){
          IntegrityCheck(slice, ActiveParticles);
        }
      }
    }
    now++;
  }//


  // Collect some data
  int clustsize = 0;
  for (int n = 0;n < InCluster.size();n++){
    if(InCluster(n))
      clustsize++;
  }
cerr << "Then Moved a cluster of size " << clustsize << endl;
  if (clustsize > clustmax)
    clustmax = clustsize;
  total += clustsize;
*///////////////////////////////////////////

  if (TimesCalled%100 == 0){
    cerr << TimesCalled << " moves; Average first cluster size is " << firsttotal/TimesCalled << endl;
//    cerr << TimesCalled << " moves; Average return cluster size is " << total/TimesCalled << endl;
    cerr << "Maximum cluster size is " << clustmax << endl;
  }
}

void GlobalFlip::Read(IOSectionClass &moveInput) {
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="GlobalFlip");
  assert(moveInput.ReadVar("name",Name));
  string actionLabel;
  assert(moveInput.ReadVar("ActionName",actionLabel));
  ActionBaseClass* newAction = PathData.Actions.GetAction(actionLabel);
  cerr << "MOLECULEBIAS ACCEPTED ACTION LABELELD " << actionLabel << " AS ACTION TYPE MOLECULEINTERACTIONSCLASS; MAKE SURE THIS IS VALID" << endl;
  assert(typeid(newAction) == typeid(MoleculeInteractionsClass*));
  MolAction = static_cast<MoleculeInteractionsClass*> (newAction);
}

void GlobalFlip::RotateMol(int slice,int mol,dVec Q){
  Array<int,1> activeParticles(5);
  AssignPtcl(mol,activeParticles);
  for(int i=0;i<5;i++){
    dVec r = PathData.Path(slice,activeParticles(i));
    dVec q = (2.0*Q) - (1.0*r);
    //    dVec q = Q;
    //    q=0.5*q;
    //    q=q  - r;
    PathData.Path.SetPos(slice,activeParticles(i),q);
  }
}

void GlobalFlip::AssignPtcl(int mol,Array<int,1>& activeParticles){
  for (int i = 0;i<5;i++)
    activeParticles(i) = mol + PathData.Mol.NumMol()*i;
}

double GlobalFlip::MolPairAction(int slice,int m,int n){
  // This function generates the activeParticles array for the two specified molecules and feeds it to the action function.
  Array<int,1> activeParticles(10);
  // load values for molecules m and n
  for (int i = 0;i<5;i++){
    activeParticles(i) = m + PathData.Mol.NumMol()*i;
    activeParticles(i+5) = n + PathData.Mol.NumMol()*i;
  }
  //return PathData.Actions.MoleculeInteractions.Action(slice, slice, activeParticles, 0);
  return MolAction->Action(slice, slice, activeParticles, 0);
}

void GlobalFlip::MakeMove()
{
  int numMol = PathData.Mol.NumMol();
  double theta = M_PI;
  Array<int,1> ActiveParticles(5);
  Array<bool,1> AlreadyMoved(numMol);
  AlreadyMoved = false;
  Array<bool,1> InCluster(numMol);
  InCluster = false;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numMol);
  Array<int,1> Stack(numMol);
  Stack = -1;
  int now = 0;
  int next = 0;

  // Generate a random vector in the box to specify point of reflection
  dVec A;
  for (int x = 0; x < 3; x++){
    A(x) = (PathData.Path.Random.Local()-0.5)*PathData.Path.GetBox()(x);
  }

  // For now make this a classical move
  int slice = 0;

  // Rotate the first molecule and add to stack
  AssignPtcl(choosemol,ActiveParticles);
  RotateMol(slice,choosemol,A);
  PathData.AcceptMove(slice,slice,ActiveParticles);
  Stack(next) = choosemol;
  next++;
  AlreadyMoved(choosemol) = true;
  InCluster(choosemol) = true;

  // Proceed through stack (indexed by now) building cluster
  while (now<Stack.size() && Stack(now)!=-1){
    for (int m = 0;m < numMol; m++){
      if (!(AlreadyMoved(m))){
        AssignPtcl(m,ActiveParticles);
        double oldAction = MolPairAction(slice,Stack(now),m);
        RotateMol(slice,m,A);
        double newAction = MolPairAction(slice,Stack(now),m);
        if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
          PathData.AcceptMove(slice,slice,ActiveParticles);
          Stack(next) = m;
          next++; 
          AlreadyMoved(m) = true;
          InCluster(m) = true;
        }
        else{
          PathData.RejectMove(slice,slice,ActiveParticles);
        }
      }
    }
    now++;
  }

  // Collect some data
  int clustsize = 0;
  for (int n = 0;n < InCluster.size();n++){
    if(InCluster(n))
      clustsize++;
  }
//cerr << "Moved a cluster of size " << clustsize << endl;
  if (clustsize > clustmax)
    clustmax = clustsize;
  total += clustsize;

  if (TimesCalled%20 == 0){
    cerr << TimesCalled << " moves; Average cluster size is " << total/TimesCalled << endl;
    cerr << "Maximum cluster size is " << clustmax << endl;
  }
}
