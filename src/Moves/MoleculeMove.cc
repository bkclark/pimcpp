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

#include "MoleculeMove.h"

void MoleculeRotate::Read(IOSectionClass &moveInput) {
  cerr << "  Rotate read in" << endl;
  //assert(moveInput.ReadVar("SetTheta",Theta));
  assert(moveInput.ReadVar("SetAngle",Theta));
  doAllSlices = false;
  if(moveInput.ReadVar("AllSlices",doAllSlices))
    cerr << "Rotate over all slices set to " << doAllSlices << endl;
  
	MolMoveClass::Read(moveInput);
  sizeOfMulti = 0;
}

void PairActionTestClass::Read(IOSectionClass &moveInput) {
  assert(moveInput.ReadVar("rstart",r0));
  assert(moveInput.ReadVar("rend",rc));
  assert(moveInput.ReadVar("File",filename));
  assert(moveInput.ReadVar("Species1",spec1));
  assert(moveInput.ReadVar("Species2",spec2));
  int s1num = PathData.Path.SpeciesNum(spec1.c_str());
  int s2num = PathData.Path.SpeciesNum(spec2.c_str());
  ptcl1 = PathData.Path.Species(s1num).FirstPtcl;
  if(s1num == s2num)
    ptcl2 = ptcl1 + 1;
  else
    ptcl2 = PathData.Path.Species(s2num).FirstPtcl;
  cerr << "PairActionTest for species " << spec1 << " " << spec2 << endl;
  cerr << "Using ptcls " << ptcl1 << " " << ptcl2 << endl;
	MolMoveClass::Read(moveInput);
}

void BondStretch::Read(IOSectionClass &moveInput) {
  cerr << "  Stretch read in" << endl;
  assert(moveInput.ReadVar("SetStrain",s));
  doAllSlices = false;
  if(moveInput.ReadVar("AllSlices",doAllSlices))
    cerr << "Stretch over all slices set to " << doAllSlices << endl;
  MolMoveClass::Read(moveInput);
}

void MoleculeTranslate::Read(IOSectionClass &moveInput) {
  cerr << "  Molecule Translate read in" << endl;
  assert(moveInput.ReadVar("SetStep",Step));
  doAllSlices = false;
  if(moveInput.ReadVar("AllSlices",doAllSlices))
    cerr << "Translate over all slices set to " << doAllSlices << endl;
	MolMoveClass::Read(moveInput);
}

void MoleculeMulti::Read(IOSectionClass &moveInput) {
  cerr << "  Multi read in" << endl;
  Rotate.Read(moveInput);
  Trans.Read(moveInput);
  Stretch.Read(moveInput);
  MolMoveClass::Read(moveInput);
}

void DimerMove::Read(IOSectionClass &moveInput) {
  cerr << "  Dimer Move read in" << endl;
  assert(moveInput.ReadVar("SetStep",Step));
	MolMoveClass::Read(moveInput);
  assert(mode == GLOBAL);
}

void ParticleTranslate::Read(IOSectionClass &moveInput) {
  cerr << "  Ptcl Trans read in" << endl;
  assert(moveInput.ReadVar("SetSigma",Sigma));
	MolMoveClass::Read(moveInput);
}

void DummyEvaluate::Read(IOSectionClass &moveInput) {
	MolMoveClass::Read(moveInput);
}

//void MoleculeTranslate::MakeMove() {
double MoleculeTranslate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //cout << "TRANSLATE::SAMPLE" << endl;
  //double step = 0.3;
  double step = Step; // Using Step from input file

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}
  else {
    cerr << "Mode " << mode << " for Translate not supported" << endl;
    exit(1);
  }
  //cout << "chose particles " << activeParticles << endl;
//cerr << counter << ", ";
//  int choosemol = counter;
//	counter = (counter+1)%PathData.Mol.NumMol();
//	activeParticles.resize(MolMembers(choosemol).size());
//	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(choosemol)(i);
//cerr << "  chose molecule " << choosemol << " of " << PathData.Mol.NumMol() << endl;
//cerr << "  Getting info from " << MolMembers(choosemol) << endl;
//cerr << "  activeParticles is " << activeParticles << endl;

//	activeParticles.resize(PathData.Path.NumParticles());
//	for(int x=0; x<activeParticles.size(); x++) activeParticles(x) = x;

  int numSlices = PathData.Path.TotalNumSlices;
  if(doAllSlices) {
    //cerr << "ALL SLICE SAMPLE" << endl;
    slice1 = 0;
    slice2 = numSlices;
	  for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
    	dVec move = TranslateMolAll(PathData.Mol.MembersOf(MoveList(activeMol)),step); 
    }
  }
  else {
    // choose a time slice to move
	  int slice =0;
    slice1 = slice;
    slice2 = slice;
    if(numSlices>1){
      int P_max = numSlices - 1;
      slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
      //slice1 = slice;
      //slice2 = slice;
      slice1 = slice-1;
      slice2 = slice+1;
    }

	  double move_mag_sq = 0.0;
	  for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
	  	//cout << "  Before move: chose slice " << slice << endl;
	  	//for(int i=0; i<activeParticles.size(); i++) cout << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
    	dVec move = TranslateMol(slice,PathData.Mol.MembersOf(MoveList(activeMol)),step); 
	  	//cout << "  After move: of " << move << endl;
	  	//for(int i=0; i<activeParticles.size(); i++) cout << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
    	//move_mag_sq += move(0)*move(0) + move(1)*move(1) + move(2)*move(2);
	  }
  }

  //// this block is a hack to output configurations
  //// in a format for QBox
	//int numMol = PathData.Mol.NumMol();
	//int output1 = 0;
	//int output2 = output1 + 1;
	//int output3 = output2 + 1;
	//if(numMoves == output1){
	//	ofstream before("H2O.64.rigid.global.1.dat");
	//	for(int i=0; i<numMol; i++){
	//		dVec RO = PathData.Path(0,i);
	//		dVec RH1 = PathData.Path(0,i+numMol);
	//		dVec RH2 = PathData.Path(0,i+2*numMol);
	//		before << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
	//		before << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
	//		before << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
	//	}
	//	before.close();
	//}
	//else if(numMoves == output2){
	//	ofstream after("H2O.64.rigid.global.2.dat");
	//	for(int i=0; i<numMol; i++){
	//		dVec RO = PathData.Path(0,i);
	//		dVec RH1 = PathData.Path(0,i+numMol);
	//		dVec RH2 = PathData.Path(0,i+2*numMol);
	//		after << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
	//		after << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
	//		after << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
	//	}
	//	after.close();
	//}
	//else if(numMoves == output3){
	//	ofstream third("H2O.64.rigid.global.3.dat");
	//	for(int i=0; i<numMol; i++){
	//		dVec RO = PathData.Path(0,i);
	//		dVec RH1 = PathData.Path(0,i+numMol);
	//		dVec RH2 = PathData.Path(0,i+2*numMol);
	//		third << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
	//		third << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
	//		third << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
	//	}
	//	third.close();
	//}
  /////////////////////////////////////////////////////////////////////////////


  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current translate ratio is " << double(numAccepted)/numMoves << " with step size " << step << endl;
    //double avg = sqrt(move_mag_sq)/numAccepted;
    //double diff = move_mag_sq/numMoves;
    //cerr << "TRANSLATE diffusion value is " << diff << " avg step is " << avg << endl;
  }
  numMoves++;

//	if(counter == 0)
//		cerr << endl << "I have average action " << UAction << "/" << numMoves << " = " << UAction/numMoves << endl;
	//cerr << "+";

	if(mode == SEQUENTIAL) Advance();

	return 1;
}

double DimerMove::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //cerr << " DimerMove::Sample ";
  // load activeParticles; this move is hard-wired to be a GLOBAL update
	activeParticles.resize(PathData.Path.NumParticles());
	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;

  double step = Step; // Using Step from input file
  // this is only written to work for a dimer!
  assert(MoveList.size() == 2);

  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
    slice1 = slice-1;
    slice2 = slice+1;
  }

  MoveDimerSeparation(slice,PathData.Mol.MembersOf(0), PathData.Mol.MembersOf(1),step); 

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current Dimer Move ratio is " << double(numAccepted)/numMoves << " with step size " << step << " really numAcc is " << numAccepted << endl;
  }
  numMoves++;

	return 1;
}

//void MoleculeRotate::MakeMove()
double MoleculeRotate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //cerr<<"CUMMING Rotating"<<endl;
  double dtheta = M_PI*Theta; // Using Theta from input file
  //cerr << " MoleculeRotate::Sample ";
  //double dtheta = 2*M_PI*0.3;

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}
  else if(mode == MULTIPLE) {
    assert(!doAllSlices);
  }
  else {
    cerr << "Mode " << mode << " for rotate not supported" << endl;
    exit(1);
  }

  Array<int,1> ActiveSlices;
  //int startS, endS;
	int numSlices = PathData.Path.TotalNumSlices;
  if(doAllSlices) {
    //cerr << "ALL SLICE SAMPLE" << endl;
    //startS = 1;
    //endS = numSlices - 1;
    slice1 = 0;
    slice2 = numSlices;
	  for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
	    double theta = 2*(PathData.Path.Random.Local()-0.5)*dtheta;
	    RotateMolXYZAll(PathData.Mol.MembersOf(MoveList(activeMol)), theta);
    }
  }
  else {
	  // choose a time slice to move
	  int slice=0;
	  slice1 = 0;
	  slice2 = 0;
	  if(numSlices>1){
	    int P_max = numSlices - 1;
	    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
      //slice1 = slice;
      //slice2 = slice;
	    slice1 = slice-1;
	    slice2 = slice+1;
	  }

    if(mode == MULTIPLE) {
      int nActiveP = 0;
      vector<int> molToAdd(0);
      vector<int> ptclToAdd(0);
    	int m = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
      molToAdd.push_back(m);
      for(int n=0; n<PathData.Mol.NumMol(); n++) {
        double r_mag;
        dVec R;
	      PathData.Path.DistDisp(slice, m, n, r_mag, R);
        if(r_mag < rc && m!=n) {
          molToAdd.push_back(n);
        }
      }
      MoveList.resize(molToAdd.size());
      for (int n = 0; n<molToAdd.size(); n++) {
        MoveList(n) = molToAdd[n];
        for(int i=0; i<PathData.Mol.MembersOf(MoveList(n)).size(); i++) {
          ptclToAdd.push_back(PathData.Mol.MembersOf(MoveList(n))(i));
          nActiveP ++;
        }
      }
	  	activeParticles.resize(nActiveP);
	  	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = ptclToAdd[i];
      sizeOfMulti += molToAdd.size();
    }
	  //cerr << "Rotate choosing slice " << slice << " and molecule " << MoveList << " wit members " << PathData.Mol.MembersOf(MoveList(0)) << endl;
	  for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
	    //activeParticles.resize(MolMembers(MoveList(activeMol)).size());
	    //for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(activeMol))(i);
	    	
	   	double theta = 2*(PathData.Path.Random.Local()-0.5)*dtheta;
        //cin >> theta;
        //cout << " theta " << theta << endl;
        //// ********************************************
        //Array<int,1> mol1Members;
        //int slice1 = 1;
        //int slice2 = 2;
        //PathData.Mol.MembersOf(mol1Members, 0);
        //dVec p1a = PathData.Path(slice1, mol1Members(1)) - PathData.Path(slice1, mol1Members(0));
        //dVec p1b = PathData.Path(slice1, mol1Members(2)) - PathData.Path(slice1, mol1Members(0));
        //dVec p2a = PathData.Path(slice2, mol1Members(1)) - PathData.Path(slice2, mol1Members(0));
        //dVec p2b = PathData.Path(slice2, mol1Members(2)) - PathData.Path(slice2, mol1Members(0));

        //// compute bisector, norm, and in-plane unit vectors
        //dVec b1 = GetBisector(p1a, p1b);
        //dVec n1 = Normalize(cross(p1a, p1b));
        //dVec u1 = Normalize(cross(b1, n1));
        //dVec b2 = GetBisector(p2a, p2b);
        //dVec n2 = Normalize(cross(p2a, p2b));
        //dVec u2 = Normalize(cross(b2, n2));
        //RotateMol(slice, PathData.Mol.MembersOf(MoveList(activeMol)), n1, theta);
        //cerr << "ROTATE " << theta << " about " << n1 << endl;
        //// **********************************************
	   	RotateMol(slice,PathData.Mol.MembersOf(MoveList(activeMol)), theta);
	    	//RotateMolXYZ(slice,PathData.Mol.MembersOf(MoveList(activeMol)), theta);
        //dVec axis;
        //cin >> axis(0);
        //cin >> axis(1);
        //cin >> axis(2);
        //cout << " A_x " << axis(0) << endl;
        //cout << " A_y " << axis(1) << endl;
        //cout << " A_z " << axis(2) << endl;
	    	//RotateMol(slice,PathData.Mol.MembersOf(MoveList(activeMol)), axis, theta);
    }
  }

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current rotate ratio is " << double(numAccepted)/numMoves << " with angle size " << dtheta << endl;
    if(mode == MULTIPLE)
      cerr << "Avg size of clusters " << sizeOfMulti/double(numMoves) << endl;
  }
  numMoves++;

	if(mode == SEQUENTIAL) Advance();

	return 1;
}

double PairActionTestClass::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  int nPoints = 200;
  double dx = (rc - r0)/float(nPoints);
  cerr << "PATest sample ptcls " << ptcl1 << " " << ptcl2 << endl;
  cerr << "PATest sample " << nPoints << " " << r0 << " " << rc << " " << dx << endl;
  activeParticles.resize(2);
  activeParticles(0) = ptcl1;
  activeParticles(1) = ptcl2;
  ofstream out(filename.c_str());
  out << "# r";
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while(actionIter != Actions.end()) {
    out << " " << (*actionIter)->GetName();
    actionIter++;
  }
  out << endl;

  int numSlices = PathData.Path.TotalNumSlices;
  assert(numSlices > 3);
  int BisectionLevel = 0;
  slice1 = 0;
  slice2 = numSlices - 1;
  //int slice = 1;
  dVec R0((0., 0., 0.));
  dVec R1((0., 0., 0.));
  cerr << "positions " << R0 << " " << R1 << endl;
  R1(0) = r0;
  cerr << "positions " << R0 << " " << R1 << endl;
  for(int s=0; s<=numSlices; s++) {
    PathData.Path.SetPos(s,ptcl1,R0);
    PathData.Path.SetPos(s,ptcl2,R1);
  }
  for(int n=0; n<nPoints; n++) {
    double r = r0 + n*dx;
    R1(0) = r;
    //for(int s=1; s<=2; s++)
    for(int s=0; s<=numSlices; s++)
      PathData.Path.SetPos(s,ptcl2,R1);
    cerr << n << " " << r << " " << R1 << endl;
    out << r;
    actionIter=Actions.begin();
    while (actionIter!=Actions.end()){
      double S = ((*actionIter)->Action(slice1, slice2, activeParticles, BisectionLevel));
      out << " " << S;
      actionIter++;
    }
    out << endl;
  }

  out.close();
  cerr << "Wrote " << filename << "; Done." << endl;
  exit(1);
}

double MoleculeMulti::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //cerr<<"CUMMING MOleculeMulti Rotate"<<endl;
  double r = Rotate.Sample(slice1, slice2, activeParticles);
  double t = Trans.Sample(slice1, slice2, activeParticles);
  double s = Stretch.Sample(slice1, slice2, activeParticles);

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current MULTI ratio is " << double(numAccepted)/numMoves << endl;
    cerr << "NOTE: please disregard acceptance output from individual moves; it will erroneously state 0" << endl;
  }

  numMoves++;

	return(r*s*t); 
}

double ParticleTranslate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  numMoves++;
	//cerr << " ParticleTranslate::Sample" << endl;

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}
  else {
    cerr << "Mode " << mode << " for Translate not supported" << endl;
    exit(1);
  }

  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
    slice1 = slice-1;
    slice2 = slice+1;
  }
	//cerr << numMoves << ": ParticleTranslate moving " << activeParticles << endl;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		for(int activeP = 0; activeP<PathData.Mol.MembersOf(MoveList(activeMol)).size(); activeP++){
			int movePtcl = PathData.Mol.MembersOf(MoveList(activeMol))(activeP);
      //if(Path.ParticleSpeciesNum(movePtcl) != 0){
			  //cerr << "Moving ptcl " <<  movePtcl << " at slice " << slice << " from " << PathData.Path(slice,movePtcl);
  		TranslatePtcl(slice, movePtcl, Sigma);
			  //cerr << " to " << PathData.Path(slice,movePtcl) << endl;
      //}
		}
	}

  if (numMoves%10000 == 0){
    cerr << numMoves << " moves; current PARTICLE translate ratio is " << double(numAccepted)/numMoves << " with step size " << Sigma << endl;
  }
//	if(counter == 0)
//		cerr << endl << "I have average action " << UAction << "/" << numMoves << " = " << UAction/numMoves << endl;
	//cerr << "+";

	if(mode == SEQUENTIAL) Advance();

	//cerr << "Goodbye Ptcl Translate Sample " << endl;
	return 1;
}

double DummyEvaluate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  numMoves++;

  if (numMoves%10000 == 0){
    cerr << numMoves << " moves; current DUMMY accept ratio is " << double(numAccepted)/numMoves << endl;
  }
  //cout << "DUMMY " << numAccepted << " acc out of " << numMoves << " = " << double(numAccepted)/numMoves << " and just for fun " << double(numAccepted)/(numMoves-1) << endl;

	return 1;
}


double BondStretch::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles){

  if(mode == SINGLE){
    int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
    MoveList(0) = choosemol;
    activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
    for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
  }
  else if(mode == SEQUENTIAL){
    activeParticles.resize(PathData.Mol.MembersOf(MoveList(0)).size());
    for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(MoveList(0))(i);
  }
  else if(mode == GLOBAL){
    activeParticles.resize(PathData.Path.NumParticles());
    for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
  }
  else {
    cerr << "Mode " << mode << " for Stretch not supported" << endl;
    exit(1);
  }

  int startSlice = 0;
  int endSlice = 0;
  int numSlices = PathData.Path.TotalNumSlices;
  if(doAllSlices) {
    slice1 = 0;
    slice2 = numSlices;
    startSlice = 0;
    endSlice = numSlices;
    cerr << "STRETCH all slice1 " << slice1 << " " << slice2 << " " << startSlice << " " << endSlice << endl;
  }
  else {
    // choose a time slice to move
    int slice=0;
    slice1 = 0;
    slice2 = 0;
    if(numSlices>1){
      int P_max = numSlices - 1;
      slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
      //slice1 = slice;
      //slice2 = slice;
      slice1 = slice-1;
      slice2 = slice+1;
      startSlice = slice;
      endSlice = slice;
    }
  }

  double theta0 = 1.910611932158; // 109.47 degrees in rad
  //cerr << "Rotate choosing slice " << slice << " and molecule " << MoveList << " wit members " << PathData.Mol.MembersOf(MoveList(0)) << endl;
  for(int molIndex=0; molIndex<MoveList.size(); molIndex++){
    int activeMol = MoveList(molIndex);

    // get relevant bonds from MoveList and PathData.Mol.MembersOf
    vector<int> bondPtcls(0);
    vector<int*> anglePairs(0);
    for(int p1=1; p1<PathData.Mol.MembersOf(activeMol).size(); p1++){
      int ptcl1 = PathData.Mol.MembersOf(activeMol)(p1);
      if(ptcl1 != PathData.Mol(ptcl1)){
        bondPtcls.push_back(ptcl1);
        //cerr << "S adding bound ptcl " << ptcl1 << endl;
      }
      for(int p2=p1+1; p2<PathData.Mol.MembersOf(activeMol).size(); p2++){
        int ptcl2 = PathData.Mol.MembersOf(activeMol)(p2);
        int pair[2];
        pair[0] = ptcl1;
        pair[1] = ptcl2;
        anglePairs.push_back(pair);
        //cerr << "S adding angle pair " << pair[0] << " & " << pair[1] << endl;
      }
    }
    // stretch bonds along their length
    for(int b=0; b<bondPtcls.size(); b++){
      double strain = 2*s*(PathData.Path.Random.Local()-0.5);
      for(int slice=startSlice; slice<=endSlice; slice++) {
        StressBond(slice, bondPtcls[b], activeMol, strain);
      }
    }
    // stretch angles by rotating
    double dtheta = 2*s*theta0*(PathData.Path.Random.Local()-0.5);
    for(int slice=startSlice; slice<=endSlice; slice++) {
      for(int a=0; a<anglePairs.size(); a++){
        dVec v1 = PathData.Path(slice,anglePairs[a][0]) - PathData.Path(slice,activeMol);
        dVec v2 = PathData.Path(slice,anglePairs[a][1]) - PathData.Path(slice,activeMol);
        Path.PutInBox(v1);
        Path.PutInBox(v2);
        //double theta0 = GetAngle(v1,v2);
        //cerr << "S going to stretch angle " << theta0*180/M_PI << " by dtheta which is in deg " << dtheta*180/M_PI << endl;
        dVec u = crossprod(v1,v2);
        //StressAngle(slice, anglePairs[a][1],Normalize(u),dtheta);

        //v1 = PathData.Path(slice,anglePairs[a][0]) - PathData.Path(slice,activeMol);
        //v2 = PathData.Path(slice,anglePairs[a][1]) - PathData.Path(slice,activeMol);
        //Path.PutInBox(v1);
        //Path.PutInBox(v2);
        //double afterTheta = GetAngle(v1,v2);
        //cerr << "S and after stress angle is " << afterTheta << " in deg " << afterTheta*180/M_PI << endl;
      }
    }
  }

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current STRETCH ratio is " << double(numAccepted)/numMoves << endl;
  }
  numMoves++;

  if(mode == SEQUENTIAL) Advance();

  return 1;
}
