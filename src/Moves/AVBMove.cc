#include "AVBMove.h"
#include "MoveUtils.h"

AVBMove::AVBMove(PathDataClass &myPathData,IOSectionClass outSection) : 
	MolMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
{
	vol = PathData.Path.GetVol();
}

void AVBMove::Read(IOSectionClass &moveInput){
	assert(moveInput.ReadVar("Rout", Rout));
	assert(moveInput.ReadVar("Rin", Rin));

	RCubeDiff = Rout*Rout*Rout - Rin*Rin*Rin;
	Vin = 4.0*M_PI/3.0*RCubeDiff;
	Vout = vol - Vin;

	MolMoveClass::Read(moveInput);
}

bool AVBMove::AreBound(int slice, int mol1, int mol2){
	bool bound = false;
	dVec O1 = GetCOM(slice,mol1);
	dVec O2 = GetCOM(slice,mol2);
  dVec Roo = O2 - O1;
	double Rmag = Mag(Roo);
	if(Rmag < Rout && Rmag>Rin) bound = true;
	return bound;
}

bool AVBMove::AreBound(int slice, int mol1, dVec coord){
	bool bound = false;
	dVec O1 = GetCOM(slice,mol1);
  dVec Roo = O1 - coord;
	double Rmag = Mag(Roo);
	if(Rmag < Rout && Rmag>Rin) bound = true;
	return bound;
}

dVec AVBMove::GenerateBound(int slice, int ToMol, int FromMol){
	dVec move, disp;
	double rmag;
	for(int i=0; i<3; i++) move(i) = (PathData.Path.Random.Local() - 0.5);
	double scale = pow((PathData.Path.Random.Local()*RCubeDiff + Rin*Rin*Rin),1.0/3);
	move = Scale(move, scale/Mag(move));
	PathData.Path.DistDisp(slice,FromMol,ToMol,rmag,disp);
	return (move + disp);
}

dVec AVBMove::GenerateUnBound(int slice, int ToMol, int FromMol){
	dVec move;
	bool isBound = true;
	while(isBound){
		for(int i=0; i<3; i++) move(i) = (PathData.Path.Random.Local() - 0.5)*PathData.Path.GetBox()(i);
		isBound = AreBound(slice,ToMol,move);
	}
	return (move - PathData.Path(slice,FromMol));
}

double AVBMove::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles){
  //double step = Step; // Using Step from input file
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciesp = PathData.Path.SpeciesNum("p");
  int speciese = PathData.Path.SpeciesNum("e");

	// choose molecule i and j for "swap" move
  int choosemoli = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
  int choosemolj = (int)floor(PathData.Path.Random.Local()*PathData.Mol.NumMol());
	activeParticles.resize(PathData.Mol.MembersOf(choosemoli).size());
	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Mol.MembersOf(choosemoli)(i);

  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    slice1 = slice-1;
    slice2 = slice+1;
  }

	// choose to move to "bound" state or not
	bool ToBound = false;
	if(PathData.Path.Random.Local() < Pbias) ToBound = true;

	// set transition probability
	double T = 1.0;
	bool areBound = AreBound(slice, choosemoli, choosemolj);
	if(areBound && !ToBound)
		T = Pbias*Vout/((1-Pbias)*Vin);
	else if(!areBound && ToBound)
		T = (1-Pbias)*Vin/(Pbias*Vout);

	dVec move;
	if(ToBound){
		move = GenerateBound(slice, choosemoli, choosemolj);
	}
	else{
		move = GenerateUnBound(slice, choosemoli, choosemolj);
	}

	for(int i=0; i<PathData.Mol.MembersOf(choosemoli).size(); i++){
		int ptcl = PathData.Mol.MembersOf(choosemoli)(i);
		dVec newPos = PathData.Path(slice,ptcl) + move;
		PathData.Path.SetPos(slice,ptcl,newPos);
	}

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current AVB ratio is " << double(numAccepted)/numMoves << endl;
	}
	numMoves++;

	return T;
}
