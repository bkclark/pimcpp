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


#include "SPS.h"


void SPSClass::Read(IOSectionClass& in)
{
  cerr << "Stochastic Potential Switching Move read..." << endl;
  NumSteps = 1;
  in.ReadVar("NumSteps",NumSteps);
  //assert(in.ReadVar("Vmax",VMAX));
  Tf2s = Ts2f = 0;
  in.ReadVar("LogTransProb",Tf2s); 
  in.ReadVar("LogReverseTransProb",Ts2f); 
  int stages = in.CountSections("Stage");
  //in.ReadVar("NumMoveStages",stages);
	//Array<string,1> methodList;
	//assert(in.ReadVar("MoveMethod",methodList));
  //assert(methodList.size() == stages);
  //Array<int,1> numActions(stages);
  //numActions = 1;
  //assert(in.ReadVar("NumActions", numActions));
  //assert((numActions.size()-1) == stages);
  //int startIndex = 0;
  for(int s=0; s<stages; s++){
    MolMoveClass* fullMoveStage;
    MolMoveClass* switchMoveStage;
    assert(in.OpenSection("Stage",s));
    string method;
    assert(in.ReadVar("MoveMethod",method));
    //int actionsToRead = numActions(s);
    // adding this; actionstoRead should be deprecated!!
    //actionsToRead = 0;
    cerr << "  Init " << s+1 << " of " << stages << " stages: " << method << endl;
	  if(method == "Translate"){
	  	cerr << "Creating new Translate move...";
    	fullMoveStage = new MoleculeTranslate(PathData, IOSection);//, actionsToRead, startIndex);
    	switchMoveStage = new MoleculeTranslate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Rotate"){
	  	cerr << "Creating new Rotate move...";
    	fullMoveStage = new MoleculeRotate(PathData, IOSection);//, actionsToRead, startIndex);
    	switchMoveStage = new MoleculeRotate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Multi"){
	  	cerr << "Creating new Multiple move...";
    	fullMoveStage = new MoleculeMulti(PathData, IOSection);//, actionsToRead, startIndex);
    	switchMoveStage = new MoleculeMulti(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Stretch"){
	  	cerr << "Creating new bond-stretching move...";
    	fullMoveStage = new BondStretch(PathData, IOSection);//, actionsToRead, startIndex);
    	switchMoveStage = new BondStretch(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ForceBias"){
	  	cerr << "Creating new Force Bias move...";
    	fullMoveStage = new MoleculeForceBiasMove(PathData, IOSection);//, actionsToRead, startIndex);
    	switchMoveStage = new MoleculeForceBiasMove(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else {
	  	cerr << "ERROR: method " << method << " is not supported." << endl;
      exit(1);
	  }
    fullMoveStage->Read(in);
    switchMoveStage->Read(in);
    FullMoveStages.push_back(fullMoveStage);
    SwitchMoveStages.push_back(switchMoveStage);
    //startIndex += numActions(s);
    in.CloseSection();
  }

  Array<string,1> ActionList;
  assert (in.ReadVar ("FullActions", ActionList));
  int numActionsToRead = ActionList.size();
  //cerr << "  Looking for " << numActionsToRead << " actions starting at index " << startIndex << endl;
  //assert ((numActionsToRead + startIndex) <= ActionList.size());
	for(int a=0; a<numActionsToRead; a++){
    cerr << "Read in actions " << ActionList << endl;
		string setAction = ActionList(a);//+startIndex);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction);
    FullActions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;
  }
  list<MolMoveClass*>::iterator stageIter=FullMoveStages.begin();
  while (stageIter!=FullMoveStages.end()){
    (*stageIter)->LoadActions(FullActions);
    stageIter++;
  }
  assert (in.ReadVar ("SwitchActions", ActionList));
  numActionsToRead = ActionList.size();
  //cerr << "  Looking for " << numActionsToRead << " actions starting at index " << startIndex << endl;
  //assert ((numActionsToRead + startIndex) <= ActionList.size());
	for(int a=0; a<numActionsToRead; a++){
    cerr << "Read in actions " << ActionList << endl;
		string setAction = ActionList(a);//+startIndex);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction);
    SwitchActions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;
  }
  stageIter=SwitchMoveStages.begin();
  while (stageIter!=SwitchMoveStages.end()){
    (*stageIter)->LoadActions(SwitchActions);
    stageIter++;
  }

  // initialize stuff
  NumPtcls = PathData.Path.NumParticles();
  NumSlices = PathData.Path.NumTimeSlices();
  mode = FULL;
  NumMakeMove = NumFMakeMove = NumSMakeMove = 0;
  NumF2SProp = NumS2FProp = NumF2SAcc = NumS2FAcc = 0;
  NumFAttempted = NumFAccepted = NumSAttempted = NumSAccepted =0;
  cout << "NumMakeMove FullMakeMove Fattempt Faccept SwitchMakeMove Sattempt Saccept F2Sproposed F2Saccept S2Fproposed S2Faccept" << endl;
}

// need to rework this!!
void SPSClass::WriteRatio()
{
  cout << NumMakeMove << " " << NumFMakeMove << " " << NumFAttempted << " " << NumFAccepted << " " << NumSMakeMove << " " << NumSAttempted << " " << NumSAccepted << " " << NumF2SProp << " " << NumF2SAcc << " " << NumS2FProp << " " << NumS2FAcc << endl;
  //list<StageClass*>::iterator stageIter=MoveStages.begin();
  //while (stageIter!=MoveStages.end()){
  //  (*stageIter)->WriteRatio();
  //  stageIter++;
  //}
  //FinalStage->WriteRatio();
  // handle it here; don't use the default
  //MoveClass::WriteRatio();
  //RatioVar.Write(double(NumFinalAccept)/NumSteps);
}

void SPSClass::GetMode(ActionMode& setMode)
{
  setMode = mode;
}

void SPSClass::GetMode(string& setMode)
{
  if(mode == FULL)
    setMode = "FULL";
  else if(mode == SWITCH)
    setMode = "SWITCH";
}

double SPSClass::ComputeEnergy(ActionMode setMode)
{
  double TotalU  = 0.0;
  list<ActionBaseClass*>::iterator i;
  if(setMode == FULL) {
    i=FullActions.begin();
    while (i!=FullActions.end()){
      TotalU += (*i)->d_dBeta(0, NumSlices, 0);
      i++;
    }
  }
  else if(setMode == SWITCH) {
    i=SwitchActions.begin();
    while (i!=SwitchActions.end()){
      TotalU += (*i)->d_dBeta(0, NumSlices, 0);
      i++;
    }
  }
  return TotalU;
}

void SPSClass::MakeMove()
{
  NumMakeMove++;
  // propose and attempt potential switching
  if(mode == FULL) {
    NumFMakeMove++;
    // switch to log!!
    //bool trySwitch = Tf2s >= PathData.Path.Random.Local();
    bool trySwitch = Tf2s >= log(PathData.Path.Random.Local());
    if(trySwitch) {
      NumF2SProp++;
      double V = ComputeEnergy(FULL);
      double V_prime = ComputeEnergy(SWITCH);
      //cout << 0 << " " << V << " " << V_prime << " " << V_prime-V << endl;
      double deltaS = PathData.Path.tau * (V_prime - V);
      //double T = Ts2f/Tf2s;
      double logT = Ts2f - Tf2s;
      //cerr << "Attempting switch V " << V << " V_prime " << V_prime << " deltaS " << deltaS << " switch prob " << log(T)-deltaS << endl;
      double q = logT - deltaS;
      cerr << "Switching w/ log(prob) " << q << " deltaS " << deltaS << endl;
      bool makeSwitch = q >= log(PathData.Path.Random.Local());
      if(makeSwitch) {
        NumF2SAcc++;
        mode = SWITCH;
      }
    }
  }
  else if (mode == SWITCH) {
    NumSMakeMove++;
    // switch to LOg!!!! bool trySwitch = Ts2f >= log(PathData.Path.Random.Local();
    bool trySwitch = Ts2f >= log(PathData.Path.Random.Local());
    if(trySwitch) {
      NumS2FProp++;
      double V = ComputeEnergy(FULL);
      double V_prime = ComputeEnergy(SWITCH);
      //cout << 1 << " " << V << " " << V_prime << " " << V_prime-V << endl;
      double deltaS = PathData.Path.tau * (V - V_prime);
      //double T = Tf2s/Ts2f;
      double logT = Tf2s - Ts2f;
      bool makeSwitch = (logT - deltaS) >= log(PathData.Path.Random.Local());
      if(makeSwitch) {
        NumS2FAcc++;
        mode = FULL;
      }
    }
  }

  // make configurational MC move
  bool toAccept=true;
  double prevActionChange = 0.0;
  Array<int,1> ActiveParticles;
  list<MolMoveClass*>::iterator stageIter;
  if(mode == FULL) {
    NumFAttempted++;
    stageIter=FullMoveStages.begin();
    while (stageIter!=FullMoveStages.end() && toAccept){
      toAccept = (*stageIter)->Attempt(Slice1,Slice2,
	  			     ActiveParticles,prevActionChange);
      stageIter++;
    }
  }
  else if(mode == SWITCH) {
    NumSAttempted++;
    stageIter=SwitchMoveStages.begin();
    while (stageIter!=SwitchMoveStages.end() && toAccept){
      toAccept = (*stageIter)->Attempt(Slice1,Slice2,
	  			     ActiveParticles,prevActionChange);
      stageIter++;
    }
  }
  
  if (toAccept) {
    Accept();
    if(mode==FULL)
      NumFAccepted++;
    else if (mode==SWITCH)
      NumSAccepted++;
  }
  else 
    Reject();
}


//// consider switch...
//  Array<int,1> allPtcls(NumPtcls);
//  for (int p=0; p<allPtcls.size(); p++)
//    allPtcls(p) = p;
//
//	int slice =0;
//  int slice1 = slice;
//  int slice2 = slice;
//  if(NumSlices>1){
//    int P_max = NumSlices - 1;
//    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
//    slice1 = slice-1;
//    slice2 = slice+1;
//  }
//
//  double V = ComputeAction(FULL, slice1, slice2, allPtcls, 0);
//  double V_prime = ComputeAction(SWITCH, slice1, slice2, allPtcls, 0);
//  double deltaV = V - V_prime;
//  if(deltaV > VMAX) {
//    cerr << "VMAX " << VMAX << " exceeded " << deltaV << " " << V << " " << V_prime << endl;
//    assert(0);
//  }
//
//  double S = exp(deltaV - VMAX);
//  int index;
//  double T;
//  ActionMode useActionMode;
//  if(PathData.Path.Random.Local() < S) {
//    useActionMode = SWITCH;
//    T = 0.0;
//    index = 1;
//  }
//  else {
//    useActionMode = FULL;
//    // check this!! should be right
//    T = -1 * log(1.0 - S);
//    index = 0;
//  }
//
//  // based on MultiStage::MakeMove
//  double prevActionChange = T;
//  Array<int,1> activeParticles;
//  for(int i=0; i<NumSteps; i++){
//    list<StageClass*>::iterator stageIter=MoveStages.begin();
//    while (stageIter!=MoveStages.end())
//    {
//      //prevActionChange = 0.0;
//      //toAccept = (*stageIter)->Attempt(Slice1,Slice2,
//	  	//		     ActiveParticles,prevActionChange);
//      SetMode (NEWMODE);
//      double sampleRatio=(*stageIter)->Sample(slice1,slice2,activeParticles);
//      SetMode(OLDMODE);
//      double oldAction=ComputeAction(useActionMode, slice1,slice2,activeParticles,0);
//      SetMode(NEWMODE);
//      double newAction =ComputeAction(useActionMode, slice1,slice2,activeParticles,0);
//      double currActionChange=newAction-oldAction;
//      double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
//      bool toAccept = logAcceptProb>=log(PathData.Path.Random.Local()); /// Accept condition
//      if (toAccept){
//        NumAccepted++;
//        Accept();
//        (*stageIter)->Accept();
//      } else {
//        //PathData.RejectMove(Slice1,Slice2,ActiveParticles);
//        Reject();
//        (*stageIter)->Reject();
//      }
//      //NumAttempted++;
//      prevActionChange=currActionChange;
//
//      stageIter++;
//    }
//  }
//}


void SPSClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2+1,ActiveParticles);
}

void SPSClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2+1,ActiveParticles);
}
