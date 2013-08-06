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


#include "PreSampling.h"

// new read
void PreSamplingClass::Read(IOSectionClass& in)
{
  cerr << "Pre-Sampling Move read..." << endl;
  in.ReadVar("NumPreSteps",TotalNumPreSteps);
  int stages = in.CountSections("PreStage");
  assert(stages>0);
  //in.ReadVar("NumPreStages",stages);
  NumPreSampleSteps.resize(stages);
  NumPreSampleAccept.resize(stages);
  for(int n=0; n<stages; n++) {
    NumPreSampleSteps(n) = 0;
    NumPreSampleAccept(n) = 0;
  }
  for(int s=0; s<stages; s++){
    in.OpenSection("PreStage",s);
	  string method;
	  assert(in.ReadVar("MoveMethod",method));
    StageClass* MoveStage;
    //cerr << "  Init " << s+1 << " of " << stages << " stages: " << method << endl;
	  if(method == "Translate"){
	  	cerr << "Creating new Translate move...";
    	MoveStage = new MoleculeTranslate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Rotate"){
	  	cerr << "Creating new Rotate move...";
    	MoveStage = new MoleculeRotate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Multi"){
	  	cerr << "Creating new Multiple move...";
    	MoveStage = new MoleculeMulti(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Stretch"){
	  	cerr << "Creating new bond-stretching move...";
    	MoveStage = new BondStretch(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ForceBias"){
	  	cerr << "Creating new Force Bias move...";
    	MoveStage = new MoleculeForceBiasMove(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else {
	  	cerr << "ERROR: method " << method << " is not supported." << endl;
      exit(1);
	  }
    MoveStage->Read(in);
    PreStages.push_back(MoveStage);
    //startIndex += numActions(s);
    in.CloseSection();
  }

  // read in actions for presampling
  Array<string,1> getPreActions;
  assert(in.ReadVar("PreActions",getPreActions));
  int numPre = getPreActions.size();
  cerr << "PREACTION read " << numPre << endl;
	for(int a=0; a<numPre; a++){
    cerr << "Read in actions " << getPreActions << endl;
		string setAction = getPreActions(a);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction); 
    cerr << "Pushing back final action " << newAction << endl;
    PreActions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;
  }

  // read in actions for final accept/reject
  Array<string,1> getFinalActions;
  assert(in.ReadVar("FinalActions",getFinalActions));
  int numFinal = getFinalActions.size();
  cerr << "FINAL read " << numFinal << endl;
	for(int a=0; a<numFinal; a++){
    cerr << "Read in actions " << getFinalActions << endl;
		string setAction = getFinalActions(a);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction); 
    cerr << "Pushing back final action " << newAction << endl;
    FinalActions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;
  }
 // string finalMethod;
 // assert(in.ReadVar("FinalStageMethod",finalMethod));
 // int actionsToRead = numActions(numActions.size()-1);
 // if (finalMethod == "Dummy"){
 //   cerr << "Creating new dummy stage to evaluate the specified action (should be preceded by an actual move stage which is evaluate with a \"cheap\" action...";
 //   FinalStage = new PreSampleDummy(PathData, IOSection, actionsToRead, startIndex);
 //   cerr << " done." << endl;
 // } else {
 //   cerr << "Actually, I'm not supporting anything other than ``Dummy'' for the final move stage for this algorithm" << endl;
 //   assert(0);
 // }
 // FinalStage->Read(in);

  // initialize stuff
  NumPtcls = PathData.Path.NumParticles();
  NumSlices = PathData.Path.NumTimeSlices();
  Slice1 = Slice2 = 0;
  if(NumSlices>2) {
    cerr << "TROUBLE IN PRESAMPLING: MULTIPLE TIME SLICES IS NOT STABLE SINCE MULTIPE CALLS TO SAMPLE WILL RESULT IN MANY SLICES BEING MANIPULATED..." << endl;
    exit(1);
  }
  InitialPath.resize(NumSlices, NumPtcls);
  cerr << "PreSamplingClass read I have " << NumPtcls << " particles and " << NumSlices-1 << " time slices" << endl;
}

void PreSamplingClass::WriteRatio()
{
  list<StageClass*>::iterator stageIter=PreStages.begin();
  while (stageIter!=PreStages.end()){
    (*stageIter)->WriteRatio();
    stageIter++;
  }
  //FinalStage->WriteRatio();
  // handle it here; don't use the default
  //MoveClass::WriteRatio();
  RatioVar.Write(double(NumFinalAccept)/NumSteps);
}

void PrintPaths(int S, int N, PathDataClass& PD)
{
  cerr << "PrintPaths " << S << " " << N << endl;
  for(int s=0; s<S-1; s++){
    for(int n=0; n<N; n++){
      cerr << s << " " << n << " ";
      SetMode(OLDMODE);
      cerr << PD.Path(s,n) << " ";
      SetMode(NEWMODE);
      cerr << PD.Path(s,n) << endl;
    }
  }
}

void PreSamplingClass::StoreInitialPath()
{
  for(int s=0; s<NumSlices; s++){
    for(int n=0; n<NumPtcls; n++){
      InitialPath(s,n) = PathData.Path(s,n);
    }
  }
}

void PreSamplingClass::AssignInitialPath()
{
  for(int s=0; s<NumSlices; s++){
    for(int n=0; n<NumPtcls; n++){
      PathData.Path.SetPos(s,n,InitialPath(s,n));
    }
  }
}

void PreSamplingClass::MakeMove()
{
  //cerr << endl << endl << "*** MAKEMOVE *** " << endl;
  StoreInitialPath();
  //cerr << "STORED PATH IS" << endl;
  //for(int n=0; n<NumPtcls; n++)
  //  cerr << 0 << " " << n << " " << InitialPath(0,n) << endl;

  //cerr << "INITIAL PATHS OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  NumSteps++;
  NumPreAccept = 0;
  PreDeltaAction = 0.0;
  bool finalAccept, toAccept;
  double prevActionChange = 0.0;
  Array<bool,1> alreadyActive(NumPtcls);
  alreadyActive = false;
  Array<int,1> myActiveP(0);

  ActiveParticles.resize(NumPtcls);
  for(int p=0; p<ActiveParticles.size(); p++)
    ActiveParticles(p) = p;

  double oldAction=StageAction(FinalActions, Slice1,Slice2,ActiveParticles);
  //cerr << "StageAction ";
  double oldPreAction=StageAction(PreActions, Slice1,Slice2,ActiveParticles);

  for(int i=0; i<TotalNumPreSteps; i++){
    int stageIndex = 0;
    list<StageClass*>::iterator stageIter=PreStages.begin();
    while (stageIter!=PreStages.end())
    {
      NumPreSampleSteps(stageIndex)++;
      prevActionChange = 0.0;
      toAccept = (*stageIter)->Attempt(Slice1,Slice2,
	  			     ActiveParticles,prevActionChange);
      //cerr << "Paths After " << i << " " << " toAcc " << toAccept << " ";
      //PrintPaths(NumSlices, NumPtcls, PathData);
      if(toAccept){
        NumPreAccept++;
        PreDeltaAction += prevActionChange;
        //cerr << prevActionChange << " from ActiveParticles " << ActiveParticles << endl;
        Accept();
        //PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
        (*stageIter)->Accept();
        //cout << setw(12) << prevActionChange;

        for(int p=0; p<ActiveParticles.size(); p++){
          int ptcl = ActiveParticles(p);
          if(!alreadyActive(ptcl)){
            alreadyActive(ptcl) = true;
            int size = myActiveP.size();
            myActiveP.resizeAndPreserve(size+1);
            myActiveP(size) = ptcl;
          }
        }
        NumPreSampleAccept(stageIndex)++;
      } else {
        //PathData.RejectMove(Slice1,Slice2,ActiveParticles);
        Reject();
        (*stageIter)->Reject();
        //cout << setw(12) << 0.0 << " ";
      }
      //cerr << "Presample action change " << PreDeltaAction << endl;

      //cout << "now activeP is " << myActiveP << endl;
      stageIter++;
      stageIndex++;
    }
    //cerr << "PreDeltAction " << PreDeltaAction << endl;
    //double predS, finaldS;
    //ActiveParticles.resize(NumPtcls);
    //for(int p=0; p<ActiveParticles.size(); p++)
    //  ActiveParticles(p) = p;
    //cerr << i << " StageAction ";
    //predS = StageAction(PreActions, Slice1,Slice2,ActiveParticles) - oldPreAction;
    //finaldS = StageAction(FinalActions, Slice1,Slice2,ActiveParticles) - oldAction;
    //cerr << i << " Computed intermediate stage deltaS PRE " << predS << " FINAL " << finaldS << endl;

  }

  ActiveParticles.resize(NumPtcls);
  for(int p=0; p<ActiveParticles.size(); p++)
    ActiveParticles(p) = p;
  //cerr << "StageAction ";
  double newAction =StageAction(FinalActions, Slice1,Slice2,ActiveParticles);
  double newPreAction=StageAction(PreActions, Slice1,Slice2,ActiveParticles);

  //cout << "############################################################################################################" << endl;
  //cerr << "Presampled Paths are OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  //cerr << "Now resetting old path" << endl;
  //SetMode(OLDMODE);
  //AssignInitialPath();
  //cout << "FINAL PATHS OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  //ActiveParticles.resize(NumPtcls);
  //for(int p=0; p<ActiveParticles.size(); p++)
  //  ActiveParticles(p) = p;
  ActiveParticles.resize(myActiveP.size());
  for(int p=0; p<ActiveParticles.size(); p++)
    ActiveParticles(p) = myActiveP(p);

  double currActionChange=newAction-oldAction;
  prevActionChange = newPreAction - oldPreAction;
  //double logAcceptProb = -currActionChange+prevActionChange;
  double logAcceptProb = -currActionChange+PreDeltaAction;
  //perr << NumSteps << " PreDeltaS " << PreDeltaAction << " " << prevActionChange;// << endl;
  //perr << " FinalDeltaS " << currActionChange << " logAccProb " << logAcceptProb << endl;
  double logRand = log(PathData.Path.Random.Local());
  finalAccept = logAcceptProb>= logRand; /// Accept condition

  //finalAccept = FinalStage->Attempt(Slice1, Slice2, ActiveParticles, PreDeltaAction);

  if(finalAccept){
    NumFinalAccept++;
    //FinalStage->Accept();
    Accept();
  } else {
    //SetMode(NEWMODE);
    //AssignInitialPath();
    //Accept();
    //FinalStage->Reject();
    SetMode(OLDMODE);
    AssignInitialPath();
    Reject();
  }
  //NumAttempted++;
  
  if(NumSteps % 1 == 0) {
    cout << NumSteps;
    for(int n=0; n<NumPreSampleSteps.size(); n++)
      cout << " " << double(NumPreSampleAccept(n))/NumPreSampleSteps(n);
    cout << " FINAL ACCEPT ";
    cout << finalAccept;
    cout << " ratio " << setw(10) << double(NumFinalAccept)/NumSteps << endl;
  }

  //cerr << "QUITTING" << endl;
  //exit(1);
}

double PreSamplingClass::StageAction(std::list<ActionBaseClass*> ActionList, int startSlice,int endSlice, const Array<int,1> &changedParticles)
{ 
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=ActionList.begin();
  while (actionIter!=ActionList.end()){
    // here we hard-wire 0 for the level (i.e. no bisection)
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     0));
    actionIter++;
  }
  return TotalAction;
}

void PreSamplingClass::Accept()
{
  //cout << "ACCEPT SLICES " << Slice1 << " " << Slice2 << " PTCLS " << ActiveParticles << endl;
  PathData.AcceptMove(Slice1,Slice2+1,ActiveParticles);
  //int N = ActiveParticles.size();
  //Array<dVec, 2> P(NumSlices, N);
  //SetMode(NEWMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    P(s,p) = PathData.Path(s,ActiveParticles(p));
  //  }
  //}

  //SetMode(OLDMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    PathData.Path.SetPos(s,ActiveParticles(p),P(s,p));
  //  }
  //}

  //NumAccepted++;
}

void PreSamplingClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2+1,ActiveParticles);
  //int N = ActiveParticles.size();
  //Array<dVec, 2> P(NumSlices, N);
  //SetMode(OLDMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    P(s,p) = PathData.Path(s,ActiveParticles(p));
  //  }
  //}

  //SetMode(NEWMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    PathData.Path.SetPos(s,ActiveParticles(p),P(s,p));
  //  }
  //}
}

//PreSampleDummy::PreSampleDummy(PathDataClass& PathData, IOSectionClass& IO, int actionsToRead, int startIndex):
//  DummyEvaluate(PathData, IO, actionsToRead, startIndex)
//{
//  toRead = actionsToRead;
//  startI = startIndex;
//}

bool PreSampleDummy::Attempt(int &slice1, int &slice2, 
			      Array<int,1> &activeParticles,
			      double &prevActionChange)
{
  double rold, rnew;
  dVec r;
  PathData.Path.DistDisp(0, 0, 1, rold, r);
  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  PathData.Path.DistDisp(0, 0, 1, rnew, r);
  SetMode(OLDMODE);
  double oldAction=StageAction(slice1,slice2,activeParticles);
  double oldPreAction=PreSampleAction(slice1,slice2,activeParticles);
  SetMode(NEWMODE);
  double newAction =StageAction(slice1,slice2,activeParticles);
  double newPreAction=PreSampleAction(slice1,slice2,activeParticles);
  //cout << "newFinal " << setw(12) << newAction;
  //cout << " oldFinal " << setw(12) << oldAction;
  //cout << " newPre " << setw(12) << newPreAction;
  //cout << " oldPre " << setw(12) << oldPreAction;
  //cout << " rold " << setw(12) << rold << " rnew " << setw(12) << rnew;
  //cout << "FINAL DELTA-ACTION " << setw(12) << newAction-oldAction;
  //cout << " PRE DELTA-ACTION " << setw(12) << newPreAction-oldPreAction;
  double currActionChange=newAction-oldAction;
  prevActionChange = newPreAction - oldPreAction;

  //cout << "PREACTINOCHANGE " << prevActionChange << endl;
  //cout << "ATTEMPT log(sample) " << log(sampleRatio) << " -currActChg " << -currActionChange << " prevActChg " << prevActionChange << endl;
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  //cout << "  so logAcceptProb is " << logAcceptProb << endl;
  //  AcceptProb=exp(logAcceptProb);
  //  OldAcceptProb=exp(log(1/sampleRatio)+currActionChange);
  //  if (AcceptProb>1.0)
  //    AcceptProb=1.0;
  //  if (OldAcceptProb>1.0)
  //    OldAcceptProb=1.0;
  double logRand = log(PathData.Path.Random.Local());
  bool toAccept = logAcceptProb>= logRand; /// Accept condition
  //cout << "  and toAccept is " << toAccept << " based on rn " << logRand << endl;
  if (toAccept){
    NumAccepted++;
  }
  NumAttempted++;
  //cout<<"Curr action change is "<<currActionChange<<endl;
  prevActionChange=currActionChange;

  return toAccept;
}

double PreSampleDummy::PreSampleAction(int startSlice,int endSlice,
				      const Array<int,1> &changedParticles)
{
  //cerr << "PreSampleAction activeP is " << changedParticles << endl;
  //cout << "PreSampleAction size of actions is " << PreActions.size() << endl;
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=PreActions.begin();
  while (actionIter!=PreActions.end()){
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     BisectionLevel));
    actionIter++;
  }
  return TotalAction;
}

void PreSampleDummy::Read (IOSectionClass &in){
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
    MoveList = PathData.Mol.MolOfType(molIndex);
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
	else{
		cerr << "Set for SINGLE particle updates." << endl;
		MoveList.resize(1);
	}
	
  cerr << "Looking for NumPreActions " << endl;
  assert(in.ReadVar("NumPreActions", numPre));
  cerr << "Looking for NumFinalActions " << endl;
  assert(in.ReadVar("NumFinalActions", numFinal));

  Array<string,1> ActionList;
  assert (in.ReadVar ("Actions", ActionList));
  cerr << "  Looking for " << toRead << " actions starting at index " << startI << endl;
  assert(toRead == numFinal);
  //assert ((toRead + startI) <= ActionList.size());
  assert ((numPre + numFinal + startI) == ActionList.size());
  cerr << "PRE looking for " << numPre << endl;
	for(int a=0; a<numPre; a++){
    cerr << "Read in PreSample Actions " << ActionList << endl;
		string setAction = ActionList(a+startI);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction);
    cerr << "Pushing back preaction " << newAction << endl;
    PreActions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;

	  //if(setAction == "MoleculeInteractions"){
		//	// read should be done in Actions now
		//	//PathData.Actions.MoleculeInteractions.Read(in);
  	//	PreActions.push_back(&PathData.Actions.MoleculeInteractions);
		//	cerr << "  Added Molecule PreActions" << endl;
		//}else if(setAction == "ST2Water"){
		//	//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
  	//	PreActions.push_back(&PathData.Actions.ST2Water);
		//	cerr << "Added ST2Water action" << endl;
		//}else if(setAction == "Kinetic"){
  	//	PreActions.push_back(&PathData.Actions.Kinetic);
		//	cerr << "Added Kinetic action" << endl;
//#ifdef USE_QMC
		//}else if(setAction == "CEIMCAction"){
		//	//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
		//	PathData.Actions.CEIMCAction.Read(in);
  	//	PreActions.push_back(&PathData.Actions.CEIMCAction);
		//	cerr << "Added CEIMC calculation of BO energy" << endl;
//#endif
		//}else if(setAction == "LongRangeCoulomb"){
  	//	PreActions.push_back(&PathData.Actions.LongRangeCoulomb);
		//	cerr << "Added long-range coulomb interaction" << endl;
		//}else if(setAction == "IonInteraction"){
  	//	PreActions.push_back(&PathData.Actions.IonInteraction);
		//	cerr << "Added intermolecular ion-ion interaction" << endl;
		//} else if(setAction == "QBoxAction"){
  	//	PreActions.push_back(&PathData.Actions.QBoxAction);
		//	cerr << "Computing action with QBox DFT code" << endl;
		//}else if(setAction == "EAM"){
  	//	PreActions.push_back(&PathData.Actions.EAM);
		//	cerr << "Added Al EAM action" << endl;
		//}else if(setAction == "PairAction"){
		//	cerr << "Added ShortRange PairAction" << endl;
    //  PreActions.push_back(&PathData.Actions.ShortRange);
		//} else
    //	cerr << "You specified " << setAction << ", which is not supported for this type of move" << endl;
	}

  cerr << "FINAL looking for " << numFinal << endl;
	for(int a=0; a<numFinal; a++){
    cerr << "Read in actions " << ActionList << endl;
		string setAction = ActionList(a+startI+numPre);
    ActionBaseClass* newAction = PathData.Actions.GetAction(setAction); 
    cerr << "PUshing back final action " << newAction << endl;
    Actions.push_back(newAction);
    cerr << "  Added action with label " << setAction << " and address " << newAction << endl;
	  
    //if(setAction == "MoleculeInteractions"){
		//	// read should be done in actions now
		//	//PathData.Actions.MoleculeInteractions.Read(in);
  	//	Actions.push_back(&PathData.Actions.MoleculeInteractions);
		//	cerr << "  Added Molecule Actions" << endl;
		//}else if(setAction == "ST2Water"){
		//	//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
  	//	Actions.push_back(&PathData.Actions.ST2Water);
		//	cerr << "Added ST2Water action" << endl;
		//}else if(setAction == "Kinetic"){
  	//	Actions.push_back(&PathData.Actions.Kinetic);
		//	cerr << "Added Kinetic action" << endl;
//#ifdef USE_QMC
		//}else if(setAction == "CEIMCAction"){
		//	//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
		//	PathData.Actions.CEIMCAction.Read(in);
  	//	Actions.push_back(&PathData.Actions.CEIMCAction);
		//	cerr << "Added CEIMC calculation of BO energy" << endl;
//#endif
		//}else if(setAction == "LongRangeCoulomb"){
  	//	Actions.push_back(&PathData.Actions.LongRangeCoulomb);
		//	cerr << "Added long-range coulomb interaction" << endl;
		//}else if(setAction == "IonInteraction"){
  	//	Actions.push_back(&PathData.Actions.IonInteraction);
		//	cerr << "Added intermolecular ion-ion interaction" << endl;
		//} else if(setAction == "QBoxAction"){
  	//	Actions.push_back(&PathData.Actions.QBoxAction);
		//	cerr << "Computing action with QBox DFT code" << endl;
		//}else if(setAction == "EAM"){
  	//	Actions.push_back(&PathData.Actions.EAM);
		//	cerr << "Added Al EAM action" << endl;
		//} else
    //	cerr << "You specified " << setAction << ", which is not supported for this type of move" << endl;
	}
}
