#include "MoleculeMoveManager.h"

void MoleculeMoveStageManagerClass::Read(IOSectionClass &in)
{
  cerr << "MolMoveMgr read..." << endl;
	Array<string,1> methodList;
	if(in.ReadVar("MoveMethod",methodList)) {
    cerr << "ERROR: This format is deprecated.  Each move stage should have a section Stage{} now. MoveMethod should be specified as a string within the Stage section.  See MoleculeMove.input for an example" << endl;
  }
  int stages = in.CountSections("Stage");
  for (int s=0;s<stages;s++){
    in.OpenSection("Stage",s);
    string method;
	  assert(in.ReadVar("MoveMethod",method));
  //Array<int,1> numActions(stages);
  //numActions = 1;
  //assert(in.ReadVar("NumActions", numActions));
  //assert(numActions.size() == stages);
  //int startIndex = 0;
  //for(int s=0; s<stages; s++){
  //  string method = methodList(s);
  //  int actionsToRead = numActions(s);
  //  cerr << "  Init " << s+1 << " of " << stages << " stages: " << method << endl;
	  if(method == "Translate"){
	  	cerr << "Creating new Translate move...";
    	MoveStage = new MoleculeTranslate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ParticleTranslate"){
	  	cerr << "Creating new INTRAmolecular displacement move...";
    	MoveStage = new ParticleTranslate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "DimerMove"){
	  	cerr << "Creating new DimerMove to change separation of TWO molecules...";
    	MoveStage = new DimerMove(PathData, IOSection);//, actionsToRead, startIndex);
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
	  } else if (method == "TestPairAction"){
	  	cerr << "Creating new Pair Action Test move...";
    	MoveStage = new PairActionTestClass(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ForceBias"){
	  	cerr << "Creating new Force Bias move...";
    	MoveStage = new MoleculeForceBiasMove(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "AVB"){
	  	cerr << "Creating new Aggregation Volume Bias move...";
    	MoveStage = new AVBMove(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  	cerr << "See B. Chen and J.I. Siepmann, J. Phys. Chem. B 104, 8725 (2000)" << endl;
	  } else if (method == "Dummy"){
	  	cerr << "Creating new dummy stage to evaluate the specified action (should be preceded by an actual move stage which is evaluate with a \"cheap\" action...";
      assert(s>0); // meant to be the second stage after pre-rejection by a cheap action
    	MoveStage = new DummyEvaluate(PathData, IOSection);//, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else {
	  	cerr << "ERROR: method " << method << " is not supported." << endl;
	  }
    MoveStage->Read(in);
    Stages.push_back(MoveStage);
    //startIndex += numActions(s);
    in.CloseSection();
  }
}
