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

#include "WormGrow.h"

void
WormGrowMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar("MaxToGrow",MaxGrowth));
  // Construct action list
  //  WormGrowStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  WormGrowStage.Read(in);
  WormGrowStage.Actions.push_back(&PathData.Actions.Kinetic);
  WormGrowStage.Actions.push_back(&PathData.Actions.ShortRange);
  WormGrowStage.Actions.push_back(&PathData.Actions.Mu);
  
  Stages.push_back(&WormGrowStage);
  
  ActiveParticles.resize(1);
}

void 
WormGrowMoveClass::MakeMove()
{
  cerr<<"Move begun"<<endl;

  Slice1=0;
  Slice2=0;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;


  if (!PathData.Path.NowOpen){ //if worm closed, then reject
    MultiStageClass::Reject();
    return;
  }

  int headSlice,headPtcl,tailSlice,tailPtcl,numEmpty,wormSize;
  PathData.WormInfo(headSlice, headPtcl,
		    tailSlice, tailPtcl,
		    numEmpty, wormSize);  
  cerr<<"Worm size is "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;

	if(wormSize==(PathData.Path.NumTimeSlices()-1)){
		count++;
		double rmag;
		dVec r;
		PathData.Path.DistDisp(headSlice,headPtcl,tailPtcl,rmag,r);
		int bin = int(rmag/(numBins*dr));
		wormBin(bin,1)++;
		if(count%writeFreq == 0){
			out.open("WormData.dat");
			out << "# r g(r)" << endl;
			for(int i=0; i<numBins; i++){
				out << wormBin(i,0) << " " << wormBin(i,1) << endl;
			}
			out.close();
		}
	}

  if (wormSize==0){
    MultiStageClass::Reject();
    cerr<<"ERROR! WORM SHOULD NEVER BE SIZE 0 BECAUSE THEN THE WORLD SHOULD BE CLOSED! BUG BUG!"<<endl;
    return;
  }
  bool moveHead = PathData.Path.Random.Local() >0.5;
  bool grow = PathData.Path.Random.Local() < 0.5;
  int changeAmount=PathData.Path.Random.LocalInt(MaxGrowth)+1;
  //  cerr<<"Change amount is "<<changeAmount<<endl;
  if (grow && numEmpty-changeAmount<PathData.Path.NumTimeSlices()){
    MultiStageClass::Reject();
    return;
  }
  
  if (!grow && wormSize-changeAmount<1){
    MultiStageClass::Reject();
    return;
  }
  WormGrowStage.MoveHead=moveHead;
  WormGrowStage.Grow=grow;
  WormGrowStage.ChangeAmount=changeAmount;

  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;

  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }


  PathData.WormInfo(headSlice, headPtcl,
		    tailSlice, tailPtcl,
		    numEmpty, wormSize);

  cerr<<"  Worm size is "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;
  //  for (int ptcl=0;ptcl<PathData.Path.Permutation.size();ptcl++){
  //    cerr<<"Worm: "<<ptcl<<" "<<PathData.Path.Permutation(ptcl)<<endl;
  //  }

//   cerr<<"Worm information: "<<headSlice<<" "<<headPtcl<<" "
//       <<tailSlice<<" "<<tailPtcl<<" "
//       <<numEmpty<<" "<<wormSize<<endl;
  assert(numEmpty>=PathData.Path.NumTimeSlices());
  assert(wormSize>0);
  

  TimesCalled++;
  if (toAccept)
    Accept();
  else 
    Reject();

  PathData.WormInfo(headSlice, headPtcl,
		    tailSlice, tailPtcl,
		    numEmpty, wormSize);
  if (wormSize==PathData.Path.NumTimeSlices()){
    cerr<<"CORRECT Worm SIZE!"<<endl;
    cerr<<"After Worm size is "<<wormSize<<" "<<numEmpty<<" "
	<<headSlice<<" "<<headPtcl<<" "
	<<tailSlice<<" "<<tailPtcl<<" "<<endl;
  }

  cerr<<"Move done"<<endl;

  

}
