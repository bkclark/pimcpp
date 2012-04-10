#include "QMCWrapper.h"
#include <sstream>

void QMCWrapperClass::DummyInit(PathDataClass& PathData){
#ifdef USE_QMC
	ostringstream filename;
	filename << "QMCDummy." << PathData.MetaWorldComm.MyProc() << "." << PathData.QMCComm.MyProc() << ".out";
	out.open(filename.str().c_str());
	initialized = true;
#endif
}

void QMCWrapperClass::QMCDummy(PathDataClass& PathData){

#if USE_QMC

	// SLICE SHOULD BE SENT, RIGHT??
	int slice = 0;

	bool newmode;
	PathData.QMCComm.Broadcast(0, newmode);
	// not sure about looping over slices; not done right now
  //for(int slice=slice1; slice<slice2; slice++) left bracket
	Array<string,1> setPtclSet(PathData.ptclSet0.size());
	for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet0(s);
	if(newmode && correlated)
		for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet1(s);
	// set ion positions for either mode
	int ptclSize;
	PathData.QMCComm.Broadcast(0,ptclSize);
	Array<int,1> SpeciesIndex(ptclSize);
	Array<int,1> OffsetList(ptclSize);
	Array<Vec3,1> CoordList(ptclSize);
	PathData.QMCComm.Broadcast(0, SpeciesIndex);
	PathData.QMCComm.Broadcast(0, OffsetList);
	PathData.QMCComm.Broadcast(0, CoordList);
	for(int i=0; i<ptclSize; i++){
	  PathData.qmc->SetPtclPos(setPtclSet(SpeciesIndex(i)), OffsetList(i), CoordList(i).data());
	}

  bool isNewDriver = true;
	if(correlated){
		if(newmode){
			if(QMCMethod == "VMC"){
	  		PathData.qmc->SetVMCMultiple(dt, walkers, steps, blocks);
			} else if (QMCMethod == "RQMC"){
	  		isNewDriver = PathData.qmc->SetRQMCMultiple(dt, chains, steps, blocks);
			} else {
				cerr << "QMCMethod " << QMCMethod << " not recognized." << endl;
				assert(0);
			}

      if(isNewDriver){
	  	  PathData.qmc->process();
	
			  EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
			  EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
			  EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
			  //WeightIndex0 = PathData.qmc->qmcDriver->addObservable("WPsi0");
			  //WeightIndex1 = PathData.qmc->qmcDriver->addObservable("WPsi1");
      }
      else{
        //PathData.qmc->qmcDriver->Estimators->resetReportSettings(false);
        PathData.qmc->qmcDriver->Estimators->reset();
      }
  		PathData.qmc->execute();
		}
	}
	else {
  	PathData.qmc->SetVMC(dt, walkers, steps, blocks);
		PathData.qmc->process();
		EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
  	PathData.qmc->execute();
	}

	vector<double> Uvalues(0);
	//vector<double> DataBuff(0);

	if(correlated){
		if(newmode){
			PathData.qmc->qmcDriver->Estimators->getData(EnergyDiffIndex,Uvalues); // get energy difference
      cerr << PathData.QMCComm.MyProc() << " collected " << Uvalues.size() << " entries" << endl;
			Array<double, 1> DataBuff(Uvalues.size());
			for(int u=0; u<DataBuff.size(); u++)
				DataBuff(u) = Uvalues[u];
      cerr << PathData.QMCComm.MyProc() << " comparing " << DataBuff.size() << " with blocks " << blocks << endl;
			assert(blocks == DataBuff.size());
			//PathData.QMCComm.Send (&DataBuff, DataBuff.size(), MPI_DOUBLE, 0, 1);
			PathData.QMCComm.Send (0, DataBuff);
		}
	}
	else{
		PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Uvalues); // get energy
		Array<double, 1> DataBuff(Uvalues.size());
		for(int u=0; u<DataBuff.size(); u++)
			DataBuff(u) = Uvalues[u];
		PathData.QMCComm.Send (0, DataBuff);
	}


#endif
}
