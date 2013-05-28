#ifndef QMC_WRAPPER_H
#define QMC_WRAPPER_H

#include "PathDataClass.h"

#ifdef USE_QMC
  #include <QMCApp/QMCInterface.h>
  #include <QMCDrivers/QMCDriver.h>
  #include "Message/Communicate.h"
  #include "Utilities/OhmmsInfo.h"
	#include <iostream>
	#include <fstream>
#endif

class QMCWrapperClass
{
	bool initialized;
	ofstream out;
	public:
#ifdef USE_QMC
	double dt;
	int walkers;
	int steps;
	int blocks;
	int chains;
	bool correlated;
	string QMCMethod;
	int EnergyDiffIndex, EnergyIndex0, EnergyIndex1;
	void DummyInit(PathDataClass& PathData);
	void QMCDummy(PathDataClass& PathData);
	QMCWrapperClass(PathDataClass& pathData){
		initialized = false;
	  dt = pathData.dt;
    walkers = pathData.walkers;
    steps = pathData.steps;
    blocks = pathData.blocks;
    correlated = pathData.correlated;
		chains = pathData.chains;
		QMCMethod = pathData.QMCMethod;
	}
#else
	QMCWrapperClass(PathDataClass& pathData){};
	void QMCDummy(PathDataClass& PathData);
	void DummyInit(PathDataClass& PathData);
#endif
};

#endif
