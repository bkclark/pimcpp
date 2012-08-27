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

#ifndef OBSERVABLE_CLASS_H
#define OBSERVABLE_CLASS_H

#include "Angular.h"

#include "Sign.h"

#include "Coupling.h"
#include "DistanceToHead.h"
#include "DynamicStructureFactor.h"
#include "Forces.h"
#include "Hexatic.h"
#include "HBond.h"
#include "JosephsonPathDump.h"
#include "ObservableCorrelation.h"
#include "ObservableRefCorrelation.h"
#include "ObservableDiffusion.h"
#include "ObservableEnergy.h"
#include "OpenOrientation.h"
#include "PairCorrelationReweighting.h"
#include "ParticleAverageLoc.h"
//#include "VacancyLocation3.h"
#include "PathDump.h"
#include "PermutationCount.h"
#include "Phik.h"
#include "Pressure.h"
#include "StructureFactor.h"     
#include "SuperfluiDrop.h"
#include "SuperfluidFraction.h"
#include "SuperfluidFractionPerLayer.h"
#include "Time.h"
#include "TimeLindenman.h"
#include "TimeHexatic.h"
//#include "VacancyLocation.h"
//#include "VacancyLocation2.h"
//#include "VacancyLocationNearby.h"
#include "Weight.h"
#include "WindingNumber.h"
// #include "VariationalPIEnergy.h"
//#include "VacancyDensity.h"
#include "PlaneDensity.h"
#include "SpecificHeat.h"
#include "SpecificHeatA.h"
/// This template class will be used to construct distributed versions
/// of many different types of observables:  scalar observables, dVec
/// observables, array observables, array of dVec observables, etc.
/// We will write one class functions which correctly manages
/// collecting observables from all processors with MPI.
// template<class T> 
// class DistributedObservableClass : public ObservableClass
// {
//   int dummy;

//   DistributedObservableClass(PathDataClass &myPathData) : ObservableClass (myPathData)
//   { /* Do nothing for now. */ }
// };








#endif
