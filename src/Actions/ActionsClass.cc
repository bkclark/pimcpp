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

#include "../Communication/Communication.h"
#include "ActionsClass.h"
#include "../PathDataClass.h"
#include "../IO/FileExpand.h"
#include "../Blitz.h"
// now include ActionBaseClass headers here, not in .h
#include "ExternalPotential.h"
#include "CummingsWaterPotential.h"
#include "DiagonalActionClass.h"
#include "ShortRangeClass.h"
#include "ShortRangeOnClass.h"
#include "ShortRangeOn_diagonal_Class.h" 
#include "ShortRangeOn_diagonal_displace_Class.h" 
#include "ShortRangeApproximateClass.h"
#include "ShortRangePrimitive.h"
#include "LongRangeClass.h"
#include "LongRangeCoulombClass.h"
#include "LongRangeRPAClass.h"
#include "ShortRangePotClass.h"
#include "LongRangePotClass.h"
#include "KineticClass.h"
#include "MoleculeInteractionsClass.h"
#include "ST2WaterClass.h"
#include "TIP5PWaterClass.h"
#include "EAMClass.h"
#include "NodalActionClass.h"
#include "FreeNodalActionClass.h"
#include "SHONodalActionClass.h"
#include "Sal.h"
// #include "GroundStateNodalActionClass.h"

#include "DavidLongRangeClassYk.h"
#include "DavidLongRangeClass.h"
#include "QMCSamplingClass.h"
#include "OpenLoopImportance.h"

#include "StructureReject.h"
// #include "KineticRotorClass.h"
// #include "KineticSphereClass.h"
// #include "KineticVibrationClass.h"



#include "Mu.h"
#include "Tether.h"
// #include "NonlocalClass.h"
//#include "ReadAction.h"
#include "BlendActions.h"

///Actionsclass. Stores all the actsion

void ActionsClass::Read(IOSectionClass &in)
{
  //cerr<<"Reading Action."<<endl;
  //verr<<"Starting Actions Read"<<endl;
  PathClass &Path = PathData.Path;

  // Reading Pair Action Files
  ReadPairActions(in);

  assert(in.ReadVar ("MaxLevels", MaxLevels));
  assert(in.ReadVar ("NumImages", NumImages));
  Kinetic.SetNumImages (NumImages);
  //   KineticSphere.SetNumImages(NumImages);
  Mu.Read(in);
  UseNonlocal = false;

  //bool checkJosephson=false;
  //in.ReadVar("Josephson",checkJosephson);
  //if (checkJosephson){
  //  Josephson.Read(in);
  //  Hermele.Read(in);
  //  DualHermele.Read(in);
  //}

  //bool usePairAction=false;
  //in.ReadVar("Paired",usePairAction);
  //if (usePairAction){
  //  PairFixedPhase.Read(in);
  //}

  Array<string,1> SpecificHeatPAFiles;
  bool readSpecificHeatFiles=false;
  if (in.ReadVar("SpecificHeatPAFiles",SpecificHeatPAFiles))
    readSpecificHeatFiles=true;
  if (readSpecificHeatFiles)
    SpecificHeatPairArray.resize(SpecificHeatPAFiles.size());

  // Read pair actions files
  verr << "declaring IOSectionClass...";
  IOSectionClass PAIO;


  if (readSpecificHeatFiles){
    cerr<<"I READ SPECIFIC HEAT FILES"<<endl;
    assert((in.ReadVar("TauValues",TauValues)));
    assert(TauValues.size()==SpecificHeatPAFiles.size());
    for (int i=0; i<SpecificHeatPAFiles.size() ; i++) {
      // Allow for tilde-expansion in these files
      string name = ExpandFileName(SpecificHeatPAFiles(i));
      assert(PAIO.OpenFile (name));
      SpecificHeatPairArray(i) = ReadPAFit (PAIO, TauValues(i), MaxLevels);
      PAIO.CloseFile();
    }
  }

  //cerr<<"CHECKING LONG RANGE"<<endl;
  //in.ReadVar("UseLongRange", UseLongRange);
  //if (HaveLongRange()) {
  //  perr << "*** Using long-range/short-range breakup. ***\n";
  //  assert (in.ReadVar("UseBackground", LongRange.UseBackground));
  //  LongRangePot.UseBackground = LongRange.UseBackground;
  //  LongRangeRPA.UseBackground = LongRange.UseBackground;
  //}

//   if (longRange){
//     LongRange.Init(in);
//     if (UseRPA)
//       LongRangeRPA.Init(in);
//   }

  //cerr << "About to Read Nodal Action." << endl;
  ReadNodalActions (in);

  if (in.OpenSection("StructureReject")) {
    StructureReject.Read(in);
    in.CloseSection();
  }
  if (in.OpenSection("Tether")) {
    Tether.Read(in);
    in.CloseSection();
  }

  if (!in.ReadVar ("UseRPA", UseRPA))
    UseRPA = false;
  if (UseRPA)
    verr << "Using RPA for long range action.\n";
  else
    verr << "Not using RPA for long range action.\n";

  ///Reading in information for David long range action
  if (PathData.Path.DavidLongRange)
    DavidLongRange.ReadYk();

  if (PathData.Path.OpenPaths)
    OpenLoopImportance.Read(in);

  // new read section
  assert(in.ReadVar("MaxLevels",MaxLevels));
  assert(in.ReadVar("NumImages",NumImages));

  ActionBaseClass* newAction;

  int numActions = in.CountSections("Action");
  vector<string> L(0);
  //cerr << "Initializing " << numActions << " action objects" << endl;
  for(int n=0; n<numActions; n++){
    in.OpenSection("Action",n);
    string type, label;
    assert(in.ReadVar("Type", type));
    assert(in.ReadVar("Name", label));
    //cerr << "  Initializing action " << label << " of type " << type;// << endl;
    // add menu of all possible ActionBase objects here
    if(type == "Kinetic"){
      newAction = new KineticClass(PathData);
    } else if (type == "ShortRange") {
      newAction = new ShortRangeClass(PathData, PairMatrix);
    } else if (type == "DiagonalAction") {
      newAction = new DiagonalActionClass(PathData,PairMatrix);
    } else if (type == "DiagonalActionOrderN") {
      newAction = new ShortRangeOn_diagonal_class(PathData,PairMatrix);
    } else if (type == "DiagonalDisplaceActionOrderN") {
      newAction = new ShortRangeOn_diagonal_displace_class(PathData,PairMatrix);
    } else if (type == "ShortRangeOrderN") {
      newAction = new ShortRangeOnClass(PathData,PairMatrix);
    } else if (type == "Sal") {
      newAction = new SalClass(PathData);
    } else if (type == "MoleculeInteractions") {
      newAction = new MoleculeInteractionsClass(PathData);
    } else if (type == "CummingsWater") {
      newAction = new CummingsWaterPotentialClass(PathData);
    } else if (type == "ST2Water") {
      newAction = new ST2WaterClass(PathData);
    } else if (type == "EAM") {
      newAction = new EAMPotentialClass(PathData);
#ifdef USE_QBOX
    } else if (type == "QboxAction") {
      newAction = new QBoxActionClass(PathData);
#endif
    // } else if (type == "FixedAxisRotor") {
    //   newAction = new FixedAxisRotorClass(PathData);
    // } else if (type == "KineticRotor") {
    //   newAction = new KineticRotorClass(PathData);
    // } else if (type == "KineticVibration") {
    //   newAction = new KineticVibrationClass(PathData);
    // } else if (type == "KineticVibrationEigenFunction") {
    //   newAction = new KineticVibrationEigenFunctionClass(PathData);
    } else if (type == "LongRangeCoulomb") {
      newAction = new LongRangeCoulombClass(PathData, PairMatrix, PairArray);
    // } else if (type == "ReadFromFile") {
    //   newAction = new ReadFromFileActionClass(PathData);
    } else if (type == "BlendActions") {
      newAction = new BlendActionsClass(PathData);
    } else if (type == "ExternalPotential") {
      newAction = new ExternalPotential(PathData);
    } else if (type == "Water") {
      newAction = new WaterClass(PathData);
    } else {
      cerr << endl << "ActionBaseClass Type " << type << " not recognized" << endl;
      exit(0);
    }

    cout << PathData.Path.CloneStr << " Added action of type " << type << endl;
    newAction->Read(in);
    //cerr << " with address " << newAction << endl;
    ActionList.push_back(newAction);
    ActionLabels.push_back(label);
    // check for duplicate labels
    L.push_back(label);
    for(int l=0; l<L.size()-1; l++) {
      if(label == L[l]) {
        cerr << "Duplicate label " << label << " conflicts with " << L[l] << " at " << l << endl;
        assert(0);
      }
    }

    in.CloseSection();
  }
}

void 
ActionsClass::ReadPairActions(IOSectionClass &in)
{
  PathClass &Path=PathData.Path;
  //cerr << "Starting to read pair action files" << endl;
  Array<string,1> PAFiles;
  assert (in.ReadVar ("PairActionFiles", PAFiles));
  int numPairActions = PAFiles.size();
  PairArray.resize(numPairActions);
  PairMatrix.resize(Path.NumSpecies(),Path.NumSpecies());
  PairIndex.resize(Path.NumSpecies(),Path.NumSpecies());
  // Initialize to a nonsense value so we can later check in the table
  // element was filled in.
  for (int i=0; i<Path.NumSpecies(); i++)
    for (int j=0; j<Path.NumSpecies(); j++)
      PairMatrix(i,j) = (PairActionFitClass*)NULL;
  // Read pair actions files
  IOSectionClass PAIO;
  for (int i=0; i<numPairActions; i++) {
    // Allow for tilde-expansion in these files
    string name = ExpandFileName(PAFiles(i));
    //cerr<<"We are about to try to read the PairActon file "<<name<<endl;
    //cerr<<"If the code aborts at this point, most likely the .PairAction or .Sampling.in or .dm file are incorrect or placed in the wrong place."<<endl;
    if (!PAIO.OpenFile(name)){
      cerr<<"We were unable to find the PairAction file "<<name<<endl;
      cerr<<"Please make sure you have specified it correctly."<<endl;
      abort();
    }
    PairArray(i) = ReadPAFit (PAIO, Path.tau, MaxLevels);
    bool paUsed=false;
    for (int spec1=0;spec1<Path.NumSpecies();spec1++)
      for (int spec2=spec1;spec2<Path.NumSpecies();spec2++) 
        if (((Path.Species(spec1).Type==PairArray(i)->Particle1.Name)&&
             (Path.Species(spec2).Type==PairArray(i)->Particle2.Name)) ||
            ((Path.Species(spec2).Type==PairArray(i)->Particle1.Name)&&
             (Path.Species(spec1).Type==PairArray(i)->Particle2.Name))) {
          if (PairMatrix(spec1,spec2) != NULL) {
            perr << "More than one pair action for species types (" 
                 << PairArray(i)->Particle1.Name << ", "
                 << PairArray(i)->Particle2.Name << ")." << endl;
            exit(-1);
          }
          cout << Path.CloneStr << " Found PAfile for pair (" 
               << Path.Species(spec1).Name << ", "
               << Path.Species(spec2).Name << ")\n";
          PairMatrix(spec1,spec2) = PairArray(i);
          PairMatrix(spec2,spec1) = PairArray(i);
          PairIndex(spec2,spec1) = i;
          PairIndex(spec1,spec2) = i;
          paUsed = true;
        }
    if (!paUsed) {
      perr << "Warning:  Pair action for species types (" 
           << PairArray(i)->Particle1.Name << ", "
           << PairArray(i)->Particle1.Name << ") not used.\n";
    }
    PAIO.CloseFile();
  }
  //cerr<<"Done reading PairAction files"<<endl;
  //cerr << "READPAIRACTION PA INIT" << endl;
  //for (int pai=0; pai<PairArray.size(); pai++) {
  //  cout << PairArray(pai) << " " << PairArray(pai)->Particle1.Name << " " << PairArray(pai)->Particle2.Name << endl;
  //  for(int i=0; i<100; i++) {
  //    double r = i*2./30;
        //    double U;
        //    U = PairArray(pai)->U(r,0,0, 0);
  //    cout << r << " " << U << endl;
  //  }
  //}
}

ActionBaseClass* ActionsClass::GetAction(string name)
{
  std::list<ActionBaseClass*>::iterator actionIt = ActionList.begin();
  std::list<string>::iterator labelIt = ActionLabels.begin();

  do {
    if(name == *labelIt){
      return(*actionIt);
    }
    actionIt++;
    labelIt++;
  } while(labelIt != ActionLabels.end());

  if(name == "") {
    cerr << "ActionsClass::GetAction WARNING returning NULL action because you specified an empty string.  Are you sure this was intended?" << endl;
    return NULL;
  }
  cerr << "ActionsClass::GetAction ERROR:" << endl;
  cerr << "Requested ActionBaseClass object with label " << name << " was not found.  Make sure it's specified in the Actions section of the input file." << endl;
  exit(1);
  return NULL;
}


/// Read in the nodal actions.
/// This should only be called after the PairActions have been read.
void
ActionsClass::ReadNodalActions(IOSectionClass &in)
{
  //std::cerr << "Reading Nodal Action." << endl;
  int numNodeSections=in.CountSections("NodalAction");
  cout << PathData.Path.CloneStr << " Found " << numNodeSections << " Nodal Actions." << endl;
  NodalActions.resize (PathData.Path.NumSpecies());
  NodalActions = NULL;
  for (int nodeSection=0; nodeSection<numNodeSections; nodeSection++) {
    in.OpenSection("NodalAction", nodeSection);
    string type, speciesString;
    assert (in.ReadVar ("Type", type));
    if (type == "FREE") {
      assert (in.ReadVar("Species", speciesString));
      int species = PathData.Path.SpeciesNum(speciesString);
      FreeNodalActionClass &nodeAction = 
        *(new FreeNodalActionClass (PathData, species));
      nodeAction.Read(in);
      NodalActions(species) = &nodeAction;
    }
    else if (type == "GROUNDSTATE") {
//       GroundStateClass &groundState = *new GroundStateClass(PathData);
//       groundState.Read (in);
//       NodalActions(groundState.UpSpeciesNum) = 
//      new GroundStateNodalActionClass 
//      (PathData, groundState, groundState.UpSpeciesNum);
//       NodalActions(groundState.DownSpeciesNum) = 
//      new GroundStateNodalActionClass 
//      (PathData, groundState, groundState.DownSpeciesNum);
//       NodalActions(groundState.IonSpeciesNum) = 
//      new GroundStateNodalActionClass 
//      (PathData, groundState, groundState.IonSpeciesNum);
    }
    else if (type == "FIXEDPHASE") {
//       FixedPhaseA = new FixedPhaseClass(PathData);
//       FixedPhaseA->Read (in);
//       if (PathData.Path.UseCorrelatedSampling()) {
//      FixedPhaseB = new FixedPhaseClass(PathData);
//        FixedPhaseB->Read (in);
//    }  else{
      //        FixedPhaseB =  FixedPhaseA;
      //    }
      // Now setup up actual actions
//       NodalActions(FixedPhaseA->UpSpeciesNum)   = new FixedPhaseActionClass 
//      (PathData, *FixedPhaseA, *FixedPhaseB, FixedPhaseA->UpSpeciesNum);
//       NodalActions(FixedPhaseA->DownSpeciesNum) = new FixedPhaseActionClass 
//      (PathData, *FixedPhaseA, *FixedPhaseB, FixedPhaseA->DownSpeciesNum);
//       NodalActions(FixedPhaseA->IonSpeciesNum) = new FixedPhaseActionClass 
//      (PathData, *FixedPhaseA, *FixedPhaseB, FixedPhaseA->IonSpeciesNum);
    }
    else if (type == "SHO") {
      //cerr << "SHO Nodal Action." << endl;
      assert (in.ReadVar("Species", speciesString));
      int species = PathData.Path.SpeciesNum(speciesString);
      //cerr << "Action Species: " << species << endl;
// AGGRESSIVE COMPILING ERROR (FIX)
      SHONodalActionClass *nodeAction = (new SHONodalActionClass (PathData, species));
      nodeAction -> Read(in);
      NodalActions(species) = nodeAction;
    }
    in.CloseSection();
  }
  verr<<"Ending actions read"<<endl;
}



void
ActionsClass::Energy (double& kinetic, double &dUShort, double &dULong, 
                      double &node, double &vShort, double &vLong,
                      double &duNonlocal)
{
  double residual;
  Energy(kinetic,dUShort,dULong,node,vShort,vLong,duNonlocal,residual);
}

void
ActionsClass::Energy (double& kinetic, double &dUShort, double &dULong, 
                      double &node, double &vShort, double &vLong,
                      double &duNonlocal,
                      double &residual)
//void ActionsClass::Energy(map<double>& Energies)
{
        //double kinetic, dUShort, dULong, node, vShort, vLong;
  bool doLongRange = HaveLongRange() && UseLongRange;
  int M = PathData.Path.NumTimeSlices()-1;
  kinetic = Kinetic.d_dBeta (0, M, 0);
  if (PathData.Path.OrderN){

    //dUShort=0.0;
    ///    dUShort=ShortRangeOn.d_dBeta(0,M,0);
    dUShort=((ShortRangeOn_diagonal_class*)(GetAction("DiagonalActionOrderN")))->d_dBeta(0,M,0);

    //    cerr<<"Energies: "<<dUShort<<" "<<dUShortp<<endl;
    //residual=((ShortRangeOn_diagonal_class*)(GetAction("DiagonalActionOrderN")))->residual_energy();
  }
  else
    dUShort = ShortRange.d_dBeta (0, M, 0);
  dULong=0.0;
  if (doLongRange){
    if (UseRPA)
      dULong = LongRangeRPA.d_dBeta (0, M, 0);
    else
      dULong = LongRange.d_dBeta (0, M, 0);
  }
  vShort = 0.0; vLong = 0.0;
  if (PathData.Path.DavidLongRange){
    dULong = DavidLongRange.d_dBeta(0,M,0);
    vLong = DavidLongRange.V(0,M,0);
  }
  node = 0.0;
  for (int species=0; species<PathData.Path.NumSpecies(); species++)
    if (NodalActions(species) != NULL)
      node += NodalActions(species)->d_dBeta(0, M, 0);


//   ///HACK! BUG! Removing V  
//   for (int slice=0; slice <= M; slice++) {
//     double factor = ((slice==0)||(slice==M)) ? 0.5 : 1.0;
//     vShort += factor * ShortRangePot.V(slice);
//     if (doLongRange)
//       vLong  += factor *  LongRangePot.V(slice);
 
//   }
//   if (UseNonlocal)
//     duNonlocal = Nonlocal.d_dBeta(0,M,0);
//   else
    duNonlocal = 0.0;
  //Energies["kinetic"] = kinetic;
  //Energies["dUShort"] = dUShort;
  //Energies["dULong"] = dULong;
  //Energies["node"] = node;
  //Energies["vShort"] = vShort;
  //Energies["vLong"] = vLong;
}



void
ActionsClass::GetActions (double& kinetic, double &UShort, double &ULong, 
                          double &node)
{
  bool doLongRange = HaveLongRange() && UseLongRange;
  Array<int,1> activePtcls(PathData.Path.NumParticles());
  for (int i=0; i<PathData.Path.NumParticles(); i++)
    activePtcls(i) = i;
  
  int M = PathData.Path.NumTimeSlices()-1;
  kinetic = Kinetic.Action (0, M, activePtcls, 0);
  UShort = ShortRange.Action (0, M, activePtcls, 0);
  ULong=0.0;
  if (doLongRange){
    if (UseRPA)
      ULong = LongRangeRPA.Action (0, M, activePtcls, 0);
    else if (PathData.Path.DavidLongRange)
      ULong = DavidLongRange.Action(0,M, activePtcls, 0);
    else
      ULong = LongRange.Action (0, M, activePtcls, 0);
  }
  node = 0.0;
  for (int species=0; species<PathData.Path.NumSpecies(); species++)
    if (NodalActions(species) != NULL)
      node += NodalActions(species)->Action(0, M, activePtcls, 0);
}


//   PathClass &Path = PathData.Path;
//   for (int link=0; link<Path.NumTimeSlices()-1; link++) {    
//     for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
//       int specNum1 = Path.ParticleSpeciesNum(ptcl1);
//       SpeciesClass &spec1 = Path.Species(specNum1);
//       if (spec1.lambda != 0.0) {
//      // Compute kinetic energy
//      /// Add constant part to kinetic part of energy
//      kinetic += (NDIM*0.5)/tau;
//      // Now do spring part
//      double fourLambdaTauInv = 1.0/(4.0*spec1.lambda*tau);
//      dVec vel;
//      vel = Path.Velocity (link, link+1, ptcl1);
//      double Z = 1.0;
//      dVec gaussSum = 0.0;
//      dVec numSum = 0.0;
//      for (int dim=0; dim<NDIM; dim++) {
//        for (int image=-NumImages; image<=NumImages; image++) {
//          double dist = vel[dim]+(double)image*Path.GetBox()[dim];
//          double dist2OverFLT = dist*dist*fourLambdaTauInv;
//          double expPart = exp(-dist2OverFLT);
//          gaussSum[dim] += expPart;
//          numSum[dim]   += dist2OverFLT*expPart/tau;
//        }
//        Z *= gaussSum[dim];
//      }
      
//      double scalarnumSum = 0.0;
//      for (int dim=0;dim<NDIM;dim++){
//        dVec numProd=1.0;
//        for (int dim2=0;dim2<NDIM;dim2++)
//          if (dim2!=dim)
//            numProd[dim] *= gaussSum[dim2];
//          else 
//            numProd[dim] *=  numSum[dim2];
//        scalarnumSum += numProd[dim];
//      }
//      kinetic += scalarnumSum/Z; 
//       }
    
//       // Now do short-range part of energy
//       for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
//      int specNum2 = Path.ParticleSpeciesNum(ptcl2);
//      dVec r, rp;
//      double rmag, rpmag;
//      Path.DistDisp(link, link+1, ptcl1, ptcl2, rmag,rpmag,r,rp); 
        
//      double s2 = dot(r-rp, r-rp);
//      double q = 0.5*(rmag+rpmag);
//      double z = (rmag-rpmag);
        
//      PairActionFitClass &pa = *PairMatrix(specNum1,specNum2);
//      duShort += pa.dU(q, z, s2, 0);
//      // Subtract off long-range part from short-range action
//      if (pa.IsLongRange())
//        duShort -= 0.5*(pa.dUlong(0)(rmag)+pa.dUlong(0)(rpmag));
//       }
//     }
//   }

//   if (UseRPA)
//     dULong = LongRangeRPA.d_dBeta (0, Path.NumTimeSlices()-1, 0);
//   else
//     dULong = LongRange.d_dBeta (0, Path.NumTimeSlices()-1, 0);


//    // Now, calculate potential
//    for (int slice=0; slice < Path.NumTimeSlices; slice++) {
//      double factor = ((slice==0)||(slice==Path.NumTimesSlices()-1)) ? 0.5 : 1.0;
//      for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
//        int specNum1 = Path.ParticleSpeciesNum(ptcl1);
//        for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
//       int specNum2 = Path.ParticleSpeciesNum(ptcl2);
//       double dist;
//       dVec disp;
//       Path.DistDisp (slice, ptcl1, ptcl2, dist, disp);
//       PairActionFitClass &pa = *PairMatrix(specNum1,specNum2);
//       vShort += factor * pa.V(dist);
//       if (pa.IsLongRange())
//         vShort -= factor * pa.Vlong(dist);
//        }
//      }
//    }

     
         

// }


Potential&
ActionsClass::GetPotential (int species1, int species2)
{
  return *(PairMatrix(species1, species2)->Pot);
}


void 
ActionsClass::ShiftData (int slicesToShift)
{
  OpenLoopImportance.ShiftData(slicesToShift);
  StructureReject.ShiftData(slicesToShift);
  ShortRange.ShiftData(slicesToShift);
  ShortRangeOn.ShiftData(slicesToShift);
  ShortRangeOnDiagonal.ShiftData(slicesToShift);
  ShortRangeApproximate.ShiftData(slicesToShift);
  ShortRangePrimitive.ShiftData(slicesToShift);
  DiagonalAction.ShiftData(slicesToShift);
  LongRange.ShiftData(slicesToShift);
  LongRangeRPA.ShiftData(slicesToShift);
  DavidLongRange.ShiftData(slicesToShift);
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i)!=NULL)
      NodalActions(i)->ShiftData(slicesToShift);
}


void ActionsClass::AcceptCopy (int startSlice, int endSlice, const Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;
  bool activeSpecies[Path.NumSpecies()];
  for (int i=0; i<Path.NumSpecies(); i++)
    activeSpecies[i] = false;
  for (int pi=0; pi<activeParticles.size(); pi++)
    activeSpecies[Path.ParticleSpeciesNum(activeParticles(pi))] = true;

  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL && activeSpecies[i])
      NodalActions(i)->AcceptCopy (startSlice, endSlice);

  //QBoxAction.AcceptCopy(startSlice, endSlice);
}


void ActionsClass::RejectCopy (int startSlice, int endSlice, const Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;
  bool activeSpecies[Path.NumSpecies()];
  for (int i=0; i<Path.NumSpecies(); i++)
    activeSpecies[i] = false;
  for (int pi=0; pi<activeParticles.size(); pi++)
    activeSpecies[Path.ParticleSpeciesNum(activeParticles(pi))] = true;

  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL && activeSpecies[i])
      NodalActions(i)->RejectCopy (startSlice, endSlice);

  //QBoxAction.RejectCopy(startSlice, endSlice);
}

void
ActionsClass::Init()
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->Init();
}


bool ActionsClass::HaveLongRange()
{
  bool longRange = false;
  for (int i=0; i<PairArray.size(); i++)
    longRange = longRange || PairArray(i)->IsLongRange();
  return (UseLongRange && longRange);
}

void 
ActionsClass::Setk(Vec3 k)
{
  for (int i=0; i<PathData.Path.NumSpecies(); i++) {
    if (NodalActions(i) != NULL)
      NodalActions(i)->Setk(k);
  }
}

void
ActionsClass::WriteInfo(IOSectionClass &out)
{
  /// If we have nodal actions, have them write any pertinent info to
  /// the output file.
  bool haveNodeActions = false;
  for (int i=0; i<PathData.Path.NumSpecies(); i++) 
    if (NodalActions(i) != NULL)
      haveNodeActions = true;
  if (haveNodeActions) {
    if (PathData.IntraComm.MyProc() == 0) {
      out.NewSection("NodalActions");
      for (int i=0; i<PathData.Path.NumSpecies(); i++) {
        out.NewSection ("NodeAction");
        if (NodalActions(i) == NULL) 
          out.WriteVar("Type", "NONE");
        else 
          NodalActions(i)->WriteInfo(out);
        
        out.CloseSection();
      }
      out.CloseSection();
    }
  }
}


void
ActionsClass::GetForces(const Array<int,1> &ptcls, 
                        Array<dVec,1> &Fshort, Array<dVec,1> &Flong)
{
  //Move the join to the end so we don't have to worry about
  //permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  PathClass &Path = PathData.Path;
  assert (Fshort.size() == ptcls.size());
  assert (Flong.size() == ptcls.size());
  Array<dVec,1> FtmpShort(Fshort.size());
  Array<dVec,1> FtmpLong(Flong.size());
  dVec zero; 
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  FtmpShort = zero;
  FtmpLong = zero;
  /// Calculate short and long-range pair actions for now.
  ShortRange.GradAction(0, Path.NumTimeSlices()-1, ptcls, 0, FtmpShort);
  LongRange.GradAction(0, Path.NumTimeSlices()-1, ptcls, 0, FtmpLong);
  /// Answer must be divided by beta.
  double beta = Path.TotalNumSlices * Path.tau;
  //  Fshort -= (1.0/beta)*FtmpShort;
  //  Flong -= (1.0/beta)*FtmpLong;
  for (int dim=0;dim<NDIM;dim++){
    Fshort[dim]-=(1.0/beta)*FtmpShort[dim];
    Flong[dim]-=(1.0/beta)*FtmpLong[dim];
    
  }
}


void
ActionsClass::GetForcesFD(const Array<int,1> &ptcls, Array<dVec,1> &F)
{
  const double eps = 1.0e-6;
  PathClass &Path = PathData.Path;
  Array<dVec,1> Ftmp(F.size());
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Ftmp = zero;
  Array<int,1> onePtcl(1);
  for (int pi=0; pi < ptcls.size(); pi++) {
    int ptcl = ptcls(pi);
    onePtcl(0) = ptcl;
    dVec savePos = Path(0, ptcl);
    dVec uPlus, uMinus;
    for (int dim=0; dim<NDIM; dim++) {
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
        Path(slice,ptcl)[dim] = savePos[dim] + eps;
      uPlus[dim] = ShortRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      uPlus[dim] += LongRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
        Path(slice,ptcl)[dim] = savePos[dim] - eps;
      uMinus[dim] = ShortRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      uMinus[dim] += LongRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
        Path(slice,ptcl)[dim] = savePos[dim];
    }
    Ftmp(pi) = (0.5/eps)*(uPlus-uMinus);
  }
  double beta = Path.TotalNumSlices * Path.tau;
  for (int dim=0;dim<NDIM;dim++){
    F[dim]-=(1.0/beta)*Ftmp[dim];
  }
  //  F -= (1.0/beta)*Ftmp;
}


void
ActionsClass::UpdateNodalActions()
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->Update();
}

void
ActionsClass::MoveJoin (int oldJoinPos, int newJoinPos)
{
  // Currently, only some nodal actions actually need their MoveJoin called.
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->MoveJoin (oldJoinPos, newJoinPos);
}

// old deprecated read
//void 
//ActionsClass::Read(IOSectionClass &in)
//{ 
//  verr<<"Starting Actions Read"<<endl;
//  PathClass &Path = PathData.Path;
//
//  assert(in.ReadVar ("MaxLevels", MaxLevels));
//  assert(in.ReadVar ("NumImages", NumImages));
//  Kinetic.SetNumImages (NumImages);
//  KineticSphere.SetNumImages(NumImages);
//  Mu.Read(in);
//      bool doMolRead = false;
//      in.ReadVar("InitMoleculeInteractions",doMolRead);
//  if(doMolRead){
//              MoleculeInteractions.Read(in);
//              MoleculeInteractions.SetNumImages(NumImages);
//      }
//
//  bool doLRCRead = false;
//      in.ReadVar("InitLongRangeCoulomb",doLRCRead);
//  if(doLRCRead){
//              LongRangeCoulomb.Read(in);
//      }
//      bool doQBoxRead = false;
//      in.ReadVar("InitQBoxAction",doQBoxRead);
//  if(doQBoxRead){
//              QBoxAction.Read(in);
//      }
//
//      bool doEAMRead = false;
//      in.ReadVar("InitEAMAction",doEAMRead);
//  if(doEAMRead){
//              EAM.Read(in);
//      }
//
//      bool doKineticRotorRead = false;
//      in.ReadVar("InitKineticRotorAction",doKineticRotorRead);
//  if(doKineticRotorRead){
//              KineticRotor.Read(in);
//      }
//
//      bool doFARRead = false;
//      in.ReadVar("InitFixedAxisRotorAction",doFARRead);
//  if(doFARRead){
//              FixedAxisRotor.Read(in);
//      }
//#ifdef ORDER_N_FERMIONS
//  VariationalPI.Read(in);
//  TruncatedInverse.Read(in);
//#endif
//  //  VariationalPI.Read(in);
//  verr << "MaxLevels = " << MaxLevels << endl;
//  bool checkJosephson=false;
//  in.ReadVar("Josephson",checkJosephson);
//  if (checkJosephson){
//    Josephson.Read(in);
//    Hermele.Read(in);
//    DualHermele.Read(in);
//  }
//  bool usePairAction=false;
//  in.ReadVar("Paired",usePairAction);
//  if (usePairAction){
//    PairFixedPhase.Read(in);
//  }
//  
//  if (!in.ReadVar ("UseRPA", UseRPA))
//    UseRPA = false;
//  if (UseRPA) 
//    verr << "Using RPA for long range action.\n";
//  else      
//    verr << "Not using RPA for long range action.\n";
//
//  verr << " going to read pair actions" << endl;
//  Array<string,1> PAFiles;
//  assert (in.ReadVar ("PairActionFiles", PAFiles));
//  Array<string,1> SpecificHeatPAFiles;
//  bool readSpecificHeatFiles=false;
//  if (in.ReadVar("SpecificHeatPAFiles",SpecificHeatPAFiles))
//    readSpecificHeatFiles=true;
//  
//
//  int numPairActions = PAFiles.size();
//  ///  cerr << "Looking for " << numPairActions << endl;
//  PairArray.resize(numPairActions);
//  if (readSpecificHeatFiles)
//    SpecificHeatPairArray.resize(SpecificHeatPAFiles.size());
//  PairMatrix.resize(Path.NumSpecies(),Path.NumSpecies());
//  // Initialize to a nonsense value so we can later check in the table
//  // element was filled in.
//  for (int i=0; i<Path.NumSpecies(); i++)
//    for (int j=0; j<Path.NumSpecies(); j++)
//      PairMatrix(i,j) = (PairActionFitClass*)NULL;
//  // Read pair actions files
//  verr << "declaring IOSectionClass...";
//  IOSectionClass PAIO;
//  cerr << " done" << endl;
// 
//  UseNonlocal = false;
//  for (int i=0; i<numPairActions; i++) {
//    // Allow for tilde-expansion in these files
//    string name = ExpandFileName(PAFiles(i));
//    verr<<"We are about to try to read the PairActon file "<<name<<endl;
//    verr<<"If the code aborts at this point, most likely the .PairAction or .Sampling.in or .dm file are incorrect or placed in the wrong place."<<endl;
//    if (!PAIO.OpenFile(name)){
//      cerr<<"We were unable to find the PairAction file "<<name<<endl;
//      cerr<<"Please make sure you have specified it correctly."<<endl;
//      abort();
//    }
//    verr<<"Done reading the PairAction file "<<name<<endl;
//    //    assert(PAIO.OpenFile (name));
//    verr << i << ": reading " << name << endl;
//    PairArray(i) = ReadPAFit (PAIO, Path.tau, MaxLevels);
//    if (PairArray(i)->Pot != NULL)
//      if (PairArray(i)->Pot->IsNonlocal())
//      UseNonlocal = true;
//    bool paUsed=false;
//    for (int spec1=0;spec1<Path.NumSpecies();spec1++)
//      for (int spec2=spec1;spec2<Path.NumSpecies();spec2++) 
//      if (((Path.Species(spec1).Type==PairArray(i)->Particle1.Name)&&
//           (Path.Species(spec2).Type==PairArray(i)->Particle2.Name)) ||
//          ((Path.Species(spec2).Type==PairArray(i)->Particle1.Name)&&
//           (Path.Species(spec1).Type==PairArray(i)->Particle2.Name))) {
//        if (PairMatrix(spec1,spec2) != NULL) {
//          perr << "More than one pair action for species types (" 
//               << PairArray(i)->Particle1.Name << ", "
//               << PairArray(i)->Particle2.Name << ")." << endl;
//          exit(-1);
//        }
//        verr << "Found PAfile for pair (" 
//             << Path.Species(spec1).Name << ", "
//             << Path.Species(spec2).Name << ")\n";
//        PairMatrix(spec1,spec2) = PairArray(i);
//        PairMatrix(spec2,spec1) = PairArray(i);
//        paUsed = true;
//      }
//    if (!paUsed) {
//      perr << "Warning:  Pair action for species types (" 
//         << PairArray(i)->Particle1.Name << ", "
//         << PairArray(i)->Particle1.Name << ") not used.\n";
//    }
//
//    PAIO.CloseFile();
//  }
//
//
//  if (readSpecificHeatFiles){
//    cerr<<"I READ SPECIFIC HEAT FILES"<<endl;
//    assert((in.ReadVar("TauValues",TauValues)));
//    assert(TauValues.size()==SpecificHeatPAFiles.size());
//    for (int i=0; i<SpecificHeatPAFiles.size() ; i++) {
//      // Allow for tilde-expansion in these files
//      string name = ExpandFileName(SpecificHeatPAFiles(i));
//      assert(PAIO.OpenFile (name));
//      SpecificHeatPairArray(i) = ReadPAFit (PAIO, TauValues(i), MaxLevels);   
//      PAIO.CloseFile();
//    }
//  }
//
//   
//  // Now check to make sure all PairActions that we need are defined.
//  for (int species1=0; species1<Path.NumSpecies(); species1++)
//    for (int species2=0; species2<Path.NumSpecies(); species2++)
//      if (PairMatrix(species1,species2) == NULL) {
//      if ((species1 != species2) || 
//          (Path.Species(species1).NumParticles > 1)) {
//        perr << "We're missing a PairAction for species1 = "
//             << Path.Species(species1).Name << " and species2 = "
//             << Path.Species(species2).Name << endl;
//        exit(1);
//      }
//      }
//  
//  in.ReadVar("UseLongRange", UseLongRange);
//  if (HaveLongRange()) {
//    perr << "*** Using long-range/short-range breakup. ***\n";
//    assert (in.ReadVar("UseBackground", LongRange.UseBackground));
//    LongRangePot.UseBackground = LongRange.UseBackground;
//    LongRangeRPA.UseBackground = LongRange.UseBackground;
//  }
//
////   if (longRange){
////     LongRange.Init(in);
////     if (UseRPA)
////       LongRangeRPA.Init(in);
////   }
//
////  Create nodal action objects
////   NodalActions.resize(PathData.Path.NumSpecies());
////   for (int spIndex=0; spIndex<PathData.Path.NumSpecies(); spIndex++) {
////     SpeciesClass &species = PathData.Path.Species(spIndex);
////     if (species.GetParticleType() == FERMION) {
////       if (species.NodeType == "FREE") 
////    NodalActions (spIndex) = new FreeNodalActionClass(PathData, spIndex);
////       else {
////    cerr << "Unrecognized node type " << species.NodeType << ".\n";
////    exit(EXIT_FAILURE);
////       }
////       NodalActions(spIndex)->Read(in);
////     }
////     else
////       NodalActions(spIndex) = NULL;
////   }
//  
//  ReadNodalActions (in);
//
//  if (UseNonlocal) {
//    if (FixedPhaseA == NULL) {
//      cerr << "Trying to use nonlocal action with a FixedPhase defined.\n";
//      abort();
//    }
//    Nonlocal.Setup(FixedPhaseA);
//  }
//
////   // Now create nodal actions for Fermions
////   NodalActions.resize(PathData.Path.NumSpecies());
////   for (int species=0; species<PathData.Path.NumSpecies(); species++) 
////     if (PathData.Path.Species(species).GetParticleType() == FERMION)
////       NodalActions(species) = new FPNodalActionClass(PathData, species);
////     else
////       NodalActions(species) = NULL;
//  
//  verr << "Finished reading the action.\n"; 
//
//  if (in.OpenSection("StructureReject")){
//    StructureReject.Read(in);
//    in.CloseSection();
//  }
//  if (in.OpenSection("Tether")){
//    Tether.Read(in);
//    in.CloseSection();
//  }
//  ///Reading in information for David long range action
//  if (PathData.Path.DavidLongRange){
//    DavidLongRange.Read(in);
//  }
//  if (PathData.Path.OpenPaths)
//    OpenLoopImportance.Read(in);
//}

