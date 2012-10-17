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

#include "PathDump.h"


void PathDumpClass::WriteBlock()
{
}

void PathDumpClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  in.ReadVar("AllClones", AllClones);
  DumpNodes = in.OpenSection("NodeDump");
  if (DumpNodes) {
    assert(in.ReadVar("Ptcl", NodePtcl));
    assert(in.ReadVar("Slice", NodeSlice));
    int nx, ny, nz;
    assert(in.ReadVar("Nx", nx));
    assert(in.ReadVar("Ny", ny));
    assert(in.ReadVar("Nz", nz));
    dVec box = PathData.Path.GetBox();
    Xgrid.Init(-0.5*box[0], 0.5*box[0], nx);
    Ygrid.Init(-0.5*box[1], 0.5*box[1], ny);
    Zgrid.Init(-0.5*box[2], 0.5*box[2], nz);
    if (PathData.Path.Communicator.MyProc()==0) {
      //      IOSection.WriteVar("NodeSlice", NodeSlice);
      //      IOSection.WriteVar("NodePtcl", NodePtcl);
      IOSection.NewSection("Xgrid");
      Xgrid.Write(IOSection);
      IOSection.CloseSection();
      IOSection.NewSection("Ygrid");
      Ygrid.Write(IOSection);
      IOSection.CloseSection();
      IOSection.NewSection("Zgrid");
      Zgrid.Write(IOSection);
      IOSection.CloseSection();
    }
    in.CloseSection();
  }
  else
    NodeSlice = 0;
  
  DumpRho = in.OpenSection("RhoDump");
  if (DumpRho) {
    int nx, ny, nz;
    assert (in.ReadVar("Nx", nx));
    assert (in.ReadVar("Ny", ny));
    assert (in.ReadVar("Nz", nz));
    Rho.resize(nx, ny, nz);
    in.CloseSection (); // "RhoDump"
  }
  
}

void PathDumpClass::Accumulate()
{
  std::cout << "Dumping Paths" << endl;
  if (!AllClones && (PathData.GetCloneNum() != 0))
    return;

  PathClass &Path= PathData.Path;
  int start, end, numProcs, myProc;
  numProcs = Path.Communicator.NumProcs();
  myProc   = Path.Communicator.MyProc();

  int refSave = PathData.Path.GetRefSlice();
  PathData.MoveRefSlice(0);

  if (DumpNodes)
    FindWorstBead(NodeSlice, NodePtcl);

  if (DumpRho) 
    if (PathData.Actions.NodalActions(0) != NULL)
      //      if (PathData.Actions.NodalActions(0)->Type() == GROUND_STATE_FP) {
      //	FixedPhaseActionClass &FP = 
      //	  *((FixedPhaseActionClass*)PathData.Actions.NodalActions(0));
      //	FP.CalcDensity(Rho);
      //	RhoVar.Write(Rho);
      //      }

  if (PathData.Path.OpenPaths){
    Array<double,1> tailLoc(NDIM);
    OpenLinkVar.Write((int)PathData.Path.OpenLink);
    for (int dim=0;dim<NDIM;dim++){
      tailLoc(dim)=PathData.Path(PathData.Path.OpenLink,
				 PathData.Path.NumParticles())[dim];
    }
    TailLocVar.Write(tailLoc);
    RefLinkVar.Write(PathData.Path.RefSlice);
    OpenLinkPtclVar.Write(PathData.Path.OpenPtcl);
  }

  Array<double,3> pathArray;
  int totalSlices = Path.TotalNumSlices;
  Path.SliceRange(numProcs-1, start,end);
  int maxShift = end-start;
  int slicesLeft = totalSlices;
  int offset = 0;
  int numPtcls = PathData.NumParticles();

  if (myProc == 0)
    pathArray.resize(numPtcls,totalSlices,NDIM);
  while (slicesLeft > maxShift) {
    int relNodeSlice = NodeSlice - offset;
    if (DumpNodes && (relNodeSlice>=0) && (relNodeSlice<Path.NumTimeSlices()))
      NodeDump ();
    // First copy
    PathData.MoveJoin(maxShift);
    if (myProc == 0)
      for (int i=0; i<maxShift; i++)
        for (int ptcl=0; ptcl < numPtcls; ptcl++)
          for (int dim=0; dim<NDIM; dim++)
            pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];

    // Now shift
    PathData.ShiftData(-maxShift);
    PathData.Join = 0;
    slicesLeft -= maxShift;
    offset += maxShift;
  }
  // Move join out of the way
  PathData.MoveJoin (Path.NumTimeSlices()-1);
  // Copy last part
  if (myProc == 0)
    for (int i=0; i<slicesLeft; i++)
      for (int ptcl=0; ptcl < numPtcls; ptcl++)
        for (int dim=0; dim<NDIM; dim++)
          pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];

  // Reset path to original position
  PathData.MoveRefSlice (refSave);

  Array<int,1> permVec(numPtcls);
  Path.TotalPermutation(permVec);

  // Now write
  PathVar.Write (pathArray);
  PermVar.Write (permVec);

  if (FirstTime && (myProc == 0)) {
    WriteInfo();
    Array<string,1> speciesNames(numPtcls);
    for (int spIndex=0;spIndex<PathData.NumSpecies();spIndex++){
      SpeciesClass &species = PathData.Path.Species(spIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    IOSection.WriteVar("Type","Path");
  }
  PathVar.Flush();
  FirstTime = false;
}

void
PathDumpClass::FindWorstBead(int &worstSlice, int &worstPtcl)
{
//   PathClass &Path = PathData.Path;
//   double logWorst = 0.0;
//   // First, find worst locally
//   for (int si=0; si<Path.NumSpecies(); si++) {
//     if ((PathData.Actions.NodalActions(si) != NULL) &&
// 	(Path.Species(si).lambda != 0.0))
//       if (PathData.Actions.NodalActions(si)->Type() == GROUND_STATE_FP) {
// 	FixedPhaseActionClass &FP = 
// 	  *((FixedPhaseActionClass*)PathData.Actions.NodalActions(si));
	
// 	int N = Path.Species(si).NumParticles;
// 	int M = Path.NumTimeSlices();
// 	int first = Path.Species(si).FirstPtcl;
// 	Array<double,2> Agrad2(M,N), Bgrad2(M,N);
// 	Path.SetIonConfig(0);
// 	Array<double,1> grad2;
// 	for (int slice=0; slice<M; slice++) {
// 	  grad2.reference(Agrad2(slice,Range::all()));
// 	  FP.CalcGrad2(slice, grad2);
// 	}
// 	Path.SetIonConfig(1);
// 	for (int slice=0; slice<M; slice++) {
// 	  grad2.reference(Bgrad2(slice,Range::all()));
// 	  FP.CalcGrad2(slice, grad2);
// 	}
// 	for (int slice=0; slice<M; slice++)
// 	  for (int ptcl=0; ptcl<N; ptcl++) {
// 	    double maxlog = fabs(Agrad2(slice,ptcl)-Bgrad2(slice,ptcl));
// 	    if (maxlog > logWorst) {
// 	      logWorst = maxlog;
// 	      worstSlice = slice;
// 	      worstPtcl  = ptcl + first;
// 	    }
// 	  }
//       }
//   }
//   Path.SetIonConfig(0);
//   // Now, compare among processors to find worst bead
//   Array<double,1> procWorst(1), allWorst(Path.Communicator.NumProcs());
//   Array<int,1>    procBeads(2), allBeads(Path.Communicator.NumProcs()*2);
//   procWorst(0) = logWorst;
//   int ts = Path.TotalNumSlices;
//   int myFirstSlice, myLastSlice;
//   Path.SliceRange(Path.Communicator.MyProc(), myFirstSlice, myLastSlice);
//   procBeads(0) = (worstSlice-Path.GetRefSlice()+myFirstSlice+ts)%ts;
//   procBeads(1) = worstPtcl;
//   Path.Communicator.AllGather(procWorst, allWorst);
//   Path.Communicator.AllGather(procBeads, allBeads);
//   for (int i=0; i<Path.Communicator.NumProcs(); i++) 
//     if (allWorst(i) > logWorst) {
//       logWorst = allWorst(i);
//       worstSlice = allBeads(2*i);
//       worstPtcl  = allBeads(2*i+1);
//     }

//   Path.Communicator.PrintSync();
//   cerr << "All worst is " << logWorst << " at (" << worstSlice << ", " 
//        << worstPtcl << ") myProc = " << Path.Communicator.MyProc() << endl;
//   cerr << "Ref slice = " << Path.GetRefSlice() << endl;
  
}

void
PathDumpClass::NodeDump()
{
  int myProc = PathData.Path.Communicator.MyProc();
  int speciesNum = PathData.Path.ParticleSpeciesNum(NodePtcl);
  if (PathData.Actions.NodalActions(speciesNum) != NULL) {
    NodeType nodeType = PathData.Actions.NodalActions(speciesNum)->Type();
    if ((nodeType==GROUND_STATE) && (myProc == 0))
      GroundStateNodeDump();
    else if (nodeType == GROUND_STATE_FP && (myProc==0)) {
      FixedPhaseNodeDump();
    }
    else if ((nodeType==FREE_PARTICLE) && (myProc==0))
      FreeParticleNodeDump();
  }
}

void
PathDumpClass::FixedPhaseNodeDump()
{
//   PathClass &Path= PathData.Path;

//   if (Path.UseCorrelatedSampling()) 
//     Path.SetIonConfig(0);

//   int relNodeSlice = NodeSlice + Path.GetRefSlice();
//   if (relNodeSlice > Path.TotalNumSlices)
//     relNodeSlice -= Path.TotalNumSlices;

//   /// HACK HACK HACK
//   /// Find ptcl closest to ion we moved.
// //   double closestDist = 1.0e100;
// //   int closestPtcl = 16;
// //   for (int ptcl=16; ptcl < 32; ptcl++) {
// //     dVec disp = Path(relNodeSlice,ptcl) - Path(relNodeSlice,0);
// //     Path.PutInBox(disp);
// //     double dist = sqrt(dot(disp,disp));
// //     if (dist < closestDist) {
// //       closestDist = dist;
// //       closestPtcl = ptcl;
// //     }
// //   }
// //   NodePtcl = closestPtcl;

//   int nx = Xgrid.NumPoints;
//   int ny = Ygrid.NumPoints;
//   int nz = Zgrid.NumPoints;

// //   Array<double,3> nodeDump(nx,ny,nz);
// //   int speciesNum = Path.ParticleSpeciesNum(NodePtcl);
// //   FixedPhaseActionClass &FP = 
// //   *((FixedPhaseActionClass*)PathData.Actions.NodalActions(speciesNum));

// //   Array<int,1> dummyPtcls;

// //   double grad2A = FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL);
// //   PathData.Path.SetIonConfig(1);
// //   double grad2B = FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL);
// //   PathData.Path.SetIonConfig(0);
// //   cerr << "Difference = " << grad2A - grad2B << endl;

// //   dVec savePos = Path(relNodeSlice,NodePtcl);

// //   // Try space warp on the whole time slice
// //   int N = Path.Species(speciesNum).NumParticles;
// //   int first = Path.Species(speciesNum).FirstPtcl;
// //   Array<dVec,1> saveSlice(N);
// //   dVec warpPos;
// //   for (int i=0; i<N; i++) {
// //     saveSlice(i) = Path(relNodeSlice,i+first);
// //     dVec tmp = saveSlice(i);
// //     if ((i+first) == NodePtcl)
// //       warpPos = tmp;
// //     // if B is happier than A
// //     if (grad2A > grad2B) {
// //       /// Warp from the happy B position to the A position.
// //       Path.WarpBtoA(tmp);
// //       Path.SetPos(relNodeSlice,i+first,tmp);
// //     }
// //     else {
// //       /// Warp from the happy A position to a hopefully better B position
// //       Path.WarpAtoB(tmp);
// //       Path.SetPos(relNodeSlice,i+first,tmp);
// //     }
// //   }
// //   if (grad2A > grad2B) {
// //     // Recalculate warped A slice
// //     Path.SetIonConfig(0);
// //     grad2A = FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL);
// //   }
// //   else {
// //     // Recalculate warped B slice.
// //     Path.SetIonConfig(1);
// //     grad2B = FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL);
// //     Path.SetIonConfig(0);
// //   }
// //   cerr << "Warped difference = " << (grad2A-grad2B) << endl;
// //   // Restore original position
// //   for (int i=0; i<N; i++)
// //     Path(relNodeSlice, i+first) = saveSlice(i);


// //   dVec newPos;
// //   for (int ix=0; ix<nx; ix++) {
// //     newPos[0] = Xgrid(ix);
// //     for (int iy=0; iy<ny; iy++) {
// //       newPos[1] = Ygrid(iy);
// //       for (int iz=0; iz<nz; iz++) {
// // 	newPos[2] = Zgrid(iz);
// // 	Path(relNodeSlice,NodePtcl) = newPos;
// // 	nodeDump(ix,iy,iz) = log10(FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL));
// //       }
// //     }
// //   }
// //   Path(relNodeSlice,NodePtcl) = savePos;
// //   NodeVarA.Write(nodeDump);
// //   if (PathData.Path.UseCorrelatedSampling()) {
// //     PathData.Path.SetIonConfig(1);
// //     for (int ix=0; ix<nx; ix++) {
// //       newPos[0] = Xgrid(ix);
// //       for (int iy=0; iy<ny; iy++) {
// // 	newPos[1] = Ygrid(iy);
// // 	for (int iz=0; iz<nz; iz++) {
// // 	  newPos[2] = Zgrid(iz);
// // 	  Path(relNodeSlice,NodePtcl) = newPos;
// // 	  nodeDump(ix,iy,iz) = log10(FP.CalcGrad2(relNodeSlice, dummyPtcls, UPDATE_ALL));
// // 	}
// //       }
// //     }
// //     Path(relNodeSlice,NodePtcl) = savePos;
// //     NodeVarB.Write(nodeDump);
// //     PathData.Path.SetIonConfig(0);
// //   }
// //   NodePtclVar.Write(NodePtcl);
// //   NodeSliceVar.Write(NodeSlice);
// //   Array<double,1> warpTmp(NDIM);
// //   for (int i=0; i<NDIM; i++)
// //     warpTmp(i) = warpPos[i];
// //   WarpPosVar.Write(warpTmp);
}


void
PathDumpClass::GroundStateNodeDump()
{
//   int nx = Xgrid.NumPoints;
//   int ny = Ygrid.NumPoints;
//   int nz = Zgrid.NumPoints;

//   PathClass &Path= PathData.Path;
//   Array<double,3> nodeDump(nx,ny,nz);
//   int speciesNum = Path.ParticleSpeciesNum(NodePtcl);
//   GroundStateNodalActionClass &GS = 
//     *((GroundStateNodalActionClass*)PathData.Actions.NodalActions(speciesNum));

//   int relNodeSlice = NodeSlice + Path.GetRefSlice();
//   if (relNodeSlice > PathData.Path.TotalNumSlices)
//     relNodeSlice -= Path.TotalNumSlices;

//   dVec savePos = Path(relNodeSlice,NodePtcl);
//   dVec newPos;
//   if (PathData.Path.UseCorrelatedSampling()) 
//     PathData.Path.SetIonConfig(0);
//   for (int ix=0; ix<nx; ix++) {
//     newPos[0] = Xgrid(ix);
//     for (int iy=0; iy<ny; iy++) {
//       newPos[1] = Ygrid(iy);
//       for (int iz=0; iz<nz; iz++) {
// 	newPos[2] = Zgrid(iz);
// 	Path(relNodeSlice,NodePtcl) = newPos;
// 	nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
//       }
//     }
//   }
//   Path(relNodeSlice,NodePtcl) = savePos;
//   NodeVarA.Write(nodeDump);
//   if (PathData.Path.UseCorrelatedSampling()) {
//     PathData.Path.SetIonConfig(1);
//     for (int ix=0; ix<nx; ix++) {
//       newPos[0] = Xgrid(ix);
//       for (int iy=0; iy<ny; iy++) {
// 	newPos[1] = Ygrid(iy);
// 	for (int iz=0; iz<nz; iz++) {
// 	  newPos[2] = Zgrid(iz);
// 	  Path(relNodeSlice,NodePtcl) = newPos;
// 	  nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
// 	}
//       }
//     }
//     Path(relNodeSlice,NodePtcl) = savePos;
//     NodeVarB.Write(nodeDump);
//     PathData.Path.SetIonConfig(0);
//   }
//   NodePtclVar.Write(NodePtcl);
}


void
PathDumpClass::FreeParticleNodeDump()
{
//   int nx = Xgrid.NumPoints;
//   int ny = Ygrid.NumPoints;
//   int nz = Zgrid.NumPoints;

//   PathClass &Path= PathData.Path;
//   Array<double,3> nodeDump(nx,ny,nz);
//   int speciesNum = PathData.Path.ParticleSpeciesNum(NodePtcl);
//   FreeNodalActionClass &GS = 
//     *((FreeNodalActionClass*)PathData.Actions.NodalActions(speciesNum));

//   int relNodeSlice = NodeSlice + Path.GetRefSlice();
//   if (relNodeSlice > Path.TotalNumSlices)
//     relNodeSlice -= Path.TotalNumSlices;

//   if (PathData.Path.UseCorrelatedSampling())
//     PathData.Path.SetIonConfig(0);
//   dVec savePos = Path(relNodeSlice,NodePtcl);
//   dVec newPos;
//   for (int ix=0; ix<nx; ix++) {
//     newPos[0] = Xgrid(ix);
//     for (int iy=0; iy<ny; iy++) {
//       newPos[1] = Ygrid(iy);
//       for (int iz=0; iz<nz; iz++) {
// 	newPos[2] = Zgrid(iz);
// 	Path(relNodeSlice,NodePtcl) = newPos;
// 	nodeDump(ix,iy,iz) = GS.Det(relNodeSlice);
//       }
//     }
//   }
//   Path(relNodeSlice,NodePtcl) = savePos;
//   NodeVarA.Write(nodeDump);
//   NodePtclVar.Write(NodePtcl);
}
