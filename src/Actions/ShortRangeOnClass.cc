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

#include "../PathDataClass.h"
#include "ShortRangeOnClass.h"
#include  "time.h"
///DO NOT USE IF YOUR CUTOFF IS SUCH THAT ALL PARTICLE WILL BE INCLUDED IN ANY DIRECTION! THERE IS A BUG THAT WILL CAUSE IT TO BREAK!
///This has to be called after pathdata knows how many
///particles it has
void ShortRangeOnClass::Read(IOSectionClass& in)
{
  DoPtcl.resize(PathData.Path.NumParticles());
}

ShortRangeOnClass::ShortRangeOnClass(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{
}

double 
ShortRangeOnClass::SingleAction (int slice1, int slice2,
				 const Array<int,1> &changedParticles,
				 int level)
{
  PathClass &Path = PathData.Path;
  int xEffect=Path.Cell.Xeffect;
  int yEffect=Path.Cell.Yeffect;
  int zEffect=Path.Cell.Zeffect;
  
  // First, sum the pair actions
  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);

  int totalParticles=0;
  Array<bool,1> todoIt(PathData.Path.NumParticles());
  todoIt=true;

  ///Check to see if the open is one of the changed particles  
  bool changingOpen=false;
  for (int i=0;i<changedParticles.size();i++)
    if (changedParticles(i)==PathData.Path.OpenPtcl)
      changingOpen=true;
  
  Array<int,1> changedPlusOpenPtcl;
  if (changingOpen || !PathData.Path.OpenPaths){
    changedPlusOpenPtcl.resize(changedParticles.size());
    for (int i=0;i<changedParticles.size();i++)
      changedPlusOpenPtcl(i)=changedParticles(i);
  }
  else{
    changedPlusOpenPtcl.resize(changedParticles.size()+1);
    for (int i=0;i<changedParticles.size();i++)
      changedPlusOpenPtcl(i)=changedParticles(i);
    changedPlusOpenPtcl(changedPlusOpenPtcl.size()-1)=PathData.Path.OpenPtcl;
  }


  for (int slice=slice1;slice<slice2;slice+=skip){
    
    for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
      Path.DoPtcl(ptcl)=true;
    //    Path.Cell.BinParticles(slice);
    //    Path.Cell.BinParticles(slice+skip);
    for (int ptcl1Index=0;ptcl1Index<changedPlusOpenPtcl.size();ptcl1Index++){
      int ptcl1 = changedPlusOpenPtcl(ptcl1Index);
      int species1=Path.ParticleSpeciesNum(ptcl1);
      Path.DoPtcl(ptcl1)=false;
      for (int ptcl2Index=ptcl1Index+1;ptcl2Index<changedPlusOpenPtcl.size();ptcl2Index++){
    	int ptcl2=changedPlusOpenPtcl(ptcl2Index);
    	PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
    	dVec r, rp;
    	double rmag, rpmag;
	if (slice==slice1)
	  totalParticles++;
    	PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
    			       rmag, rpmag, r, rp);
    	double s2 = dot (r-rp, r-rp);
   	double q = 0.5 * (rmag + rpmag); 
    	double z = (rmag - rpmag);
    	double U;
	//	U = ((DavidPAClass*)&PA)->UDiag_exact(q, level);

	U = PA.U(q,z,s2, level);
	// Subtract off long-range part from short-range action
	//	  if (PA.IsLongRange())
	//	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
    	TotalU += U;
      }
    }

    for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){

      todoIt=true;
      int ptcl1 = changedParticles(ptcl1Index);
      Path.DoPtcl(ptcl1) = false;
      int species1=Path.ParticleSpeciesNum(ptcl1);
      int xBox,yBox,zBox;
#if NDIM==3
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox,zBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);

	
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
#endif
#if NDIM==2
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
#endif

	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2)){ //I think this is ok
	    todoIt(ptcl2)=false;
	    //		Path.DoPtcl(ptcl2)=false;
	    
	    //		cerr<<"Particles: "<<ptcl1<<" "<<ptcl2<<endl;
	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r, rp;
	    double rmag, rpmag;
	    if (slice==slice1)
	      totalParticles++;
	    PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				   rmag, rpmag, r, rp);
	    double s2 = dot (r-rp, r-rp);
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    double U;
	    ///	U = ((DavidPAClass*)&PA)->UDiag_exact(q, level);
	    U = PA.U(q,z,s2, level);
	    // Subtract off long-range part from short-range action
	    //	  if (PA.IsLongRange())
	    //	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
	    //		//		if (TotalU>10000){
	    //		//		  cerr<<TotalU<<" "<<ptcl1<<" "<<ptcl2<<" "<<slice<<" "<<endl;
	    //		//		  cerr<<Path(slice,ptcl1)<<" "<<Path(slice,ptcl2)<<" "<<changedParticles<<endl;
	    //		//		}
	    TotalU += U;
	  }
	}
      } 
      if (PathData.Path.OpenPaths &&
	  slice+skip==PathData.Path.OpenLink && 
	  ptcl1==PathData.Path.OpenPtcl)
#if NDIM==3
	Path.Cell.FindBox(Path(slice+skip,PathData.Path.NumParticles()),
			  xBox,yBox,zBox);
      else
	Path.Cell.FindBox(Path(slice+skip,ptcl1),xBox,yBox,zBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
#endif	
#if NDIM==2
	Path.Cell.FindBox(Path(slice+skip,PathData.Path.NumParticles()),
			  xBox,yBox);
      else
	Path.Cell.FindBox(Path(slice+skip,ptcl1),xBox,yBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
#endif	

	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2) && todoIt(ptcl2)){ //I think this is ok
	    //		cerr<<"Particles: "<<ptcl1<<" "<<ptcl2<<endl;
	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r, rp;
	    double rmag, rpmag;
	    if (slice==slice1)
	      totalParticles++;
	    PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				       rmag, rpmag, r, rp);
	    double s2 = dot (r-rp, r-rp);
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    double U;
	    U = PA.U(q,z,s2, level);
	    //	    U = ((DavidPAClass*)&PA)->UDiag_exact(q, level);
	    // Subtract off long-range part from short-range action
	    //	  if (PA.IsLongRange())
	    //	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
	    TotalU += U;
	  }
	}
      }
    



      
    }
    

    
  } //end slice loop
  //  cerr<<"My total number of particles is "<<totalParticles<<endl;
  return (TotalU);
}



double 
ShortRangeOnClass::d_dBeta(int slice1, int slice2,int level)
{
  PathClass &Path = PathData.Path;
  int xEffect=Path.Cell.Xeffect;
  int yEffect=Path.Cell.Yeffect;
  int zEffect=Path.Cell.Zeffect;
  
  // First, sum the pair actions
  double dU = 0.0;
  
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  
  int totalParticles=0;
  Array<bool,1> todoIt(PathData.Path.NumParticles());
  todoIt=true;
  
  
  for (int slice=slice1;slice<slice2;slice+=skip){
    
    for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
      Path.DoPtcl(ptcl)=true;
    for (int ptcl1=0; ptcl1<PathData.Path.NumParticles(); ptcl1++){
      
      todoIt=true;
      Path.DoPtcl(ptcl1) = false;
      int species1=Path.ParticleSpeciesNum(ptcl1);
      int xBox,yBox,zBox;
#if NDIM==3
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox,zBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
	
	
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
#endif
#if NDIM==2
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
#endif

	
	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2)){ //I think this is ok
	    todoIt(ptcl2)=false;
	    //		Path.DoPtcl(ptcl2)=false;
	    
	    //		cerr<<"Particles: "<<ptcl1<<" "<<ptcl2<<endl;
	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r, rp;
	    double rmag, rpmag;
	    if (slice==slice1)
	      totalParticles++;
	    PathData.Path.DistDispFast(slice, slice+skip, ptcl1, ptcl2,
				   rmag, rpmag, r, rp);
	    double s2 = dot (r-rp, r-rp);
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    double U;
	    dU += PA.dU(q,z,s2, level);
	    
	  }
	}
      }
      if (PathData.Path.OpenPaths &&
	  slice+skip==PathData.Path.OpenLink && 
	  ptcl1==PathData.Path.OpenPtcl)
#if NDIM==3
	Path.Cell.FindBox(Path(slice+skip,PathData.Path.NumParticles()),
			  xBox,yBox,zBox);
      else
	Path.Cell.FindBox(Path(slice+skip,ptcl1),xBox,yBox,zBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
#endif
#if NDIM==2
	Path.Cell.FindBox(Path(slice+skip,PathData.Path.NumParticles()),
			  xBox,yBox);
      else
	Path.Cell.FindBox(Path(slice+skip,ptcl1),xBox,yBox);
      //      cerr<<"Beginning"<<endl;
      for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
	int rxbox,rybox,rzbox;
	rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	//	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
#endif
	
	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2) && todoIt(ptcl2)){ //I think this is ok
	    //		cerr<<"Particles: "<<ptcl1<<" "<<ptcl2<<endl;
	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r, rp;
	    double rmag, rpmag;
	    if (slice==slice1)
	      totalParticles++;
	    PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				   rmag, rpmag, r, rp);
	    double s2 = dot (r-rp, r-rp);
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    double U;
	    dU += PA.dU(q,z,s2, level);
	  }
	}
      }
      
      
      
      
      
    }
    

    
  } //end slice loop
  //  cerr<<"My total number of particles is "<<totalParticles<<endl;
  return (dU);
}



// double ShortRangeOnClass::d_dBeta (int slice1, int slice2,
// 				 int level)
// {
//   double levelTau=Path.tau;
//   int skip = 1<<level;
//   //  int slice2 = slice1 + (1<<level);
//   // Add constant part.  Note: we should really check the number of
//   // dimensions. 
//   double dU = 0.0;
//   for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
//     int species1=Path.ParticleSpeciesNum(ptcl1);
//     for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
//       for (int slice=slice1;slice<slice2;slice+=skip){
// 	dVec r, rp;
// 	double rmag, rpmag;
// 	PathData.Path.DistDisp(slice,slice+skip,ptcl1,ptcl2,rmag,rpmag,r,rp);
	
// 	double s2 = dot(r-rp, r-rp);
// 	double q = 0.5*(rmag+rpmag);
// 	double z = (rmag-rpmag);
	
// 	PairActionFitClass& pa=
// 	  *(PairMatrix(species1, PathData.Path.ParticleSpeciesNum(ptcl2)));
// 	dU += pa.dU(q, z, s2, level);
// 	// Subtract off long-range part from short-range action
// 	if (pa.IsLongRange())
// 	  dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
//       }
//     }
//   }
//   return dU;
// }


string 
ShortRangeOnClass::GetName()
{
  return "ShortRangeOn";
}
