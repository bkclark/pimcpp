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

#include "../PathDataClass.h"
#include "ShortRangeApproximateClass.h"

///This has to be called after pathdata knows how many
///particles it has
void ShortRangeApproximateClass::Read(IOSectionClass& in)
{
  DoPtcl.resize(PathData.Path.NumParticles());
}

ShortRangeApproximateClass::ShortRangeApproximateClass(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{
}

double 
ShortRangeApproximateClass::SingleAction (int slice1, int slice2,
					  const Array<int,1> &changedParticles,
					  int level)
{
  PathClass &Path = PathData.Path;
  //HACK! HACK! HACK!
  //  return 0.0;
  // First, sum the pair actions
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;

  double TotalU = 0.0;
  Array<double,1> TotalUArray(PathData.Path.NumTimeSlices());
  TotalUArray=0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);

    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
      if (Path.DoPtcl(ptcl2)){
	PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	for (int slice=slice1;slice<slice2;slice+=skip){
	  dVec r, rp;
	  double rmag, rpmag;

	  PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				 rmag, rpmag, r, rp);
	  double s2 = dot (r-rp, r-rp);
	  double q = 0.5 * (rmag + rpmag);
	  double z = (rmag - rpmag);

	  double U;
	  //	  if (ptcl1<3 && ptcl2<3)
	  //	    U=0.0; //*PA.U(q,z,s2,level);
	  //	  else
	    U=(PA.U(rmag,0,0,level)+PA.U(rpmag,0,0,level))/2.0;
	  //	  if (ptcl2==4 && ptcl1==0 && slice==8)
	  //	    cerr<<"MY approximate U is "<<U<<endl;

	  //	  U = PA.U(q,z,s2, level);
	  // Subtract off long-range part from short-range action
	  if (PA.IsLongRange())
	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
	  TotalU += U;
	  //	  TotalUArray(slice)+=U;
	  TotalUArray(slice)+=0.5*PA.U(rmag,0,0,level);
	  TotalUArray(slice+skip)+=0.5*PA.U(rpmag,0,0,level);
	}
      }
    }
  }

//   if (PathData.Path.LongRange){
//     // Now add in the long-range part of the action
//     // In primitive form, end slices get weighted by 1/2.  Others by 1.
//     for (int slice=slice1;slice<=slice2;slice+=skip) 
//       if ((slice==slice1) || (slice == slice2))
// 	TotalU += 0.5*LongRange_U (slice, level);
//       else
// 	TotalU += LongRange_U (slice, level);
//   }
//   return (TotalU);

//   if (PathData.Path.LongRange){
//     // Now add in the long-range part of the action
//     // In primitive form, end slices get weighted by 1/2.  Others by 1.
//     if (UseRPA) {
//       for (int slice=slice1;slice<=slice2;slice+=skip) 
// 	if ((slice==slice1) || (slice == slice2))
// 	  TotalU += 0.5*LongRange_U_RPA (slice, level);
// 	else
// 	  TotalU += LongRange_U_RPA (slice, level);
//     }
//     else {
//       for (int slice=slice1;slice<=slice2;slice+=skip) 
// 	if ((slice==slice1) || (slice == slice2))
// 	  TotalU += 0.5*LongRange_U (slice, level);
// 	else
// 	  TotalU += LongRange_U (slice, level);
//     }
//   }
//  for (int counter=0;counter<TotalUArray.size();counter++){
    //    cerr<<"My approximate link "<<counter+1<<" "<<counter+2<<" is "
  //    cerr <<TotalUArray(counter)<<endl;
  //  }

  return (TotalU);
}



double ShortRangeApproximateClass::d_dBeta (int slice1, int slice2,
				 int level)
{
  double levelTau=Path.tau;
  //  int slice2 = slice1 + (1<<level);
  for (int i=0; i<level; i++) 
    levelTau *= 2.0;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double dU = 0.0;
  const int NumImage=1;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      dVec r, rp;
      double rmag, rpmag;
      PathData.Path.DistDisp(slice1, slice2, ptcl1, ptcl2, rmag,rpmag,r,rp); 
      
      double s2 = dot(r-rp, r-rp);
      double q = 0.5*(rmag+rpmag);
      double z = (rmag-rpmag);

      PairActionFitClass& pa=*(PairMatrix(species1, PathData.Path.ParticleSpeciesNum(ptcl2)));
      dU += pa.dU(q, z, s2, level);
      // Subtract off long-range part from short-range action
      if (pa.IsLongRange())
	dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
    }
  }
  
  return dU;
}


string 
ShortRangeApproximateClass::GetName()
{
  return "ShortRangeApproximate";
}
