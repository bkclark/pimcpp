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
#include <omp.h>
#include "../PathDataClass.h"
#include "ShortRangeOn_diagonal_Class.h"
#include  "sys/time.h"
#include <omp.h>
#include <vector>
///DO NOT USE IF YOUR CUTOFF IS SUCH THAT ALL PARTICLE WILL BE INCLUDED IN ANY DIRECTION! THERE IS A BUG THAT WILL CAUSE IT TO BREAK!
///This has to be called after pathdata knows how many
///particles it has
void ShortRangeOn_diagonal_class::AcceptCopy (int slice1, int slice2)
{
  distancesNew.clear();
  distancesOld.clear();
  factorArrayNew.clear();
  factorArrayOld.clear();

}
void ShortRangeOn_diagonal_class::RejectCopy (int slice1, int slice2)
{
  distancesNew.clear();
  factorArrayNew.clear();

  distancesOld.clear();
  factorArrayOld.clear();


}

double 
ShortRangeOn_diagonal_class::residual_energy()
{
  int M=PathData.Path.NumTimeSlices()-1;
  int numParticles=PathData.Path.NumParticles();
  for (int i=0;i<PathData.Path.NumParticles();i++)
    ptcls(i)=i;
  gradVec=0.0;
  gradVecSquared=0.0;
  GradAction_help(0,M,ptcls,0,gradVec,gradVecSquared);
  double totalEnergy=0.0;
  for (int i=0;i<PathData.Path.NumParticles();i++){
    //    cerr<<gradVec(i)<<endl;
    double lambda=PathData.Path.ParticleSpecies(i).lambda;
    totalEnergy-=lambda/3.0*dot(gradVec(i),gradVec(i));
    totalEnergy+=1.0/3.0*gradVecSquared(i);
  }
  return totalEnergy;
}


//Gradient of the action only works in 3d!!!!!!!
void
ShortRangeOn_diagonal_class::GradAction_help(int slice1, int slice2, 
					     const Array<int,1> &ptcls, 
					     int level,
					     Array<dVec,1> &gradVec,
					     Array<double,1> &gradVecSquared)
{

  PathClass &Path = PathData.Path;
  int skip = (1<<level);
  assert (gradVec.size() == ptcls.size());
  assert(gradVecSquared.size()==ptcls.size());
  for (int pi=0; pi<ptcls.size(); pi++) {
    int ptcl1 = ptcls(pi);
    int species1 = Path.ParticleSpeciesNum(ptcl1);
    double lambda_i=PathData.Path.ParticleSpecies(ptcl1).lambda;
    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
      double lambda_j=PathData.Path.ParticleSpecies(ptcl2).lambda;
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      PairActionFitClass &PA=*(PathData.Actions.PairMatrix(species1,species2));
      if (ptcl1 != ptcl2) {
	for (int slice=slice1; slice<slice2; slice += skip) {
	  //	  double check=0.0;
	  dVec r, rp;
	  double rmag, rpmag, du_dq, du_dz;
	  Path.DistDisp(slice, slice+skip, ptcl1, ptcl2, rmag, rpmag, r, rp);
	  //	  check += dUdR(slice,ptcl1,ptcl2,level)+dUdR(slice+skip,ptcl1,ptcl2,level);
	  double q = 0.5*(rmag+rpmag);
	  double z = (rmag-rpmag);
	  double s2 = dot (r-rp, r-rp);
	  PA.Derivs(q,z,s2,level,du_dq, du_dz);
	  dVec rhat  = (1.0/rmag)*r;
	  dVec rphat = (1.0/rpmag)*rp;
	  
	  double g1 = 1.0;
	  double g2 = 1.0;

	  dVec grad_ij = -1.0* (g1*(0.5*du_dq + du_dz)*rhat + 
			  g2*(0.5*du_dq - du_dz)*rphat);
// 	  cerr<<"ptcls: "<<ptcl1<<" "<<ptcl2<<" "<<slice<<" "<<slice+skip<<endl;
// 	  cerr<<"rmags "<<rhat<<" "<<rphat<<" "<<rmag<<" "<<rpmag<<endl;
// 	  cerr<<"pre-vals "<<q<<" "<<z<<" "<<s2<<" "<<du_dq<<" "<<du_dz<<endl;
// 	  cerr<<"grad_ij is "<<grad_ij<<endl;
	  if (ptcl1!=ptcl2)
	    gradVec(pi) += grad_ij;
	  if (ptcl1<ptcl2)
	    gradVecSquared(pi) += dot(grad_ij,grad_ij)*(lambda_i+lambda_j);
	  
  
	  //	  cerr<<"VALS: "<<(g1*(0.5*du_dq + du_dz)*rhat + 
	  //			   g2*(0.5*du_dq - du_dz)*rphat)<<" "<<check<<endl;
	  // gradVec(pi) -= (0.5*du_dq*(rhat+rphat) + du_dz*(rhat-rphat));
	  /// Now, subtract off long-range part that shouldn't be in
	  /// here 
	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	    gradVec(pi) += 0.5*(PA.Ulong(level).Deriv(rmag)*g1*rhat+
				PA.Ulong(level).Deriv(rpmag)*g2*rphat);
	}
      }
    }
  }

  //  cerr<<"Gradient not working in 2d"<<endl;
  //  assert(1==2);

}


void ShortRangeOn_diagonal_class::Read(IOSectionClass& in)
{
  DoPtcl.resize(PathData.Path.NumParticles());
  todoIt.resize(PathData.Path.NumParticles());
  //  distancesNew.reserve(PathData.Path.NumParticles()*10*8);
  //  factorArrayNew.reserve(PathData.Path.NumParticles()*10*8);
  //  distancesOld.reserve(PathData.Path.NumParticles()*10*8);
  //  factorArrayOld.reserve(PathData.Path.NumParticles()*10*8);

  ptcls.resize(PathData.Path.NumParticles());
  gradVec.resize(PathData.Path.NumParticles());
  gradVecSquared.resize(PathData.Path.NumParticles());
}

ShortRangeOn_diagonal_class::ShortRangeOn_diagonal_class(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{
}


double 
ShortRangeOn_diagonal_class::SingleAction_fast (int slice1, int slice2,
				 const Array<int,1> &changedParticles,
				 int level)
{
#if NDIM==2
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  PathClass &Path = PathData.Path;
  // First, sum the pair actions
  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  
  int totalParticles=0;
  int startSlice=slice1;
  if (slice1==0 && slice2==PathData.Path.NumTimeSlices()-1)
    startSlice=slice1;
  else 
    startSlice=slice1+skip;
  for (int slice=startSlice;slice<slice2;slice+=skip){
    double factor=1.0;
    if (slice==slice1  || slice==slice2)
      factor=0.5;
    Path.DoPtcl=true;
    
    //     ///changed particles loop against themselvs
    for (int ptcl1Index=0;ptcl1Index<numChangedPtcls;ptcl1Index++){
      int ptcl1 = changedParticles(ptcl1Index);
      int species1=Path.ParticleSpeciesNum(ptcl1);
      Path.DoPtcl(ptcl1)=false;
      for (int ptcl2Index=ptcl1Index+1;ptcl2Index<numChangedPtcls;ptcl2Index++){
    	int ptcl2=changedParticles(ptcl2Index);
    	PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
    	dVec r, rp;
    	double rmag, rpmag;
	PathData.Path.DistDispFast(slice, ptcl1, ptcl2,
				   rmag, r);
    	double U;
	U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
    	TotalU += U;
      }
    }


    for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
      int ptcl1 = changedParticles(ptcl1Index);
      Path.DoPtcl(ptcl1) = false;
      int species1=Path.ParticleSpeciesNum(ptcl1);
      int xBox,yBox;
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox);
      //      cerr<<"Affected cells"<<Path.Cell.AffectedCells.size()<<endl;
      int numAffectedCells=Path.Cell.AffectedCells.size();

 {

      for (int cellVal=0;cellVal<numAffectedCells;cellVal++){ 
	int rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	int rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2)){ //I think this is ok
	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r;
	    double rmag;
	    PathData.Path.DistDispFast(slice, ptcl1,ptcl2,
				       rmag, r);
	    //	    double U=0.0;

	    double U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	    //  double U = factor*PA.Udiag(rmag, level);
	    TotalU += U;
	  }
	}
      }
    } //parallel end
    }
    

    
  } //end slice loop
  //  cerr<<"My total number of particles is "<<totalParticles<<endl;
  gettimeofday(&end,   &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);
  //  cerr<<"Time spent in diagonal class is "<<TimeSpent<<endl;
  double checkU=SingleAction_slow(slice1,slice2,changedParticles,level);
  cerr<<"CHECK: "<<checkU<<" "<<TotalU<<endl;
  return (TotalU);
#else
  cerr<<"not implemented short range o(n) diagonal class in 2d"<<endl;
#endif
}




double 
ShortRangeOn_diagonal_class::SingleAction (int slice1, int slice2,
				 const Array<int,1> &changedParticles,
				 int level)
{
#if NDIM==2
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  PathClass &Path = PathData.Path;
  // First, sum the pair actions
  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  if (GetMode()==NEWMODE){
    distances=&distancesNew;
    factorArray=&factorArrayNew;
  }
  else if (GetMode()==OLDMODE){
    distances=&distancesOld;
    factorArray=&factorArrayOld;

  }

  int totalParticles=0;
  int startSlice=slice1+skip;
  skip=skip*2;
  
  //  distances->clear();
  //  factorArray->clear();


  for (int slice=startSlice;slice<slice2;slice+=skip){
    //for (int slice=startSlice;slice<=slice2;slice+=skip){
    double factor=1.0;
    if (slice==slice1  || slice==slice2)
      factor=0.5;
    Path.DoPtcl=true;
    
    //     ///changed particles loop against themselvs
    for (int ptcl1Index=0;ptcl1Index<numChangedPtcls;ptcl1Index++){
      int ptcl1 = changedParticles(ptcl1Index);
      int species1=Path.ParticleSpeciesNum(ptcl1);
      Path.DoPtcl(ptcl1)=false;
      for (int ptcl2Index=ptcl1Index+1;ptcl2Index<numChangedPtcls;ptcl2Index++){
    	int ptcl2=changedParticles(ptcl2Index);
	//    	PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
    	dVec r, rp;
    	double rmag, rpmag;
	PathData.Path.DistDispFast(slice, ptcl1, ptcl2,
				   rmag, r);

	  distances->push_back(rmag);
	  factorArray->push_back(factor);

	//    	double U;
	//	U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	//    	TotalU += U;
      }
    }


    for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
      int ptcl1 = changedParticles(ptcl1Index);
      Path.DoPtcl(ptcl1) = false;
      int species1=Path.ParticleSpeciesNum(ptcl1);
      int xBox,yBox;
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox);
      //      cerr<<"Affected cells"<<Path.Cell.AffectedCells.size()<<endl;
      int numAffectedCells=Path.Cell.AffectedCells.size();

 {

      for (int cellVal=0;cellVal<numAffectedCells;cellVal++){ 
	int rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	int rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2)){ //I think this is ok
	    ///	    PairActionFitClass &PA = *(PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2)));
	    dVec r;
	    double rmag;
	    PathData.Path.DistDispFast(slice, ptcl1,ptcl2,
				       rmag, r);
	    //	    rmag=2.0;
	    //	    double U=0.0;

	      distances->push_back(rmag);
	      factorArray->push_back(factor);

	    //	    double U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	    //  double U = factor*PA.Udiag(rmag, level);
	    //	    TotalU += U;
	  }
	}
      }
    } //parallel end
  }
    

    
  } //end slice loop
  PairActionFitClass &PA = *(PairMatrix(0, 0));
  DavidPAClass *DPA=((DavidPAClass*)&PA);
  int loopSize=distances->size();
  double TotalUp=0.0;
  //   cerr<<"START"<<endl;
//   cerr<<"Diagonal end is "<<((DPA->UdiagSpline)(0)).End()<<endl;
//   cerr<<"Diagonal start is "<<((DPA->UdiagSpline)(0)).Start()<<endl;
//   cerr<<"dx is "<<((DPA->UdiagSpline)(4)).dx<<" "<<((DPA->UdiagSpline)(4)).dxInv<<endl;
//   cerr<<"second is "<<((DPA->UdiagSpline)(4)).Data(10)<<endl;
//   cerr<<"third  is "<<((DPA->UdiagSpline)(4)).Data.size()<<endl;
//   double deltar=((DPA->UdiagSpline(0)).End()-(DPA->UdiagSpline(0)).Start())/(1999);
//   cerr<<"deltar is "<<deltar<<endl;
//   for (double r=0.01;r<7;r+=deltar){
//     cerr<<r<<" ";
//     cerr<<DPA->Udiag(r, 0)<<" "; 
//     cerr<<DPA->UDiag_exact(r, 0)<<endl; 
    
//     double id = ((DPA->UdiagSpline)(4)).dxInv*(r-((DPA->UdiagSpline)(4)).GridStart);
//     double ipart;
//     double frac = modf (id, &ipart);

//     int i = (int) ipart;
//     cerr<<"ipart is "<<ipart<<" "<<frac<<" "<<i<<" "<<id<<endl;
//     cerr<<((1.0-frac)*((DPA->UdiagSpline)(4)).Data(i) + frac*((DPA->UdiagSpline)(4)).Data(i+1))<<endl;





//   }
//   exit(1);

  for (int i=0;i<loopSize;i++){
    //    cerr<<(*factorArray)[i]*DPA->Udiag((*distances)[i], level)<<" "; 
    //    cerr<<(*factorArray)[i]*DPA->UDiag_exact((*distances)[i], level)<<endl; 

      TotalU+=(*factorArray)[i]*DPA->UDiag_exact((*distances)[i], level); 
    //    TotalU+=DPA->Udiag((*distances)[i], level); 
    //    TotalUp+=(*factorArray)[i]*DPA->UDiag_exact((*distances)[i], level); 

//     double id = ((DPA->UdiagSpline)(level+4)).dxInv*((*distances)[i])-((DPA->UdiagSpline)(level+4)).GridStart;
//     double ipart;
//     double frac = modf (id, &ipart);

//     int i = (int) ipart;
//     //cerr<<"ipart is "<<ipart<<" "<<frac<<" "<<i<<" "<<id<<endl;
//     TotalU+=(*factorArray)[i]*(((1.0-frac)*((DPA->UdiagSpline)(level+4)).Data(i) + frac*((DPA->UdiagSpline)(level+4)).Data(i+1)));



  }

  //  cerr<<"Diff: "<<TotalU<<" "<<TotalUp<<endl;
  //  cerr<<setprecision(12)<<endl;
  //  cerr<<"Diff: "<<TotalU-TotalUp<<endl;

  //  cerr<<"My total number of particles is "<<totalParticles<<endl;
  gettimeofday(&end,   &tz);
  //  cerr<<"OUT IS "<< (double)(end.tv_sec-start.tv_sec) +
  //    1.0e-6*(double)(end.tv_usec-start.tv_usec)<<endl;

  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);
  //  cerr<<"Time spent in diagonal class is "<<TimeSpent<<endl;
  //  double checkU=SingleAction_slow(slice1,slice2,changedParticles,level);
  //  cerr<<"CHECK: "<<checkU<<" "<<TotalU<<endl;
  return (TotalU);
#else
  cerr<<"not implemented short range o(n) diagonal class in 2d"<<endl;
#endif
} 


double 
ShortRangeOn_diagonal_class::SingleAction_slow (int slice1, int slice2,
				 const Array<int,1> &changedParticles,
				 int level)
{
  struct timeval start, end; 
  struct timezone tz;
  gettimeofday(&start, &tz);

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


  for (int slice=slice1;slice<=slice2;slice+=skip){
    double factor=1.0;
    if (slice==slice1  || slice==slice2)
      factor=0.5;
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
    	PathData.Path.DistDispFast(slice, ptcl1, ptcl2,
				   rmag, r);
	//    	double s2 = dot (r-rp, r-rp);
	//   	double q = 0.5 * (rmag + rpmag); 
	//    	double z = (rmag - rpmag);
    	double U;
	//	gettimeofday(&start, &tz);
  
	U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	//	gettimeofday(&end,   &tz);
	//	TimeSpent += (double)(end.tv_sec-start.tv_sec) +
	//	  1.0e-6*(double)(end.tv_usec-start.tv_usec);
	
	//	U =factor*PA.Udiag(rmag,level);
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
	    PathData.Path.DistDispFast(slice, ptcl1,ptcl2,
				       rmag, r);
	    //	    double s2 = dot (r-rp, r-rp);
	    //	    double q = 0.5 * (rmag + rpmag);
	    //	    double z = (rmag - rpmag);
	    double U;
	    //	gettimeofday(&start, &tz);

	    U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	    //  gettimeofday(&end,   &tz);
	    //  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
	    //    1.0e-6*(double)(end.tv_usec-start.tv_usec);

	    //	    U = factor*PA.Udiag(rmag, level);
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
	    PathData.Path.DistDispFast(slice, ptcl1, ptcl2,
				       rmag, r);
	    //	    double s2 = dot (r-rp, r-rp);
	    //	    double q = 0.5 * (rmag + rpmag);
	    //	    double z = (rmag - rpmag);
	    double U;
	    //	    U = factor*PA.Udiag(rmag, level);
	    //	gettimeofday(&start, &tz);

	    U = factor*((DavidPAClass*)&PA)->UDiag_exact(rmag, level);
	    //  gettimeofday(&end,   &tz);
	    //  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
	    //    1.0e-6*(double)(end.tv_usec-start.tv_usec);

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
  gettimeofday(&end,   &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);
  //  cerr<<"Time spent in diagonal class is "<<TimeSpent<<endl;
  return (TotalU);
}




 double 
 ShortRangeOn_diagonal_class::d_dBeta (int slice1, int slice2,
				       int level)
{

#if NDIM==2
  struct timeval start, end;
  struct timezone tz;
  PathClass &Path = PathData.Path;


  double TotalU = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  if (GetMode()==NEWMODE){
    distances=&distancesNew;
    factorArray=&factorArrayNew;
  }
  else if (GetMode()==OLDMODE){
    distances=&distancesOld;
    factorArray=&factorArrayOld;

  }
  distances->clear();
  factorArray->clear();
  int totalParticles=0;
  int startSlice=slice1;
  for (int slice=startSlice;slice<=slice2;slice+=skip){
    double factor=1.0;
    if (slice==slice1  || slice==slice2)
      factor=0.5;
    else
      factor=1.0; 
    Path.DoPtcl=true;
    for (int ptcl1=0; ptcl1<PathData.Path.NumParticles(); ptcl1++){
      Path.DoPtcl(ptcl1) = false;
      int species1=Path.ParticleSpeciesNum(ptcl1);
      int xBox,yBox;
      Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox);
      int numAffectedCells=Path.Cell.AffectedCells.size();
      {
	
	for (int cellVal=0;cellVal<numAffectedCells;cellVal++){ 
	  int rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
	  int rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
	list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox).Particles(slice);
	for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
	  int ptcl2=*i;
	  if (Path.DoPtcl(ptcl2)){ //I think this is ok
	    dVec r;
	    double rmag;
	    PathData.Path.DistDispFast(slice, ptcl1,ptcl2,
				       rmag, r);
	    distances->push_back(rmag);
	    factorArray->push_back(factor);
	  }
	}
	}
      } //parallel end
    }
    
  } //end slice loop
  PairActionFitClass &PA = *(PairMatrix(0, 0));
  DavidPAClass *DPA=((DavidPAClass*)&PA);
  int loopSize=distances->size();
  double TotalUp=0.0;
  for (int i=0;i<loopSize;i++){
    TotalU+=(*factorArray)[i]*DPA->dUdiag((*distances)[i], 0); 
    //TotalU+=(*factorArray)[i]*DPA->dUdiag_fast((*distances)[i], 0); 
  }
  distances->clear();
  factorArray->clear();
  return (TotalU);
#else
  cerr<<"not implemented short range o(n) diagonal class in 2d"<<endl;
#endif
} 



    

double 
ShortRangeOn_diagonal_class::d_dBeta_slow(int slice1, int slice2,int level)
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
    



// double ShortRangeOn_diagonal_class::d_dBeta (int slice1, int slice2,
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
ShortRangeOn_diagonal_class::GetName()
{
  return "ShortRangeOn_diagonal";
}
