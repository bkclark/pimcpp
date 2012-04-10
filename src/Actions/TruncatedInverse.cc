
#define ORDER_N_FERMIONS
// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_blas.h>
#include "../PathDataClass.h"
#include <algorithm>
#include "TruncatedInverse.h"
#include <Common/MatrixOps/MatrixOps.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
// #ifdef ORDER_N_FERMIONS
extern "C"{
#include "../det_calc_uekt.h"
}
/// #endif
// int testme();


//We really want to know what the changed particle is here. Currently
//I'm going to try to just store it when we are doing the singleactionbelow
void 
TruncatedInverseClass::AcceptCopy(int slice1, int slice2)
{
//   cout<<"I'm being accepted! YEA!"<<ChangedColumn<<endl;
//   SetMode(OLDMODE);
//   cerr<<"Pre-checking"<<endl;
//   CheckDeterminantMatrix();
//   SetMode(NEWMODE);
//   cerr<<"post checking"<<endl;
//   for (int theRow=0;theRow<DetMatrix.extent(0);theRow++){
//     DetMatrix(ChangedColumn,theRow)=newCol(theRow);
//   }
//   for (int col=0;col<DetMatrix.extent(0);col++){
//     if (col!=ChangedColumn)
//       DetMatrix(col,ChangedColumn)=newCol(col);
//   }
  
  //  CheckDeterminantMatrix();
  //  cerr<<DetMatrix<<endl;
  //  cerr<<"My new determinant is "<<Determinant(DetMatrix)<<endl;
  //  BuildDeterminantMatrix();
}

struct intdouble
{
  double distance;
  int ptcl;
};

bool operator<(const doubleint& a, const doubleint& b) {
  return a.distance < b.distance;
}

bool operator<(const intdouble& a, const intdouble& b) {
  if (a.ptcl!=b.ptcl)
    return a.ptcl < b.ptcl;
  else 
    return a.distance<b.distance;
}


void 
TruncatedInverseClass::RejectCopy(int slice1, int slice2)
{
  //  cerr<<"I'm being rejected! YEA!"<<endl;
  
}

double 
TruncatedInverseClass::MinDistance(dVec oldDvec, dVec newDvec, int ptcl)
{
  dVec disp1=PathData.Path(ptcl,0)-oldDvec;
  double dist1=dot(disp1,disp1);
  dVec disp2=PathData.Path(ptcl,0)-newDvec;
  double dist2=dot(disp2,disp2);
  return min(dist1,dist2);
  

}

// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   Array<double,2> A12;
//   Array<double,2> A22;
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;

//   vector<doubleint> determinantPtcl;

//   int slice=0;
//   int ptcl1=ChangedColumn;
//   doubleint toPushBack;

//   dVec sortFrom=olddvec-newdvec;
//   PathData.Path.PutInBox(sortFrom);
//   sortFrom=sortFrom/2+newdvec;
  
//   for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
//     toPushBack.ptcl=ptcl;
//     toPushBack.distance=MinDistance(olddvec,newdvec,ptcl);
//     determinantPtcl.push_back(toPushBack);
//   }    
  
//   sort(determinantPtcl.begin(),determinantPtcl.end());
  
//   for (int n=1;n<100;n++){
//     SmallDetMatrix.resize(n,n);
//     for (int counter1=0;counter1<n;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<n;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
//     A12.resize(n,PathData.Path.NumParticles()-n);
//     for (int counter1=0;counter1<n;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=n;counter2<PathData.Path.NumParticles();counter2++){
// 	int ptcl2=determinantPtcl[counter2-n].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	A12(counter1,counter2-n)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
//     A22.resize(PathData.Path.NumParticles()-n,PathData.Path.NumParticles()-n);
//     for (int counter1=n;counter1<PathData.Path.NumParticles();counter1++){
//       int ptcl1=determinantPtcl[counter1-n].ptcl;
//       for (int counter2=n;counter2<PathData.Path.NumParticles();counter2++){
// 	int ptcl2=determinantPtcl[counter2-n].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	A12(counter1-n,counter2-n)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
    
//     DeterminantList(n-1)=Determinant(SmallDetMatrix);
//     Array<double,2> A11Inverse=Inverse(SmallDetMatrix);
//     //    OtherInfo(N-1,0)=
//     Array<double,2> tempMatrix=A11Inverse*A12;
//     Transpose(A12);
//     tempMatrix=A12*tempMatrix;
//     tempMatrix=A22-tempMatrix;
//     OtherInfo(n-1,0)=Determinant(tempMatrix);
    
//   }
  
//   //  cerr<<SmallDetMatrix<<endl;



//   //  cerr<<SmallDetMatrix<<endl;
//   cerr<<"out of it"<<endl;
// }


 
// void
// TruncatedInverseClass::GenerateNearbyParticleList(dVec fromLocation,list<doubleint> ptclList)
// {
//   Path.Cell.FindBox(fromLocation,xBox,yBox,zBox);  
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-sortFrom;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=0;
//       dVec a=dot(disp/sqrt(dot(disp,disp)),sortFrom/(sqrt(dot(sortFrom,sortFrom))))*sqrt(dot(disp,disp))*sortFrom/(sqrt(dot(sortFrom,sortFrom)));
//       dVec m=disp-a;
//       double zDist2=dot(a,a);
//       double xyDist2=dot(m,m);
//       double normDist=dot(disp,disp);
//       cerr<<zDist2<<" "<<xyDist2<<" "<<" "<<xyDist2+zDist2<<" "<<normDist<<endl;
//       assert(abs(normDist-(zDist2+xyDist2))<1e-3);
//       ////Use the below line if you are just doing normal
//       //distance. Use the line below that if you wish to have an
//       //ellipsoid instead of a circle.
//       toPushBack.distance=sqrt(dot(disp,disp));
//       //      toPushBack.distance=sqrt(1.1*zDist2+xyDist2);
//       ptclList.push_back(toPushBack);
//     }
//   }


// }  

// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   ModeType currMode=GetMode();
//   SetMode(NEWMODE);
//   ///Initializes some paramaters
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;
//   Array<bool,1> doneAlready(PathData.Path.NumParticles());

//   vector<doubleint> determinantPtcl;
//   vector<intdouble> tempdeterminantPtcl;
//   int xBox,yBox,zBox;
//   int slice=0;
//   int ptcl1=ChangedColumn;
//   intdouble toPushBack;
//   ///Initilialization done
//   doneAlready=false;
//   determinantPtcl.clear();
//   tempdeterminantPtcl.clear();

//   ///sortFrom should be the point in which you are sorting from.
//   ///If you wish to do min(dist(new,ptcl),dist(old,ptcl) that requires
//   ///a little more work.
//   dVec sortFrom=olddvec-newdvec;

//   //Loops over the particles andd puts them into a list with respect
//   //to the distance from the other particles. This only loops over
//   //particles in some radius from the changed particle which is what
//   //this is somewhat complicated.
//   PathData.Path.PutInBox(sortFrom);
//   sortFrom=sortFrom/2+newdvec;
//   Path.Cell.FindBox(sortFrom,xBox,yBox,zBox);  
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-sortFrom;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=0;
//       dVec a=dot(disp/sqrt(dot(disp,disp)),sortFrom/(sqrt(dot(sortFrom,sortFrom))))*sqrt(dot(disp,disp))*sortFrom/(sqrt(dot(sortFrom,sortFrom)));
//       dVec m=disp-a;
//       double zDist2=dot(a,a);
//       double xyDist2=dot(m,m);
//       double normDist=dot(disp,disp);
//       cerr<<zDist2<<" "<<xyDist2<<" "<<" "<<xyDist2+zDist2<<" "<<normDist<<endl;
//       assert(abs(normDist-(zDist2+xyDist2))<1e-3);
//       ////Use the below line if you are just doing normal
//       //distance. Use the line below that if you wish to have an
//       //ellipsoid instead of a circle.
//       toPushBack.distance=sqrt(dot(disp,disp));
//       //      toPushBack.distance=sqrt(1.1*zDist2+xyDist2);
//       tempdeterminantPtcl.push_back(toPushBack);
//     }
//   }
//   ////Looping over particels and addint them to the list is now done.

//   ///Remove duplicates from tempdeterminantPtcl list and put them in determinantPtcl
//   sort(tempdeterminantPtcl.begin(),tempdeterminantPtcl.end());
//   for (int counter=0;counter<tempdeterminantPtcl.size();counter++){
//     doubleint temp;
//     temp.ptcl=tempdeterminantPtcl[counter].ptcl;
//     temp.distance=tempdeterminantPtcl[counter].distance;
//     determinantPtcl.push_back(temp);
//     if (counter+1!=tempdeterminantPtcl.size()-1 && 
// 	counter!=tempdeterminantPtcl.size() &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+2].ptcl)
//       counter=counter+2;
//     else if (counter!=tempdeterminantPtcl.size()-1 && 
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl)
//       counter=counter+1;
//   }
//   ////Removing duplicates done


//   ///Sort the particles by their distance. Distance can be defined however
//   ///one desires.
//   sort(determinantPtcl.begin(),determinantPtcl.end());
//   for (int N=1;N<=determinantPtcl.size();N++){
//     SmallDetMatrix.resize(N,N);
//     for (int counter1=0;counter1<N;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<N;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
//     ///At this point, SmallDetMatrix has the matrix of size n
//     ///in it for the NEW or OLD configuration.  You can now accumulate things 
//     ///into DeterminantList or otherInfo
//     DeterminantList(N-1)=Determinant(SmallDetMatrix);
//   }
//   int N=20;
//   Array<double,1> A12mp1;
//   A12mp1.resize(N);
//   Array<double,2> smallInverse;
//   smallInverse.resize(N,N);
//   for (int i=1;i<N;i++){
//     int ptcl1=determinantPtcl[i].ptcl;
//     int ptcl2=N;
//     dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//     PathData.Path.PutInBox(disp);
//     double dist2=dot(disp,disp);
//     A12mp1(i)=exp(-T*dist2/(4.0*lambda));
//     smallInverse=Inverse(SmallDetMatrix);
//   }
//   Array<double,1> z;
//   z.resize(N);
//   MatVecProd(smallInverse,A12mp1,z);    
//   SetMode(OLDMODE);

//   ///And do it again for \tilde z

//   doneAlready=false;
//   determinantPtcl.clear();
//   tempdeterminantPtcl.clear();

//   ///sortFrom should be the point in which you are sorting from.
//   ///If you wish to do min(dist(new,ptcl),dist(old,ptcl) that requires
//   ///a little more work.
//   sortFrom=olddvec-newdvec;

//   //Loops over the particles andd puts them into a list with respect
//   //to the distance from the other particles. This only loops over
//   //particles in some radius from the changed particle which is what
//   //this is somewhat complicated.
//   PathData.Path.PutInBox(sortFrom);
//   sortFrom=sortFrom/2+newdvec;
//   Path.Cell.FindBox(sortFrom,xBox,yBox,zBox);  
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-sortFrom;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=0;
//       dVec a=dot(disp/sqrt(dot(disp,disp)),sortFrom/(sqrt(dot(sortFrom,sortFrom))))*sqrt(dot(disp,disp))*sortFrom/(sqrt(dot(sortFrom,sortFrom)));
//       dVec m=disp-a;
//       double zDist2=dot(a,a);
//       double xyDist2=dot(m,m);
//       double normDist=dot(disp,disp);
//       cerr<<zDist2<<" "<<xyDist2<<" "<<" "<<xyDist2+zDist2<<" "<<normDist<<endl;
//       assert(abs(normDist-(zDist2+xyDist2))<1e-3);
//       ////Use the below line if you are just doing normal
//       //distance. Use the line below that if you wish to have an
//       //ellipsoid instead of a circle.
//       toPushBack.distance=sqrt(dot(disp,disp));
//       //      toPushBack.distance=sqrt(1.1*zDist2+xyDist2);
//       determinantPtcl.push_back(toPushBack);
//     }
//   }
//   ///Sort the particles by their distance. Distance can be defined however
//   ///one desires.
//   sort(determinantPtcl.begin(),determinantPtcl.end());


//   ///At this point you can now have a list of particles that are near 
//   for (int N=1;N<=determinantPtcl.size();N++){
//     SmallDetMatrix.resize(N,N);
//     for (int counter1=0;counter1<N;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<N;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
//     ///At this point, SmallDetMatrix has the matrix of size n
//     ///in it for the NEW or OLD configuration.  You can now accumulate things 
//     ///into DeterminantList or otherInfo
//     DeterminantList(N-1)=Determinant(SmallDetMatrix);
//   }
//   N=20;
//   A12mp1;
//   A12mp1.resize(N);

//   for (int i=1;i<N;i++){
//     int ptcl1=determinantPtcl[i].ptcl;
//     int ptcl2=N;
//     dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//     PathData.Path.PutInBox(disp);
//     double dist2=dot(disp,disp);
//     A12mp1(i)=exp(-T*dist2/(4.0*lambda));
//     smallInverse=Inverse(SmallDetMatrix);
//   }
//   Array<double,1> zTilde;
//   zTilde.resize(N);
//   MatVecProd(smallInverse,A12mp1,zTilde);    

//   cerr<<"printing the difference"<<" ";
//   for (int i=0;i<z.size();i++)
//     cerr<<z(i)-zTilde(i)<<" ";
//   cerr<<endl;
//   SetMode(currMode);

// }

 


////////

double
TruncatedInverseClass::GetDistance(dVec sortFrom,int ptcl2)
{
  dVec disp=PathData.Path(0,ptcl2)-sortFrom;
  double distance=0;
  PathData.Path.PutInBox(disp);
  dVec a=dot(disp/sqrt(dot(disp,disp)),sortFrom/(sqrt(dot(sortFrom,sortFrom))))*sqrt(dot(disp,disp))*sortFrom/(sqrt(dot(sortFrom,sortFrom)));
  dVec m=disp-a;
  double zDist2=dot(a,a);
  double xyDist2=dot(m,m);
  double normDist=dot(disp,disp);
  //  cerr<<zDist2<<" "<<xyDist2<<" "<<" "<<xyDist2+zDist2<<" "<<normDist<<endl;
  assert(abs(normDist-(zDist2+xyDist2))<1e-3);
  ////Use the below line if you are just doing normal
  //distance. Use the line below that if you wish to have an
  //ellipsoid instead of a circle.
  distance=sqrt(dot(disp,disp));
  //      distance=sqrt(1.1*zDist2+xyDist2);
  return distance;
}

void 
TruncatedInverseClass::BuildOldMatrix(vector<doubleint> determinantPtcl,int size,
					   Array<double,2> &SmallDetMatrix)
{
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;

  SetMode(OLDMODE);
  int N=size;
  SmallDetMatrix.resize(N,N);
  for (int counter1=0;counter1<N;counter1++){
    int ptcl1=determinantPtcl[counter1].ptcl;
    for (int counter2=0;counter2<N;counter2++){
      int ptcl2=determinantPtcl[counter2].ptcl;
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
    }
  }
  
  

}


void 
TruncatedInverseClass::BuildNewMatrix(vector<doubleint> determinantPtcl,int size,
					   Array<double,2> &SmallDetMatrix)
{
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  SetMode(NEWMODE);
  int N=size;
  SmallDetMatrix.resize(N,N);
  for (int counter1=0;counter1<N;counter1++){
    int ptcl1=determinantPtcl[counter1].ptcl;
    for (int counter2=0;counter2<N;counter2++){
      int ptcl2=determinantPtcl[counter2].ptcl;
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
    }
  }
  
}


void 
TruncatedInverseClass::BuildSmallDeterminantMatrix()
{
  ModeType currMode=GetMode();
  ///Initializes some paramaters
  vector<doubleint> determinantPtcl;
  int xBox,yBox,zBox;
  int slice=0;
  int ptcl1=ChangedColumn;
  doubleint toPushBack;

  ///sortFrom should be the point in which you are sorting from.
  ///If you wish to do min(dist(new,ptcl),dist(old,ptcl) that requires
  ///a little more work.
  dVec sortFrom=olddvec-newdvec;
  PathData.Path.PutInBox(sortFrom);
  sortFrom=sortFrom/2+newdvec;

  //Loops over the particles andd puts them into a list with respect
  //to the distance from the other particles. 
  ///If the A lines are uncommented this only loops over
  //particles in some radius from the changed particle.
  //Otherwise it loops over all the particles!
  //A  Path.Cell.FindBox(sortFrom,xBox,yBox,zBox);  
  //A  for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
  //A    int rxbox,rybox,rzbox;
  //A    rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
  //A    rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
  //A    rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
  //A    list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
  //A    for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
  for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
    //A      int ptcl2=*i;
    toPushBack.ptcl=ptcl2;
    toPushBack.distance=GetDistance(sortFrom,ptcl2);
    determinantPtcl.push_back(toPushBack);
  }
//A  }

  ///Sort the particles by their distance. Distance can be defined however
  ///one desires.
  sort(determinantPtcl.begin(),determinantPtcl.end());
  //  cerr<<"ptcl 1 is "<<ptcl1<<endl;
  //  for (int i=0;i<determinantPtcl.size();i++){
  //    cerr<<determinantPtcl[i].ptcl<<" "<<determinantPtcl[i].distance<<endl;
  //  }
  BuildOldMatrix(determinantPtcl,determinantPtcl.size(),SmallDetMatrix);
  BuildNewMatrix(determinantPtcl,determinantPtcl.size(),SmallDetMatrixOld);
  Array<double,2> A11OldLarge = SmallDetMatrixOld(Range(0,SmallDetMatrixOld.extent(0)-1), Range(0,SmallDetMatrixOld.extent(0)-1));
  Array<double,2> A11NewLarge = SmallDetMatrix(Range(0,SmallDetMatrixOld.extent(0)-1), Range(0,SmallDetMatrixOld.extent(0)-1));
  
  ///EDIT ME! BELOW HERE!
  int A11size=30;
  Array<double,2> A11Old = SmallDetMatrixOld(Range(0,A11size), Range(0,A11size));
  Array<double,2> A11New = SmallDetMatrix(Range(0,A11size), Range(0,A11size));
  
  Array<double,2> A12Old = SmallDetMatrixOld(Range(A11size+1,determinantPtcl.size()),Range(0,A11size));
  Array<double,2> A12New=  SmallDetMatrix(Range(A11size+1,determinantPtcl.size()),Range(0,A11size));
  
  Array<double,2> A22Old = SmallDetMatrixOld(Range(A11size+1,determinantPtcl.size()),Range(A11size+1,determinantPtcl.size()));
  Array<double,2> A22New=  SmallDetMatrix(Range(A11size+1,determinantPtcl.size()),Range(A11size+1,determinantPtcl.size()));
  

  ///Check size check is indicative of bad configurations
  cerr<<"DET: "<<" "<<Determinant(SmallDetMatrix)/Determinant(SmallDetMatrixOld)<<" ";
  for (int sizeCheck=25;sizeCheck<30;sizeCheck++){
    Array<double,2> A11Newer=SmallDetMatrix(Range(0,sizeCheck), Range(0,sizeCheck));
    Array<double,2> A11Older=SmallDetMatrixOld(Range(0,sizeCheck), Range(0,sizeCheck));
    cerr<<Determinant(A11Newer)/Determinant(A11Older)<<" ";
  }
  
    
    
  cerr<<endl;



	 ///End of back configurations
  double oldDeterminant=Determinant(A11Old);
  double newDeterminant=Determinant(A11New);
//  cerr<<"The truncated method with "<<A11size<<" particles included is "<<newDeterminant/oldDeterminant<<endl;
//  cerr<<"The truncated method with all particles included is "<<Determinant(A11NewLarge)/Determinant(A11OldLarge)<<endl;
  cerr<<newDeterminant/oldDeterminant<<" ";
  DeterminantList(0)=newDeterminant/oldDeterminant;
////BBB  DeterminantList(0)=Determinant(A11NewLarge)/Determinant(A11OldLarge);
////BBB  cerr<<Determinant(A11NewLarge)/Determinant(A11OldLarge)<<" ";
  int oldZero=0;
  int newZero=0;
  int totalElements=PathData.Path.NumParticles()*PathData.Path.NumParticles();
  ///BBB  cerr<<"Actual det is "<<Determinant(A11NewLarge)<<" "
  ///BBB      <<Determinant(A11OldLarge)<<endl;
  Array<double,2> TestDetNew;
  TestDetNew.resize(A11NewLarge.shape());
  Array<double,2> TestDetOld;
  TestDetOld.resize(A11OldLarge.shape());

//  cerr<<A11NewLarge<<endl;
  double prob=10.0;
  double totAverage=0.0;
  double numTimes=0.0;
  int np=PathData.Path.NumParticles();
  for (int loop=0;loop<2;loop++){

    double factor=2.0*(double)np;
      newZero=0;
      oldZero=0;
      totalElements=PathData.Path.NumParticles()*PathData.Path.NumParticles();
      for (int counter=0;counter<A11OldLarge.extent(0);counter++)
	  for (int counter2=0;counter2<A11OldLarge.extent(1);counter2++){
	      TestDetNew(counter,counter2)=A11NewLarge(counter,counter2);
	      TestDetOld(counter,counter2)=A11OldLarge(counter,counter2);
	      //	      if (A11OldLarge(counter,counter2)<1.0/(factor*10.0) && 
	      //		  A11NewLarge(counter,counter2)<1.0/(factor*10.0)){
	      //	      if (A11OldLarge(counter,counter2)>1e-5 && 
	      //		  A11NewLarge(counter,counter2)>1e-5){
	      if (A11OldLarge(counter,counter2)<1.0/(factor*10.0) && 
		  A11NewLarge(counter,counter2)<1.0/(factor*10.0)){


		  double e1=A11OldLarge(counter,counter2);
		  double e2=A11NewLarge(counter,counter2);
		  
		  if (PathData.Path.Random.Local()>factor*e1)
		      TestDetOld(counter,counter2)=0.0;
		  else 
		    TestDetOld(counter,counter2)=1.0/factor;
		  if (PathData.Path.Random.Local()>factor*e2){
		    TestDetNew(counter,counter2)=0.0;
		    totalElements--;
		    oldZero++;
		  }
		  else 
		    TestDetNew(counter,counter2)=1.0/factor;
		  
	      }
	  }





//   for (int counter=0;counter<A11OldLarge.extent(0);counter++)
//       for (int counter2=0;counter2<A11OldLarge.extent(1);counter2++){
// 	  if (A11OldLarge(counter,counter2)<0.99){
// 	      if (PathData.Path.Random.Local()>1.0/0.1*A11OldLarge(counter,counter2)){
// 		  A11OldLarge(counter,counter2)=0.0;
// 		  oldZero++;
// 	      }
// 	      else {
// 		  A11OldLarge(counter,counter2)=0.1;
// 	      }
// 	  }
//       }
//   newZero=0;
//   for (int counter=0;counter<A11NewLarge.extent(0);counter++)
//       for (int counter2=0;counter2<A11NewLarge.extent(1);counter2++){
// 	  TestDet(counter,counter2)=A11NewLarge(counter,counter2);
// 	  if (A11NewLarge(counter,counter2)<0.99){
// 	      if (PathData.Path.Random.Local()>(1.0/0.1)*A11NewLarge(counter,counter2)){
// 		  TestDet(counter,counter2)=0.0;
// 		  newZero++;
// 	      }
// 	      else {
// 		  TestDet(counter,counter2)=0.1;
// 	      }
// 	  }
//       }
//  cerr<<Determinant(TestDet)<<" "<<PathData.Path.NumParticles()*PathData.Path.NumParticles()-newZero<<endl;
  totAverage=totAverage+Determinant(TestDetNew);      
  numTimes=numTimes+1.0;
  cerr<<"TOPRINTA: "<<Determinant(TestDetNew)<<" "<<Determinant(TestDetOld)<<" "
      <<Determinant(TestDetNew)/Determinant(TestDetOld)<<" "<<totalElements<<" "<<totAverage/numTimes<<endl;
  

  }

  // cerr<<500*500-oldZero<<" "<<500*500-newZero<<" "
  ////A  cerr<<Determinant(A11NewLarge)<<" "<<Determinant(A11OldLarge)<<" "
  ////A      <<Determinant(A11NewLarge)/Determinant(A11OldLarge)<<" "<<totalElements<<endl;
  

  
  SetMode(NEWMODE);
}


//////

  
//   ///At this point you can now have a list of particles that are near 
//   for (int N=1;N<=determinantPtcl.size();N++){
//     SmallDetMatrix.resize(N,N);
//     for (int counter1=0;counter1<N;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<N;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }
//     ///At this point, SmallDetMatrix has the matrix of size n
//     ///in it for the NEW or OLD configuration.  You can now accumulate things 
//     ///into DeterminantList or otherInfo
//     DeterminantList(N-1)=Determinant(SmallDetMatrix);
//   }
//   N=20;
//   A12mp1;
//   A12mp1.resize(N);

//   for (int i=1;i<N;i++){
//     int ptcl1=determinantPtcl[i].ptcl;
//     int ptcl2=N;
//     dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//     PathData.Path.PutInBox(disp);
//     double dist2=dot(disp,disp);
//     A12mp1(i)=exp(-T*dist2/(4.0*lambda));
//     smallInverse=Inverse(SmallDetMatrix);
//   }
//   Array<double,1> zTilde;
//   zTilde.resize(N);
//   MatVecProd(smallInverse,A12mp1,zTilde);    

//   cerr<<"printing the difference"<<" ";
//   for (int i=0;i<z.size();i++)
//     cerr<<z(i)-zTilde(i)<<" ";
//   cerr<<endl;
//   SetMode(currMode);
  

// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   cerr<<"Into it"<<endl;
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;
//   Array<bool,1> doneAlready(PathData.Path.NumParticles());
//   doneAlready=false;
//   vector<doubleint> determinantPtcl;
//   vector<intdouble> tempdeterminantPtcl;
//   int xBox,yBox,zBox;
//   int slice=0;
//   int ptcl1=ChangedColumn;
//   doneAlready=false;
//   intdouble toPushBack;
//   toPushBack.distance=0.0;
//   toPushBack.ptcl=ptcl1;
//   tempdeterminantPtcl.push_back(toPushBack);
//   doneAlready(ptcl1)=true;

//   dVec sortFrom=olddvec-newdvec;
//   PathData.Path.PutInBox(sortFrom);
//   sortFrom=sortFrom/2+newdvec;
//   Path.Cell.FindBox(sortFrom,xBox,yBox,zBox);  
//   //      cerr<<"Beginning"<<endl;
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-sortFrom;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=0;
//       dVec a=dot(disp/sqrt(dot(disp,disp)),sortFrom/(sqrt(dot(sortFrom,sortFrom))))*sqrt(dot(disp,disp))*sortFrom/(sqrt(dot(sortFrom,sortFrom)));
//       dVec m=disp-a;
//       double zDist2=dot(a,a);
//       double xyDist2=dot(m,m);
//       double normDist=dot(disp,disp);
//       cerr<<zDist2<<" "<<xyDist2<<" "<<" "<<xyDist2+zDist2<<" "<<normDist<<endl;
//       assert(abs(normDist-(zDist2+xyDist2))<1e-3);
//       //toPushBack.distance=sqrt(dot(disp,disp));
//       toPushBack.distance=sqrt(1.1*zDist2+xyDist2);
//       if (!doneAlready(ptcl2))
// 	tempdeterminantPtcl.push_back(toPushBack);
//       doneAlready(ptcl2)=true;
//     }
//   }



// //   Path.Cell.FindBox(newdvec,xBox,yBox,zBox);
// //   //      cerr<<"Beginning"<<endl;
// //   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
// //     int rxbox,rybox,rzbox;
// //     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
// //     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
// //     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
// //     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
// //     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
// //       int ptcl2=*i;
// //       //      if (!doneAlready(ptcl2))
// //       //	tempdeterminantPtcl.push_back(ptcl2);
// //       toPushBack.ptcl=ptcl2;
// //       dVec disp=PathData.Path(0,ptcl2)-newdvec;
// //       PathData.Path.PutInBox(disp);
// //       toPushBack.distance=sqrt(dot(disp,disp));
// //       tempdeterminantPtcl.push_back(toPushBack);
// //       doneAlready(ptcl2)=true;
// //     }
// //   }

//   sort(tempdeterminantPtcl.begin(),tempdeterminantPtcl.end());

//   for (int counter=0;counter<tempdeterminantPtcl.size();counter++){
//     doubleint temp;
//     temp.ptcl=tempdeterminantPtcl[counter].ptcl;
//     temp.distance=tempdeterminantPtcl[counter].distance;
//     determinantPtcl.push_back(temp);
//     if (counter+1!=tempdeterminantPtcl.size()-1 && 
// 	counter!=tempdeterminantPtcl.size() &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+2].ptcl)
//       counter=counter+2;
//     else if (counter!=tempdeterminantPtcl.size()-1 && 
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl)
//       counter=counter+1;
//   }
//   sort(determinantPtcl.begin(),determinantPtcl.end());
//   for (int counter=0;counter<determinantPtcl.size();counter++)
//     cerr<<"Distance: "<<counter<<" "<<determinantPtcl[counter].distance
// 	<<" "<<Path(0,determinantPtcl[counter].ptcl)<<endl;
//   //  for (int counter=0;counter<determinantPtcl.size();counter++)
//   //    cerr<<determinantPtcl[counter].ptcl<<" "<<determinantPtcl[counter].distance<<endl;


// //   int N=determinantPtcl.size();
// //   cerr<<"Size is "<<N<<endl;
// //   SmallDetMatrix.resize(N,N);
// //   int counter1=-1;

// //   vector<doubleint>::const_iterator iter1;
// //   vector<doubleint>::const_iterator iter2;
// //   for (iter1=determinantPtcl.begin();iter1!=determinantPtcl.end();iter1++){
// //     int ptcl1=(*iter1).ptcl;
// //     counter1++;
// //     int counter2=-1;
// //     for (iter2=determinantPtcl.begin();iter2!=determinantPtcl.end();iter2++){
// //       counter2++;
// //       int ptcl2=(*iter2).ptcl;
// //       dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// //       PathData.Path.PutInBox(disp);
// //       double dist2=dot(disp,disp);
// //       //      cerr<<ptcl1<<" "<<ptcl2<<endl;
// //       SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
// //     }
// //   }      



// //  int N=determinantPtcl.size();
// //  cerr<<"Size is "<<N<<endl;
  
//   for (int N=1;N<=determinantPtcl.size();N++){
//     SmallDetMatrix.resize(N,N);
//     for (int counter1=0;counter1<N;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<N;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }

//     DeterminantList(N-1)=Determinant(SmallDetMatrix);
//   }
  
//   //  cerr<<SmallDetMatrix<<endl;



//   //  cerr<<SmallDetMatrix<<endl;
//   cerr<<"out of it"<<endl;
// }



// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   cerr<<"Into it"<<endl;
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;
//   Array<bool,1> doneAlready(PathData.Path.NumParticles());
//   doneAlready=false;
//   vector<doubleint> determinantPtcl;
//   vector<intdouble> tempdeterminantPtcl;
//   int xBox,yBox,zBox;
//   int slice=0;
//   int ptcl1=ChangedColumn;
//   intdouble toPushBack;
//   toPushBack.distance=0.0;
//   toPushBack.ptcl=ptcl1;
//   tempdeterminantPtcl.push_back(toPushBack);
//   doneAlready(ptcl1)=true;
//   Path.Cell.FindBox(olddvec,xBox,yBox,zBox);
//   //      cerr<<"Beginning"<<endl;
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       //      if (!doneAlready(ptcl2))

//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-olddvec;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=sqrt(dot(disp,disp));
//       tempdeterminantPtcl.push_back(toPushBack);
//       doneAlready(ptcl2)=true;
//     }
//   }



//   Path.Cell.FindBox(newdvec,xBox,yBox,zBox);
//   //      cerr<<"Beginning"<<endl;
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       //      if (!doneAlready(ptcl2))
//       //	tempdeterminantPtcl.push_back(ptcl2);
//       toPushBack.ptcl=ptcl2;
//       dVec disp=PathData.Path(0,ptcl2)-newdvec;
//       PathData.Path.PutInBox(disp);
//       toPushBack.distance=sqrt(dot(disp,disp));
//       tempdeterminantPtcl.push_back(toPushBack);
//       doneAlready(ptcl2)=true;
//     }
//   }

//   sort(tempdeterminantPtcl.begin(),tempdeterminantPtcl.end());

//   for (int counter=0;counter<tempdeterminantPtcl.size();counter++){
//     doubleint temp;
//     temp.ptcl=tempdeterminantPtcl[counter].ptcl;
//     temp.distance=tempdeterminantPtcl[counter].distance;
//     determinantPtcl.push_back(temp);
//     if (counter+1!=tempdeterminantPtcl.size()-1 && 
// 	counter!=tempdeterminantPtcl.size() &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl &&
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+2].ptcl)
//       counter=counter+2;
//     else if (counter!=tempdeterminantPtcl.size()-1 && 
// 	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl)
//       counter=counter+1;
//   }
//   sort(determinantPtcl.begin(),determinantPtcl.end());
//   for (int counter=0;counter<determinantPtcl.size();counter++)
//     cerr<<"Distance: "<<counter<<" "<<determinantPtcl[counter].distance
// 	<<" "<<Path(0,determinantPtcl[counter].ptcl)<<endl;
//   //  for (int counter=0;counter<determinantPtcl.size();counter++)
//   //    cerr<<determinantPtcl[counter].ptcl<<" "<<determinantPtcl[counter].distance<<endl;


// //   int N=determinantPtcl.size();
// //   cerr<<"Size is "<<N<<endl;
// //   SmallDetMatrix.resize(N,N);
// //   int counter1=-1;

// //   vector<doubleint>::const_iterator iter1;
// //   vector<doubleint>::const_iterator iter2;
// //   for (iter1=determinantPtcl.begin();iter1!=determinantPtcl.end();iter1++){
// //     int ptcl1=(*iter1).ptcl;
// //     counter1++;
// //     int counter2=-1;
// //     for (iter2=determinantPtcl.begin();iter2!=determinantPtcl.end();iter2++){
// //       counter2++;
// //       int ptcl2=(*iter2).ptcl;
// //       dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// //       PathData.Path.PutInBox(disp);
// //       double dist2=dot(disp,disp);
// //       //      cerr<<ptcl1<<" "<<ptcl2<<endl;
// //       SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
// //     }
// //   }      



// //  int N=determinantPtcl.size();
// //  cerr<<"Size is "<<N<<endl;
  
//   for (int N=1;N<=determinantPtcl.size();N++){
//     SmallDetMatrix.resize(N,N);
//     for (int counter1=0;counter1<N;counter1++){
//       int ptcl1=determinantPtcl[counter1].ptcl;
//       for (int counter2=0;counter2<N;counter2++){
// 	int ptcl2=determinantPtcl[counter2].ptcl;
// 	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
// 	PathData.Path.PutInBox(disp);
// 	double dist2=dot(disp,disp);
// 	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//       }
//     }

//     DeterminantList(N-1)=Determinant(SmallDetMatrix);
//   }
  
//   //  cerr<<SmallDetMatrix<<endl;



//   //  cerr<<SmallDetMatrix<<endl;
//   cerr<<"out of it"<<endl;
// }


// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   cerr<<"Into it"<<endl;
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;

//   list<int> determinantPtcl;
//   int xBox,yBox,zBox;
//   int slice=0;
//   int ptcl1=ChangedColumn;
//   //  determinantPtcl.push_back(ptcl1);
//   Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox,zBox);
//   //      cerr<<"Beginning"<<endl;
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       determinantPtcl.push_back(ptcl2);
//     }
//   }
//   int N=determinantPtcl.size();
//   cerr<<"Size is "<<N<<endl;
//   SmallDetMatrix.resize(N,N);
//   int counter1=-1;

//   list<int>::const_iterator iter1;
//   list<int>::const_iterator iter2;
//   for (iter1=determinantPtcl.begin();iter1!=determinantPtcl.end();iter1++){
//     int ptcl1=*iter1;
//     counter1++;
//     int counter2=-1;
//     for (iter2=determinantPtcl.begin();iter2!=determinantPtcl.end();iter2++){
//       counter2++;
//       int ptcl2=*iter2;
//       //      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//       //      PathData.Path.PutInBox(disp);
//       //      double dist2=dot(disp,disp);
//       //      cerr<<ptcl1<<" "<<ptcl2<<endl;
//       SmallDetMatrix(counter1,counter2)=DetMatrix(ptcl1,ptcl2); // exp(-T*dist2/(4.0*lambda));
//     }
//   }      
//   //  cerr<<SmallDetMatrix<<endl;
//   cerr<<"out of it"<<endl;
// }

void 
TruncatedInverseClass::BuildDeterminantMatrix()
{
  cerr<<"Original buildling"<<endl;
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      DetMatrix(ptcl1,ptcl2)=exp(-T*dist2/(4.0*lambda));
    }
  }
  //  CheckDeterminantMatrix();
}


void 
TruncatedInverseClass::CheckDeterminantMatrix()
{

  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      if (DetMatrix(ptcl1,ptcl2)-exp(-T*dist2/(4.0*lambda))>=1e-3){
	cerr<<ptcl1<<" "<<ptcl2<<endl;
	cerr<<DetMatrix(ptcl1,ptcl2)<<endl;
	cerr<<exp(-T*dist2/(4.0*lambda))<<endl;
      }
      assert(DetMatrix(ptcl1,ptcl2)-exp(-T*dist2/(4.0*lambda))<=1e-3);
    }
  }
}




TruncatedInverseClass::TruncatedInverseClass (PathDataClass &pathData) :
					    
  NodalActionClass (pathData), 
  Path (pathData.Path)

{
}

void 
TruncatedInverseClass::Read (IOSectionClass &in)
{
  // Do nothing for now
  //  SetupFreeActions();
  cerr<<"Initializing"<<endl;
  NumTimes=0;
  int N=Path.NumParticles();
  DetMatrix.resize(N,N);
  SmallDetMatrix.resize(N,N);
  DeterminantList.resize(N);
  OtherInfo.resize(N,4);
  CutoffAverage=0.0;
}


// double 
// TruncatedInverseClass::SingleAction (int startSlice, int endSlice,
// 				  const Array<int,1> &changePtcls,
// 				  int level)
// {
// #ifdef ORDER_N_FERMIONS
//   ChangedColumn=changePtcls(0);
//   SetMode(OLDMODE);
//   olddvec=Path(0,ChangedColumn);
//   SetMode(NEWMODE);
//   newdvec=Path(0,ChangedColumn);
//   SetMode(OLDMODE);
//   //  BuildDeterminantMatrix();
//   BuildSmallDeterminantMatrix();
//   cerr<<"Matrix size is "<<SmallDetMatrix.extent(0)<<endl;
//   ofstream infile;
// //   infile.open("bigMatrixA");
// //   for (int i=0;i<DetMatrix.extent(0);i++){
// //     for (int j=0;j<DetMatrix.extent(1);j++){
// //       infile<<DetMatrix(i,j);
// //       infile<<" ";
// //     }
// //     infile<<endl;
// //   }
// //   infile.close();

// //   infile.open("smallMatrixA");
// //   for (int i=0;i<SmallDetMatrix.extent(0);i++){
// //     for (int j=0;j<SmallDetMatrix.extent(1);j++){
// //       infile<<SmallDetMatrix(i,j);
// //       infile<<" ";
// //     }
// //     infile<<endl;
// //   }
// //   infile.close();

//   //  cerr<<"Det matrix is "<<DetMatrix<<endl;
//   //  cerr<<"Small det matrix is "<<SmallDetMatrix<<endl;
//   //  //  //  //  Array<double,2> smallInverse=Inverse(SmallDetMatrix);
//   //  double det_old=Determinant(DetMatrix);
//   double det_old=1.0;
//   double small_det_old=Determinant(SmallDetMatrix);
//   int oldSize=SmallDetMatrix.extent(0);
//   SetMode(NEWMODE);
//   //  BuildDeterminantMatrix();
//   BuildSmallDeterminantMatrix();
//   cerr<<"Matrix size is "<<SmallDetMatrix.extent(0)<<endl;
//   //  //  //  //  Array<double,2> MultToDet;
//   //  //  //  //  MultToDet.resize(SmallDetMatrix.extent(0),SmallDetMatrix.extent(1));
//   //  //  //  //  MatMult(smallInverse,SmallDetMatrix,MultToDet);
//   //  MultToDet*=10;
//   //  cerr<<smallInverse<<endl;

// //   infile.open("bigMatrixA2");
// //   for (int i=0;i<DetMatrix.extent(0);i++){
// //     for (int j=0;j<DetMatrix.extent(1);j++){
// //       infile<<DetMatrix(i,j);
// //       infile<<" ";
// //     }
// //     infile<<endl;
// //   }
// //   infile.close();

// //   infile.open("smallMatrixA2");
// //   for (int i=0;i<SmallDetMatrix.extent(0);i++){
// //     for (int j=0;j<SmallDetMatrix.extent(1);j++){
// //       infile<<SmallDetMatrix(i,j);
// //       infile<<" ";
// //     }
// //     infile<<endl;
// //   }
// //   infile.close();

//   //  sleep(100);
//   //  cerr<<"Determinant: "<<Determinant(MultToDet)<<endl;
//   //  double det_new=Determinant(DetMatrix);
//   double det_new=1.0;
//   double small_det_new=Determinant(SmallDetMatrix);
//   int newSize=SmallDetMatrix.extent(0);
//   cerr<<det_old<<" "<<det_new<<" "<<small_det_old<<" "<<small_det_new<<endl;
//   cerr<<det_new/det_old<<" "<<small_det_new/small_det_old<<endl;
//   cerr<<oldSize<<" "<<newSize<<endl;
  
//   int cutOff=-1;
//   SetMode(OLDMODE);
//   double last_det_old=DeterminantList(newSize-1);
//   SetMode(NEWMODE);
//   double last_det_new=DeterminantList(newSize-1);
//   cerr<<"last det new"<<last_det_new<<" "<<last_det_old<<endl;
//   for (int i=0;i<newSize;i++){
//     SetMode(OLDMODE);
//     small_det_old=DeterminantList(i);
//     SetMode(NEWMODE);
//     small_det_new=DeterminantList(i);
//     cerr<<"Using "<<i<<" "<<small_det_new/small_det_old<<" "
// 	<<last_det_new/last_det_old<<endl;
//     if (abs(small_det_new/small_det_old-last_det_new/last_det_old)>1e-3)
//       cutOff=-1;
//     else if (cutOff==-1)
//       cutOff=i;
//   }
//   NumTimes++;
//   CutoffAverage+=cutOff;
//   cerr<<"Cutoff is "<<cutOff<<" "<<CutoffAverage/(double)NumTimes<<endl;
  
//   return 1.0;

// #else
//   cerr << "PIMC++ was not configured with --enable-on-fermions.\n"
//        << "This function doesn't work.  Please reconfigure.\n";
//   abort();
//   return 0.0;
// #endif

// }

double 
TruncatedInverseClass::SingleAction (int startSlice, int endSlice,
				  const Array<int,1> &changePtcls,
				  int level)
{
#ifdef ORDER_N_FERMIONS
  ModeType currMode=GetMode();
  //Initialize things
  ChangedColumn=changePtcls(0);
  SetMode(OLDMODE);
  olddvec=Path(0,ChangedColumn);
  SetMode(NEWMODE);
  newdvec=Path(0,ChangedColumn);
  //done initialization


  BuildSmallDeterminantMatrix();
  

  NumTimes++;
  SetMode(currMode);
  return DeterminantList(0);
//  return 1.0;


#else
  cerr << "PIMC++ was not configured with --enable-on-fermions.\n"
       << "This function doesn't work.  Please reconfigure.\n";
  abort();
  return 0.0;
#endif

}


//just here sot hatt it will compile
bool
TruncatedInverseClass::IsPositive(int x)
{
  return true;
}

double 
TruncatedInverseClass::d_dBeta (int slice1, int slice2, int level)
{ 
  return 0.0;
}

NodeType TruncatedInverseClass::Type()
{
  return FREE_PARTICLE;
}


bool
TruncatedInverseClass::IsGroundState()
{
  return (false);
}

void
TruncatedInverseClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "FREE_PARTICLE");
}

string
TruncatedInverseClass::GetName()
{
  return "TruncatedInverse";
}


