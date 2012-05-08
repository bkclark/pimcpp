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



// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_blas.h>

#include "../PathDataClass.h"
#include "VariationalPI.h"
#include "../MatrixOps/MatrixOps.h"
#include <cstdlib>
#include <cstdio>
#include <string>

#ifdef ORDER_N_FERMIONS
extern "C"{
#include "../det_calc_uekt.h"
}
#endif
// int testme();


//We really want to know what the changed particle is here. Currently
//I'm going to try to just store it when we are doing the singleactionbelow
void 
VariationalPIClass::AcceptCopy(int slice1, int slice2)
{
   cout<<"I'm being accepted! YEA!"<<ChangedColumn<<endl;
//   SetMode(OLDMODE);
//   cerr<<"Pre-checking"<<endl;
//   CheckDeterminantMatrix();
//   SetMode(NEWMODE);
//   cerr<<"post checking"<<endl;
  for (int theRow=0;theRow<DetMatrix.extent(0);theRow++){
    DetMatrix(ChangedColumn,theRow)=newCol(theRow);
  }
  for (int col=0;col<DetMatrix.extent(0);col++){
    if (col!=ChangedColumn)
      DetMatrix(col,ChangedColumn)=newCol(col);
  }
  
  CheckDeterminantMatrix();
  //  cerr<<DetMatrix<<endl;
  //  cerr<<"My new determinant is "<<Determinant(DetMatrix)<<endl;
  //  BuildDeterminantMatrix();
}

void 
VariationalPIClass::RejectCopy(int slice1, int slice2)
{
    cerr<<"I'm being rejected! YEA!"<<endl;
  
}


void 
VariationalPIClass::BuildDeterminantMatrix()
{
  //HACK!
//  cerr<<"Initializing"<<endl;
  int N=Path.NumParticles();
//  cerr<<"the number of particles is"<<N<<endl;
  DetMatrix.resize(N,N);
  u.resize(N);
  newCol.resize(N);
  //HACK
//  cerr<<"Original buildling"<<endl;
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
VariationalPIClass::CheckDeterminantMatrix()
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




VariationalPIClass::VariationalPIClass (PathDataClass &pathData) :
					    
  NodalActionClass (pathData), 
  Path (pathData.Path)

{
}

void 
VariationalPIClass::Read (IOSectionClass &in)
{
  // Do nothing for now
  //  SetupFreeActions();
  cerr<<"Initializing Iterative Eigenvalue solver"<<endl;
  int N=Path.NumParticles();
  cerr<<"the number of particles is"<<N<<endl;
  DetMatrix.resize(N,N);
  u.resize(N);
  newCol.resize(N);
  //  testme();
  BuildDeterminantMatrix();

}


void 
VariationalPIClass::calc_u(const Array<int,1> &changePtcls)
{
  double lambda=PathData.Path.Species(0).lambda;
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  assert(changePtcls.size()==1);
  assert(changePtcls(0)<PathData.Path.NumParticles());
  int ptcl1=changePtcls(0);
  //  cerr<<"ptcl1 is "<<ptcl1<<" "<<u.size()<<endl;
  for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
    dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
    PathData.Path.PutInBox(disp);
    double dist2=dot(disp,disp);
    newCol(ptcl2)=exp(-T*dist2/(4.0*lambda));
    u(ptcl2)=newCol(ptcl2)-DetMatrix(ptcl1,ptcl2);
    
  }
}

int matvec(double *x,double *b,int N,void *A){
  
  //  cerr<<"Inside matvec"<<endl;
  //  cerr<<A<<endl;
  //  cerr<<((double**)A)[0]<<endl;
  //  cerr<<b<<endl;
  //  cerr<<x<<endl;
  for (int counter=0;counter<N;counter++){
    b[counter]=0.0;
    for (int counter2=0;counter2<N;counter2++){
      b[counter]=b[counter]+((double*)A)[counter2*N+counter]*x[counter2];
    }
  }
  //  cerr<<"out of matvec"<<endl;
  return 0;
}

// int matvec_blitz(double *x, double *b, void *A){
  
//   int N=((Array<double,2>*)A).extent(0);
//   Array<double,1> xBlitz(x, shape(N), neverDeleteData);
//   Array<double,1> AxBlitz(b, shape(N), neverDeleteData);
//   MatVecProd ((Array<double,2>*)A, xBlitz, AxBlitz);
//   return 0;
// }

// void set_gsl_vector_from_double(gsl_vector *V, double* x, int N){
//   V->size=N;V->stride=1;V->data=x;V->block=NULL;V->owner=0;  
// }


// int matvec_gsl_matrix(double *x,double *b,int N,void *params){
//   gsl_vector Xgsl,Bgsl;
//   set_gsl_vector_from_double(&Xgsl,x,N);
//   set_gsl_vector_from_double(&Bgsl,b,N);
//   //typedef struct
//   //{
//   //  size_t size;
//   //  size_t stride;
//   //  double *data;
//   //  gsl_block *block;
//   //  int owner;
//   //}
//   // gsl_vector;

//   //int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)  
//   int toReturn= gsl_blas_dgemv(CblasNoTrans,1,(gsl_matrix*)params,&Xgsl,0,&Bgsl);  
//   //  cerr<<"I am returning "<<toReturn;
//   return toReturn;
// }/*end matvec_gsl_matrix*/


double 
VariationalPIClass::SingleAction (int startSlice, int endSlice,
				  const Array<int,1> &changePtcls,
				  int level)
{
#ifdef ORDER_N_FERMIONS
  //  cerr<<"Calling single action"<<endl;
  //  ModeType currMode=PathData.Path.Path.GetMode();
//    SetMode(NEWMODE); //a
//    BuildDeterminantMatrix();//a
//    double new_det=Determinant(DetMatrix);//a

   SetMode(OLDMODE);
//   BuildDeterminantMatrix();//a
//   double old_det=Determinant(DetMatrix);//a
   SetMode(NEWMODE);
//   cerr<<old_det<<" "<<new_det<<" "<<old_det/new_det<<" ";//a

  //  cerr<<"Built determinant matrix"<<endl;
  //  cerr<<"My determinanat is "<<Determinant(DetMatrix)<<endl;
  //  DetOrderN DetOrderNInstance;
  drc_uekt_vanilla_parms parms={1e-3,2000};
  
//   SetMode(OLDMODE);
//   cerr<<"Checking"<<endl;
//   CheckDeterminantMatrix();
//   SetMode(NEWMODE);


  calc_u(changePtcls);
  //  cerr<<"calculated change "<<changePtcls(0)<<endl;
  double det_ratio;
  // testme();
//   cerr<<"The address you shoudl look at is "<<DetMatrix.data()<<" "
//       <<u.data()<<endl;
  
//    gsl_matrix *myMatrix;
//    myMatrix=gsl_matrix_calloc(DetMatrix.extent(0),DetMatrix.extent(1));
//    for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
//     for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
//        gsl_matrix_set(myMatrix,ptcl1,ptcl2,DetMatrix(ptcl1,ptcl2));
//      }
//    }
    
  ChangedColumn=changePtcls(0);
#ifdef ORDER_N_FERMIONS
  ////BUG!!  det_ratio_calculator_uekt_symmetric_vanilla_value(//matvec_gsl_matrix,
  ///BUG!						    matvec,
  ///BUG!						    DetMatrix.data(), 
  ///BUG!						    //(void*)myMatrix,
  ///BUG!						    u.size(), 
  ///BUG!						    u.data(), changePtcls(0),
  //BUG!    						    0, &det_ratio, 
  ///BUG!						    (void*)&parms);
#endif
 //  BuildDeterminantMatrix();
  //  cerr<<"My old determinanat is "<<Determinant(DetMatrix)<<endl;
  //  cerr<<"Called ratio calculator"<<endl;
  ////////  return new_det/old_det;
  //Somewhat of a hack
  if (isnan(det_ratio)){
    cerr<<"My det ratio is "<<0.0<<endl;
    return 0.0;
  }
  //  SetMode(currMode);
  
//  cerr<<"The iterative solver answer is "<<1.0/abs(det_ratio)<<endl;
  cerr<<1.0/abs(det_ratio)<<endl;
  
  return 1.0/abs(det_ratio);
  //  return 0.0;
#else
  cerr << "PIMC++ was not configured with --enable-on-fermions.\n"
       << "This function doesn't work.  Please reconfigure.\n";
  abort();
  return 0.0;
#endif

}


//just here sot hatt it will compile
bool
VariationalPIClass::IsPositive(int x)
{
  return true;
}

double 
VariationalPIClass::d_dBeta (int slice1, int slice2, int level)
{ 
  return 0.0;
}

NodeType VariationalPIClass::Type()
{
  return FREE_PARTICLE;
}


bool
VariationalPIClass::IsGroundState()
{
  return (false);
}

void
VariationalPIClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "FREE_PARTICLE");
}


string
VariationalPIClass::GetName()
{
  return "VariationalPI";
}

