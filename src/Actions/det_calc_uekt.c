#ifdef ORDER_N_FERMIONS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "det_calc_uekt.h"
#ifdef ORDER_N_FERMIONS
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

/* External Functions */
double ddot_(const int* N, const double *X, const int *incX, const double *Y, const int *incY);
int cg_(int *n, double* b, double* x, double* work, int *ldw, int *iter, double* resid, 
	int (*matvec)(), int (*psolve)(), int *info);
int gmres_(const int * N, const double *b, double * x, const int*restrt, double *work, const int *ldw, double *h, const int *ldh, int *iter, double  *resid,int (*matvec)(), int (*psolve)(), int *info);


struct drc_uekt_common_block{
  int N;
  int (*matvec)(double*,double*,int,void*);
  void *matrix_parameters; 
};


struct drc_uekt_common_block CDATA;

int drc_uekt_identity_preconditioner(double* x, double*y){
  memcpy(x,y,CDATA.N*sizeof(double));
  return 1;
}/*end identity preconditioner*/

int drc_uekt_templates_style_fcall(double *alpha, double* x, double *beta, double*y){
  int i,N=CDATA.N;
  double* tmp=(double*)calloc(N,sizeof(double));

  /* Compute y := alpha*A*x + beta*y */
  (*(CDATA.matvec))(x,tmp,N,CDATA.matrix_parameters);
  for(i=0;i<N;i++) y[i]=(*alpha)*tmp[i] + (*beta) * y[i];

  /* Cleanup */
  free(tmp);
  return 1;
}/*end int*/


/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/

int testme()
{
/* Allocate Space */
  
  return 2;
}
  

void det_ratio_calculator_uekt_symmetric_vanilla_value(int (*func)(double*,double*,int,void*),
				       void * data, int N, double * u, int col_id, 
				      double die_roll, double* det_ratio,  void *params){

  static int firstTime=0;
  static double *mv1;
  static double *mv2;
  static double *ekt;
  static double *work;
  if (firstTime==0){
    mv1=(double*)calloc(N,sizeof(double));
    mv2=(double*)calloc(N,sizeof(double));
    ekt=(double*)calloc(N,sizeof(double));
    //  double *work=(double*)malloc(4*N*sizeof(double));
    work=(double*)calloc(4*N,sizeof(double));
    firstTime=1;
  }
/*   /\* Allocate Space *\/ */
/*   double *mv1=(double*)calloc(N,sizeof(double)); */
/*   double *mv2=(double*)calloc(N,sizeof(double)); */
/*   double *ekt=(double*)calloc(N,sizeof(double)); */
/*   //  double *work=(double*)malloc(4*N*sizeof(double)); */
/*   double *work=(double*)calloc(4*N,sizeof(double)); */
  double lambda[2];
  double tol;
  int iters[2];
  int info,rv,one=1;

  /* Parse Parameters */
  struct drc_uekt_vanilla_parms* OPTS = (struct drc_uekt_vanilla_parms*) params;
  CDATA.N=N; CDATA.matvec=func;
  CDATA.matrix_parameters=data;
  iters[0]=OPTS->max_its;iters[1]=OPTS->max_its;

   int counter;
  // for (counter=0;counter<4*N;counter++){ 
/*     printf("AAA: %6.4e\n",work[counter]); */
  //  work[counter]=0.0; 
  //  } 
  for (counter=0;counter<N;counter++){ 
    mv1[counter]=0.0;
    mv2[counter]=0.0;
    ekt[counter]=0.0;
  }
  ekt[col_id]=1;
  /* Call CG twice */
  //int cg_(n, b, x, work, ldw, iter, resid, matvec, psolve, info)
  tol=OPTS->tolerance;
  cg_(&N,u,mv1,work,&N,iters,&tol,(int(*)())drc_uekt_templates_style_fcall,(int(*)())drc_uekt_identity_preconditioner,&info);
  if(info<0) fprintf(stderr,"Error: CG Failed with Error Code %d\n",info);
  tol=OPTS->tolerance;
  cg_(&N,ekt,mv2,work,&N,&(iters[1]),&tol,(int(*)())drc_uekt_templates_style_fcall,(int(*)())drc_uekt_identity_preconditioner,&info);
  if(info<0) fprintf(stderr,"Error: CG Failed with Error Code %d\n",info);
  OPTS->its_used=iters[0]+iters[1];
    printf("opts used is %i\n",OPTS->its_used);
  /* Calculate the lambdas (eigenvalues) */
  lambda[0]=ddot_(&N,ekt,&one,mv1,&one);
  lambda[1]=lambda[0] - ddot_(&N,ekt,&one,mv2,&one) * ddot_(&N,u,&one,mv1,&one)  / (1+lambda[0]);

  /* Calculate ratio, acceptance */
  (*det_ratio)=(1.0+lambda[0])*(1.0+lambda[1]);
  ///////  printf("%6.4e %6.4e %6.4e\n",lambda[0],lambda[1],*det_ratio);
  //  rv=(*det_ratio) > die_roll;
  printf("-- det=%6.4e die=%6.4e\n",*det_ratio,die_roll);
  
  /* Cleanup */
  //////  free(mv1);free(mv2);free(ekt);free(work);
  //  return (*det_Ratio);
}/* end det_ratio_calculator_uekt_symmetric_vanilla*/


// int det_ratio_calculator_uekt_symmetric_vanilla(int (*func)(double*,double*,int,void*),
// 				       void * data, int N, double * u, int col_id, 
// 				      double die_roll, double* det_ratio,  void *params){

//   /* Allocate Space */
//   double *mv1=(double*)calloc(N,sizeof(double));
//   double *mv2=(double*)calloc(N,sizeof(double));
//   double *ekt=(double*)calloc(N,sizeof(double));ekt[col_id]=1;
//   double *work=(double*)malloc(4*N*sizeof(double));
//   double lambda[2];
//   double tol;
//   int iters[2];
//   int info,rv,one=1;
  
//   /* Parse Parameters */
//   struct drc_uekt_vanilla_parms* OPTS = (struct drc_uekt_vanilla_parms*) params;
//   CDATA.N=N; CDATA.matvec=func;
//   CDATA.matrix_parameters=data;
//   iters[0]=OPTS->max_its;iters[1]=OPTS->max_its;
  
//   /* Call CG twice */
//   //int cg_(n, b, x, work, ldw, iter, resid, matvec, psolve, info)
//   tol=OPTS->tolerance;
//   cg_(&N,u,mv1,work,&N,iters,&tol,drc_uekt_templates_style_fcall,drc_uekt_identity_preconditioner,&info);
//   if(info<0) fprintf(stderr,"Error: CG Failed with Error Code %d\n",info);
//   tol=OPTS->tolerance;
//   cg_(&N,ekt,mv2,work,&N,&(iters[1]),&tol,drc_uekt_templates_style_fcall,drc_uekt_identity_preconditioner,&info);
//   if(info<0) fprintf(stderr,"Error: CG Failed with Error Code %d\n",info);
//   OPTS->its_used=iters[0]+iters[1];
  
//   /* Calculate the lambdas (eigenvalues) */
//   lambda[0]=ddot_(&N,ekt,&one,mv1,&one);
//   lambda[1]=lambda[0] - ddot_(&N,ekt,&one,mv2,&one) * ddot_(&N,u,&one,mv1,&one)  / (1+lambda[0]);

//   /* Calculate ratio, acceptance */
//   (*det_ratio)=(1+lambda[0])*(1+lambda[1]);
//   rv=(*det_ratio) > die_roll;
//   //  printf("-- det=%6.4e die=%6.4e\n",*det_ratio,die_roll);
  
//   /* Cleanup */
//   free(mv1);free(mv2);free(ekt);free(work);
//   return rv;
// }/* end det_ratio_calculator_uekt_symmetric_vanilla*/


// /********************************************************************************************/
// /********************************************************************************************/
// /********************************************************************************************/

// int det_ratio_calculator_uekt_nonsymmetric_col_vanilla(int (*func)(double*,double*,int,void*),
//                                                     void * data, int N, double * u, int col_id, 
//                                                     double die_roll, double* det_ratio,  void *params){
  
//   /* Allocate Space */
//   int restart=MIN(MAX(N/100,100),N-1);
//   double *mv1=(double*)calloc(N,sizeof(double));
//   double *ekt=(double*)calloc(N,sizeof(double));ekt[col_id]=1;
//   double *work=(double*)malloc(N*(restart+4)*sizeof(double));
//   double *hwork=(double*)malloc(N*(restart+2)*sizeof(double));
//   double lambda;
//   double tol;
//   int iters;
//   int info,rv,one=1;
  
//   /* Parse Parameters */
//   struct drc_uekt_vanilla_parms* OPTS = (struct drc_uekt_vanilla_parms*) params;
//   CDATA.N=N; CDATA.matvec=func;
//   CDATA.matrix_parameters=data;
//   iters=OPTS->max_its;
  
//   /* Call GMRES  */
//   //int gmres_(const int * N, const double *b, double * x, const int*restrt, double *work, const int *ldw, double *h, const int *ldh, int *iter, double  *resid,int (*matvec)(), int (*psolve)(), int *info);
//   tol=OPTS->tolerance;
//   gmres_(&N,u,mv1,&restart,work,&N,hwork,&N,&iters,&tol,drc_uekt_templates_style_fcall,drc_uekt_identity_preconditioner,&info);
//   if(info<0) fprintf(stderr,"Error: GMRES Failed with Error Code %d\n",info);
//   OPTS->its_used=iters;
  
//   /* Calculate the lambda (eigenvalue) */
//   lambda=ddot_(&N,ekt,&one,mv1,&one);
  
//   /* Calculate ratio, acceptance */
//   (*det_ratio)=1+lambda;
//   rv=(*det_ratio) > die_roll;
//   //  printf("-- det=%6.4e die=%6.4e\n",*det_ratio,die_roll);
  
//   /* Cleanup */
//   free(mv1);free(ekt);free(work);free(hwork);
//   return rv;
// }/* end det_ratio_calculator_uekt_symmetric_vanilla*/


#endif
#endif
