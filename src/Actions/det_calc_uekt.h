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
#ifdef ORDER_N_FERMIONS


#ifndef _DET_CALC_UEKT_H_
#define _DET_CALC_UEKT_H_

struct drc_uekt_vanilla_parms{
  double tolerance;
  int max_its;
  int its_used;
};


  int testme();
  //  int drc_uekt_identity_preconditioner(double* x, double*y);
    
  //  int det_ratio_calculator_uekt_symmetric_vanilla(int (*func)(double*,double*,int,void*),
  //						  void * data, int N, double * u, int col_id, 
  //						  double die_roll, double* det_ratio,  void *params);
  
  
  void det_ratio_calculator_uekt_symmetric_vanilla_value(int (*func)(double*,double*,int,void*),
							 void * data, int N, double * u, int col_id, 
							 double die_roll, double* det_ratio,  void *params);
  
  //  int drc_uekt_templates_style_fcall(double *alpha, double* x, double *beta, double*y);

  
  //  int det_ratio_calculator_uekt_nonsymmetric_col_vanilla(int (*func)(double*,double*,int,void*),
  //							 void * data, int N, double * u, int col_id, 
  //							 double die_roll, double* det_ratio,  void *params);



#endif
 #endif

