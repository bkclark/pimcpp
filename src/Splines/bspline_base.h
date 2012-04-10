#ifndef BSPLINE_BASE_H
#define BSPLINE_BASE_H

// Conventions:
// Postfixes:  
// s:  single precision real
// d:  double precision real
// c:  single precision complex
// z:  double precision complex

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Basic type declarations               ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

typedef enum { PERIODIC, DERIV1, DERIV2, FLAT, NATURAL } bc_code;
typedef enum { U1D, U2D, U3D, NU1D, NU2D, NU3D } spline_code;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal, rVal;
} BCtype_s;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal, rVal;
} BCtype_d;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_c;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_z;


typedef struct
{
  double start, end;
  int num;

  // private
  double delta, delta_inv;
} Ugrid;

typedef struct
{
  spline_code code;
} Bspline;



#endif
