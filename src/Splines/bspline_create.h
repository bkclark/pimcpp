#ifndef BSPLINE_CREATE_H
#define BSPLINE_CREATE_H

#include "bspline_base.h"

#ifdef __SSE2__
#include "bspline_structs_sse.h"
#else
#include "bspline_structs_std.h"
#endif


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Spline creation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/////////////////////////////////////
// Uniform, single precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_s *
create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_s *
create_UBspline_2d_s (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_s   xBC, BCtype_s   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_s *
create_UBspline_3d_s (Ugrid x_grid,   Ugrid y_grid,   Ugrid z_grid,
		      BCtype_s  xBC,  BCtype_s   yBC, BCtype_s   zBC,
		      float *data);


/////////////////////////////////////
// Uniform, double precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_d *
create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_d *
create_UBspline_2d_d (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_d   xBC, BCtype_d   yBC,
		      double *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_d *
create_UBspline_3d_d (Ugrid x_grid,   Ugrid   y_grid,   Ugrid z_grid,
		      BCtype_d  xBC,  BCtype_d   yBC, BCtype_d   zBC,
		      double *data);


///////////////////////////////////////
// Uniform, single precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_c *
create_UBspline_1d_c (Ugrid x_grid, BCtype_c xBC, complex_float *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_c *
create_UBspline_2d_c (Ugrid   x_grid, Ugrid   y_grid,
		      BCtype_c   xBC, BCtype_c   yBC,
		      complex_float *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_c *
create_UBspline_3d_c (Ugrid  x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype_c  xBC, BCtype_c yBC, BCtype_c zBC,
		      complex_float *data);

 
///////////////////////////////////////
// Uniform, double precision, complex//
///////////////////////////////////////
// Create 1D uniform double-precision, complex Bspline
UBspline_1d_z *
create_UBspline_1d_z (Ugrid x_grid, BCtype_z xBC, complex_double *data);

// Create 2D uniform double-precision, complex Bspline
UBspline_2d_z *
create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
		      BCtype_z   xBC, BCtype_z   yBC,
		      complex_double *data);

// Create 3D uniform double-precision, complex Bspline
UBspline_3d_z *
create_UBspline_3d_z (Ugrid  x_grid, Ugrid   y_grid, Ugrid z_grid,
		      BCtype_z  xBC, BCtype_z   yBC, BCtype_z zBC,
		      complex_double *data);


#endif
