#ifndef BSPLINE_H
#define BSPLINE_H

#include "bspline_base.h"
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////           Bspline structure definitions            ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
#ifdef __SSE2__
#include "bspline_structs_sse.h"
#include "bspline_eval_sse_s.h"
#include "bspline_eval_sse_c.h"
#include "bspline_eval_sse_d.h"
#include "bspline_eval_sse_z.h"
#else
#include "bspline_structs_std.h"
#include "bspline_eval_std_s.h"
#include "bspline_eval_std_c.h"
#include "bspline_eval_std_d.h"
#include "bspline_eval_std_z.h"
#endif

#include "bspline_create.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////            Spline evaluation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// // Value only
// inline void
// eval_UBspline_1d_s     (UBspline_1d_s * restrict spline, double x, 
// 			float* restrict val);
// // Value and gradient
// inline void
// eval_UBspline_1d_s_vg  (UBspline_1d_s* restrict spline, double x, 
// 			float* restrict val, float* restrict grad);
// // Value, gradient, and Laplacian
// inline void
// eval_UBspline_1d_s_vgl (UBspline_1d_s* restrict spline, double x, 
// 			float* restrict val, float* restrict grad, float* restrict lapl);
// // Value, gradient, and Hessian
// inline void
// eval_UBspline_1d_s_vgh (UBspline_1d_s* restrict spline, double x, 
// 			float *val, float *grad, float *hess);

// inline void
// eval_UBspline_2d_s     (UBspline_2d_s* restrict spline, 
// 			double x, double y,
// 			float* restrict val);
// inline void
// eval_UBspline_2d_s_vg  (UBspline_2d_s* restrict spline, 
// 			double x, double y, 
// 		        float* restrict val, float* restrict grad);
// inline void
// eval_UBspline_2d_s_vgl (UBspline_2d_s* restrict spline, 
// 			double x, double y,
// 			float* restrict val, float* restrict grad, 
// 			float* restrict lapl);
// inline void
// eval_UBspline_2d_s_vgh (UBspline_2d_s* restrict spline, 
// 			double x, double y,
// 			float* restrict val, float* restrict grad, 
// 			float* restrict hess);

// inline void
// eval_UBspline_3d_s     (UBspline_3d_s* restrict spline, 
// 			double x, double y, double z,
// 			float* restrict val);
// inline void
// eval_UBspline_3d_s_vg  (UBspline_3d_s* restrict spline, 
// 			double x, double y, double z,
// 			float* restrict val, float* restrict grad);
// inline void
// eval_UBspline_3d_s_vgl (UBspline_3d_s* restrict spline,
// 			double x, double y, double z,
// 			float* restrict val, float* restrict grad, 
// 			float* restrict lapl);
// inline void
// eval_UBspline_3d_s_vgh (UBspline_3d_s* restrict spline, 
// 			double x, double y, double z,
// 			float* restrict val, float* restrict grad, 
// 			float* restrict hess);

// Similarly for the rest of the types.

void
destroy_Bspline (Bspline *ptr);

#endif
