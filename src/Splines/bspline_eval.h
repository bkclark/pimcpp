#ifndef BSPLINE_EVAL_H
#define BSPLINE_EVAL_H

#ifdef __SSE3__
  #include "bspline_eval_sse.h"
#else
  #include "bspline_eval_std.h"
#endif
#endif
