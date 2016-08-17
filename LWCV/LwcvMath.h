#ifndef LWCVMATH_H
#define LWCVMATH_H
#include<xmmintrin.h>
#include<emmintrin.h>
#include<immintrin.h>
namespace SSE{
	typedef __m128 v4sf;  // vector of 4 float (sse1)
	typedef __m128i v4si; // vector of 4 int (sse2)
	typedef __m64 v2si;   // vector of 2 int (mmx)

	v4sf lwcv_log_ps(v4sf x);
	v4sf lwcv_exp_ps(v4sf x);
	v4sf lwcv_sin_ps(v4sf x);
	v4sf lwcv_cos_ps(v4sf x);
	void lwcv_sincos_ps(v4sf x, v4sf *s, v4sf *c);
}

namespace AVX{
	typedef __m256 v8sf;
	typedef __m256i v8si;

	v8sf lwcv_log_ps(v8sf x);
	v8sf lwcv_exp_ps(v8sf x);
	v8sf lwcv_sin_ps(v8sf x);
	v8sf lwcv_cos_ps(v8sf x);
	void lwcv_sincos_ps(v8sf x, v8sf *s, v8sf *c);
}

#endif