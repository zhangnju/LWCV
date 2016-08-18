// LWCV.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "LWCV.h"
#include <assert.h>
#include <malloc.h>
#include <new>


#include "LWCVVersion.h"
#include "LWCVConversion.h"
#include "LwcvMath.h"

static CPUInfo Info = { 0 };
static bool bChecked=false;
LWCV_API void LwcvVersion()
{
	if (bChecked)
		return;

	Info=CheckCPUInfo();
	bChecked = true;
}

LWCV_API void * LwcvAllocate(size_t size, size_t align)
{
	void* ptr;
	try
	{
		ptr = _aligned_malloc(size, align);
	}
	catch (const std::bad_alloc& e)
	{
		ptr = NULL;
	}
	return ptr;
}

LWCV_API void LwcvFree(void * ptr)
{
	_aligned_free(ptr);
}

LWCV_API size_t LwcvAlign(size_t size, size_t align)
{
	return (size + align - 1) & ~(align - 1);
}

LWCV_API void LwcvBgrToHsv(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * hsv, size_t hsvStride)
{
	if (!bChecked)
	{
		LwcvVersion();
	}

	Base::BgrToHsv(bgr,width,height,bgrStride,hsv,hsvStride);
}

LWCV_API void LwcvBgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
{

}

LWCV_API void LwcvBgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
{

}

LWCV_API void LwcvBgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
{

}

LWCV_API void LwcvLogPS(float* data, uint8_t size, float* res)
{
	assert(size%4==0);
	int offset = 0;
	while (size > 0)
	{
		if (size >= 8)
		{
			AVX::v8sf src = _mm256_load_ps(data+offset);
			AVX::v8sf dst = AVX::lwcv_log_ps(src);
			_mm256_store_ps(res+offset, dst);
			size -= 8;
			offset += 8;
		}
		else 
		{
			SSE::v4sf src = _mm_load_ps(data+offset);
			SSE::v4sf dst = SSE::lwcv_log_ps(src);
			_mm_store_ps(res+offset, dst);
			size -= 4;
			offset += 4;
		}
		
	}
}

LWCV_API void LwcvExpPS(float* data, uint8_t size, float* res)
{
	assert(size % 4 == 0);
	int offset = 0;
	while (size > 0)
	{
		if (size >= 8)
		{
			AVX::v8sf src = _mm256_load_ps(data + offset);
			AVX::v8sf dst = AVX::lwcv_exp_ps(src);
			_mm256_store_ps(res + offset, dst);
			size -= 8;
			offset += 8;
		}
		else
		{
			SSE::v4sf src = _mm_load_ps(data + offset);
			SSE::v4sf dst = SSE::lwcv_exp_ps(src);
			_mm_store_ps(res + offset, dst);
			size -= 4;
			offset += 4;
		}

	}
}

LWCV_API void LwcvSinPS(float* data, uint8_t size, float* res)
{
	assert(size % 4 == 0);
	int offset = 0;
	while (size > 0)
	{
		if (size >= 8)
		{
			AVX::v8sf src = _mm256_load_ps(data + offset);
			AVX::v8sf dst = AVX::lwcv_sin_ps(src);
			_mm256_store_ps(res + offset, dst);
			size -= 8;
			offset += 8;
		}
		else
		{
			SSE::v4sf src = _mm_load_ps(data + offset);
			SSE::v4sf dst = SSE::lwcv_sin_ps(src);
			_mm_store_ps(res + offset, dst);
			size -= 4;
			offset += 4;
		}

	}
}

LWCV_API void LwcvCosPS(float* data, uint8_t size, float* res)
{
	assert(size % 4 == 0);
	int offset = 0;
	while (size > 0)
	{
		if (size >= 8)
		{
			AVX::v8sf src = _mm256_load_ps(data + offset);
			AVX::v8sf dst = AVX::lwcv_cos_ps(src);
			_mm256_store_ps(res + offset, dst);
			size -= 8;
			offset += 8;
		}
		else
		{
			SSE::v4sf src = _mm_load_ps(data + offset);
			SSE::v4sf dst = SSE::lwcv_cos_ps(src);
			_mm_store_ps(res + offset, dst);
			size -= 4;
			offset += 4;
		}

	}
}

LWCV_API void LwcvSinCosPS(float* data, uint8_t size,float* res_sin, float* res_cos)
{
	assert(size % 4 == 0);
	int offset = 0;
	while (size > 0)
	{
		if (size >= 8)
		{
			AVX::v8sf src = _mm256_load_ps(data + offset);
			AVX::v8sf dst_sin, dst_cos;
			AVX::lwcv_sincos_ps(src,&dst_sin,&dst_cos);
			_mm256_store_ps(res_sin + offset, dst_sin);
			_mm256_store_ps(res_cos + offset, dst_cos);
			size -= 8;
			offset += 8;
		}
		else
		{
			SSE::v4sf src = _mm_load_ps(data + offset);
			SSE::v4sf dst_sin, dst_cos;
			SSE::lwcv_sincos_ps(src, &dst_sin, &dst_cos);
			_mm_store_ps(res_sin + offset, dst_sin);
			_mm_store_ps(res_cos + offset, dst_cos);
			size -= 4;
			offset += 4;
		}

	}
}





