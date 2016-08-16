// LWCV.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "LWCV.h"
#include <malloc.h>
#include <new>


#include "LWCVVersion.h"
#include "LWCVConversion.h"

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

