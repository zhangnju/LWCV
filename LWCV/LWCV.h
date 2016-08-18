// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the LWCV_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// LWCV_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

#ifdef LWCV_EXPORTS
#define LWCV_API __declspec(dllexport)
#else
#define LWCV_API __declspec(dllimport)
#endif

#ifndef LWCV_H
#define LWCV_H

#include <stdint.h>

LWCV_API void LwcvVersion();

LWCV_API void * LwcvAllocate(size_t size, size_t align);

LWCV_API void LwcvFree(void * p);

LWCV_API size_t LwcvAlign(size_t size, size_t align);

LWCV_API void LwcvBgrToHsv(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * hsv, size_t hsvStride);

LWCV_API void LwcvBgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

LWCV_API void LwcvBgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

LWCV_API void LwcvBgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

LWCV_API void LwcvLogPS(float* data, uint8_t size, float* res);

LWCV_API void LwcvExpPS(float* data, uint8_t size, float* res);

LWCV_API void LwcvSinPS(float* data, uint8_t size, float* res);

LWCV_API void LwcvCosPS(float* data, uint8_t size, float* res);

LWCV_API void LwcvSinCosPS(float* data, uint8_t size, float* res_sin, float* res_cos);

#endif
