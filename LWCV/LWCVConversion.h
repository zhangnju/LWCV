#ifndef LWCV_CONV_H
#define LWCV_CONV_H
#include <stdint.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<immintrin.h>
#include<tmmintrin.h>
/*
Y =  0.257 * R + 0.504 * G + 0.098 * B + 16;
U = -0.148 * R - 0.291 * G + 0.439 * B + 128;
V =  0.439 * R - 0.368 * G - 0.071 * B + 128;

B =  1.164 * (Y - 16) +  2.018 * (U - 128);
G =  1.164 * (Y - 16) -  0.391 * (U - 128) - 0.813 * (V - 128);
R =  1.164 * (Y - 16)                      + 1.596 * (V - 128);
*/
namespace Base
{
	void BgrToHsv(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * hsv, size_t hsvStride);

	void BgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void Yuv420pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv422pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv444pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);
}

namespace SSE
{
	void BgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void Yuv420pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv422pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv444pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);
}

namespace AVX
{
	void BgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void BgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride);

	void Yuv420pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv422pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);

	void Yuv444pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride);
}
#endif