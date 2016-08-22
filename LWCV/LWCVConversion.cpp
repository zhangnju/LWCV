#include "stdafx.h"
#include "LWCVConversion.h"
#include <assert.h>

const float KF_255_DIV_6 = 255.0f / 6.0f;

const int Y_ADJUST = 16;
const int UV_ADJUST = 128;
const int YUV_TO_BGR_AVERAGING_SHIFT = 13;
const int YUV_TO_BGR_ROUND_TERM = 1 << (YUV_TO_BGR_AVERAGING_SHIFT - 1);
const int Y_TO_RGB_WEIGHT = int(1.164*(1 << YUV_TO_BGR_AVERAGING_SHIFT) + 0.5);
const int U_TO_BLUE_WEIGHT = int(2.018*(1 << YUV_TO_BGR_AVERAGING_SHIFT) + 0.5);
const int U_TO_GREEN_WEIGHT = -int(0.391*(1 << YUV_TO_BGR_AVERAGING_SHIFT) + 0.5);
const int V_TO_GREEN_WEIGHT = -int(0.813*(1 << YUV_TO_BGR_AVERAGING_SHIFT) + 0.5);
const int V_TO_RED_WEIGHT = int(1.596*(1 << YUV_TO_BGR_AVERAGING_SHIFT) + 0.5);

const int BGR_TO_YUV_AVERAGING_SHIFT = 14;
const int BGR_TO_YUV_ROUND_TERM = 1 << (BGR_TO_YUV_AVERAGING_SHIFT - 1);
const int BLUE_TO_Y_WEIGHT = int(0.098*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int GREEN_TO_Y_WEIGHT = int(0.504*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int RED_TO_Y_WEIGHT = int(0.257*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int BLUE_TO_U_WEIGHT = int(0.439*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int GREEN_TO_U_WEIGHT = -int(0.291*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int RED_TO_U_WEIGHT = -int(0.148*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int BLUE_TO_V_WEIGHT = -int(0.071*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int GREEN_TO_V_WEIGHT = -int(0.368*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);
const int RED_TO_V_WEIGHT = int(0.439*(1 << BGR_TO_YUV_AVERAGING_SHIFT) + 0.5);

namespace Base
{
	__forceinline int Min(int a, int b)
	{
		return a < b ? a : b;
	}

	__forceinline int Max(int a, int b)
	{
		return a > b ? a : b;
	}

	__forceinline int RestrictRange(int value, int min = 0, int max = 255)
	{
		return Max(min, Min(max, value));
	}

	__forceinline int Average(int a, int b)
	{
		return (a + b + 1) >> 1;
	}

	__forceinline int Average(int a, int b, int c, int d)
	{
		return (a + b + c + d + 2) >> 2;
	}

	void BgrToHsv(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * hsv, size_t hsvStride)
	{
		for (int row = 0; row < (int)height; ++row)
		{
			const uint8_t * pBgr = bgr + row*bgrStride;
			uint8_t * pHsv = hsv + row*hsvStride;
			for (const uint8_t * pBgrEnd = pBgr + width * 3; pBgr < pBgrEnd; pBgr += 3, pHsv += 3)
			{
				int max = Max(pBgr[2], Max(pBgr[1], pBgr[0]));
				int min = Min(pBgr[2], Min(pBgr[1], pBgr[0]));
				int range = max - min;

				if (range)
				{
					int dividend;

					if (pBgr[2] == max)
						dividend = pBgr[1] - pBgr[0] + 6 * range;
					else if (pBgr[1] == max)
						dividend = pBgr[0] - pBgr[2] + 2 * range;
					else
						dividend = pBgr[2] - pBgr[1] + 4 * range;

					hsv[0] = int(KF_255_DIV_6*dividend / range);
				}
				else
					hsv[0] = 0;

				hsv[1] = max ? 255 * range / max : 0;

				hsv[2] = max;
			}
		}
	}
	__forceinline int BgrToY(int blue, int green, int red)
	{
		return RestrictRange(((BLUE_TO_Y_WEIGHT*blue + GREEN_TO_Y_WEIGHT*green + RED_TO_Y_WEIGHT*red +
			BGR_TO_YUV_ROUND_TERM) >> BGR_TO_YUV_AVERAGING_SHIFT) + Y_ADJUST);
	}

	__forceinline int BgrToU(int blue, int green, int red)
	{
		return RestrictRange(((BLUE_TO_U_WEIGHT*blue + GREEN_TO_U_WEIGHT*green + RED_TO_U_WEIGHT*red +
			BGR_TO_YUV_ROUND_TERM) >> BGR_TO_YUV_AVERAGING_SHIFT) + UV_ADJUST);
	}

	__forceinline int BgrToV(int blue, int green, int red)
	{
		return RestrictRange(((BLUE_TO_V_WEIGHT*blue + GREEN_TO_V_WEIGHT*green + RED_TO_V_WEIGHT*red +
			BGR_TO_YUV_ROUND_TERM) >> BGR_TO_YUV_AVERAGING_SHIFT) + UV_ADJUST);
	}
	__forceinline void BgrToYuv420p(const uint8_t * bgr0, size_t bgrStride, uint8_t * y0, size_t yStride, uint8_t * u, uint8_t * v)
	{
		const uint8_t * bgr1 = bgr0 + bgrStride;
		uint8_t * y1 = y0 + yStride;

		y0[0] = BgrToY(bgr0[0], bgr0[1], bgr0[2]);
		y0[1] = BgrToY(bgr0[3], bgr0[4], bgr0[5]);
		y1[0] = BgrToY(bgr1[0], bgr1[1], bgr1[2]);
		y1[1] = BgrToY(bgr1[3], bgr1[4], bgr1[5]);

		int blue = Average(bgr0[0], bgr0[3], bgr1[0], bgr1[3]);
		int green = Average(bgr0[1], bgr0[4], bgr1[1], bgr1[4]);
		int red = Average(bgr0[2], bgr0[5], bgr1[2], bgr1[5]);

		u[0] = BgrToU(blue, green, red);
		v[0] = BgrToV(blue, green, red);
	}

	void BgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
	{
		assert((width % 2 == 0) && (height % 2 == 0) && (width >= 2) && (height >= 2));

		for (size_t row = 0; row < height; row += 2)
		{
			for (size_t colUV = 0, colY = 0, colBgr = 0; colY < width; colY += 2, colUV++, colBgr += 6)
			{
				BgrToYuv420p(bgr + colBgr, bgrStride, y + colY, yStride, u + colUV, v + colUV);
			}
			y += 2 * yStride;
			u += uStride;
			v += vStride;
			bgr += 2 * bgrStride;
		}
	}

	__forceinline void BgrToYuv422p(const uint8_t * bgr, uint8_t * y, uint8_t * u, uint8_t * v)
	{
		y[0] = BgrToY(bgr[0], bgr[1], bgr[2]);
		y[1] = BgrToY(bgr[3], bgr[4], bgr[5]);

		int blue = Average(bgr[0], bgr[3]);
		int green = Average(bgr[1], bgr[4]);
		int red = Average(bgr[2], bgr[5]);

		u[0] = BgrToU(blue, green, red);
		v[0] = BgrToV(blue, green, red);
	}

	void BgrToYuv422p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
	{
		assert((width % 2 == 0) && (width >= 2));

		for (size_t row = 0; row < height; ++row)
		{
			for (size_t colUV = 0, colY = 0, colBgr = 0; colY < width; colY += 2, colUV++, colBgr += 6)
				BgrToYuv422p(bgr + colBgr, y + colY, u + colUV, v + colUV);
			y += yStride;
			u += uStride;
			v += vStride;
			bgr += bgrStride;
		}
	}

	__forceinline void BgrToYuv444p(const uint8_t * bgr, uint8_t * y, uint8_t * u, uint8_t * v)
	{
		const int blue = bgr[0], green = bgr[1], red = bgr[2];
		y[0] = BgrToY(blue, green, red);
		u[0] = BgrToU(blue, green, red);
		v[0] = BgrToV(blue, green, red);
	}

	void BgrToYuv444p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride, uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
	{
		for (size_t row = 0; row < height; ++row)
		{
			for (size_t col = 0, colBgr = 0; col < width; ++col, colBgr += 3)
				BgrToYuv444p(bgr + colBgr, y + col, u + col, v + col);
			y += yStride;
			u += uStride;
			v += vStride;
			bgr += bgrStride;
		}
	}

	__forceinline int YuvToBlue(int y, int u)
	{
		return RestrictRange((Y_TO_RGB_WEIGHT*(y - Y_ADJUST) + U_TO_BLUE_WEIGHT*(u - UV_ADJUST) +
			YUV_TO_BGR_ROUND_TERM) >> YUV_TO_BGR_AVERAGING_SHIFT);
	}

	__forceinline int YuvToGreen(int y, int u, int v)
	{
		return RestrictRange((Y_TO_RGB_WEIGHT*(y - Y_ADJUST) + U_TO_GREEN_WEIGHT*(u - UV_ADJUST) +
			V_TO_GREEN_WEIGHT*(v - UV_ADJUST) + YUV_TO_BGR_ROUND_TERM) >> YUV_TO_BGR_AVERAGING_SHIFT);
	}

	__forceinline int YuvToRed(int y, int v)
	{
		return RestrictRange((Y_TO_RGB_WEIGHT*(y - Y_ADJUST) + V_TO_RED_WEIGHT*(v - UV_ADJUST) +
			YUV_TO_BGR_ROUND_TERM) >> YUV_TO_BGR_AVERAGING_SHIFT);
	}

	__forceinline void YuvToBgr(int y, int u, int v, uint8_t * bgr)
	{
		bgr[0] = YuvToBlue(y, u);
		bgr[1] = YuvToGreen(y, u, v);
		bgr[2] = YuvToRed(y, v);
	}

	__forceinline void Yuv422pToBgr(const uint8_t *y, int u, int v, uint8_t * bgr)
	{
		YuvToBgr(y[0], u, v, bgr);
		YuvToBgr(y[1], u, v, bgr + 3);
	}

	void Yuv420pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride)
	{
		assert((width % 2 == 0) && (height % 2 == 0) && (width >= 2) && (height >= 2));

		for (size_t row = 0; row < height; row += 2)
		{
			for (size_t colUV = 0, colY = 0, colBgr = 0; colY < width; colY += 2, colUV++, colBgr += 6)
			{
				int u_ = u[colUV];
				int v_ = v[colUV];
				Yuv422pToBgr(y + colY, u_, v_, bgr + colBgr);
				Yuv422pToBgr(y + yStride + colY, u_, v_, bgr + bgrStride + colBgr);
			}
			y += 2 * yStride;
			u += uStride;
			v += vStride;
			bgr += 2 * bgrStride;
		}
	}

	void Yuv422pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride)
	{
		assert((width % 2 == 0) && (width >= 2));

		for (size_t row = 0; row < height; ++row)
		{
			for (size_t colUV = 0, colY = 0, colBgr = 0; colY < width; colY += 2, colUV++, colBgr += 6)
				Yuv422pToBgr(y + colY, u[colUV], v[colUV], bgr + colBgr);
			y += yStride;
			u += uStride;
			v += vStride;
			bgr += bgrStride;
		}
	}

	void Yuv444pToBgr(const uint8_t * y, size_t yStride, const uint8_t * u, size_t uStride, const uint8_t * v, size_t vStride,
		size_t width, size_t height, uint8_t * bgr, size_t bgrStride)
	{
		for (size_t row = 0; row < height; ++row)
		{
			for (size_t col = 0, colBgr = 0; col < width; col++, colBgr += 3)
				YuvToBgr(y[col], u[col], v[col], bgr + colBgr);
			y += yStride;
			u += uStride;
			v += vStride;
			bgr += bgrStride;
		}
	}

}


namespace SSE{
	const __m128i K8_SHUFFLE_BGR0_TO_BLUE = _mm_set_epi8(0x0, 0x3, 0x6, 0x9, 0xC, 0xF, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR1_TO_BLUE = _mm_set_epi8(-1, -1, -1, -1, -1, -1, 0x2, 0x5, 0x8, 0xB, 0xE, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR2_TO_BLUE = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0x1, 0x4, 0x7, 0xA, 0xD);

	const __m128i K8_SHUFFLE_BGR0_TO_GREEN = _mm_set_epi8(0x1, 0x4, 0x7, 0xA, 0xD, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR1_TO_GREEN = _mm_set_epi8(-1, -1, -1, -1, -1, 0x0, 0x3, 0x6, 0x9, 0xC, 0xF, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR2_TO_GREEN = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0x2, 0x5, 0x8, 0xB, 0xE);

	const __m128i K8_SHUFFLE_BGR0_TO_RED = _mm_set_epi8(0x2, 0x5, 0x8, 0xB, 0xE, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR1_TO_RED = _mm_set_epi8(-1, -1, -1, -1, -1, 0x1, 0x4, 0x7, 0xA, 0xD, -1, -1, -1, -1, -1, -1);
	const __m128i K8_SHUFFLE_BGR2_TO_RED = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0x0, 0x3, 0x6, 0x9, 0xC, 0xF);

    __forceinline void LoadBgr(const __m128i * p, __m128i & blue, __m128i & green, __m128i & red)
	{
		__m128i bgr[3];
		bgr[0] = _mm_load_si128(p + 0);
		bgr[1] = _mm_load_si128(p + 1);
		bgr[2] = _mm_load_si128(p + 2);
		
		//get blue data
		blue  = _mm_or_si128(_mm_shuffle_epi8(bgr[0], K8_SHUFFLE_BGR0_TO_BLUE), 
			         _mm_or_si128(_mm_shuffle_epi8(bgr[1], K8_SHUFFLE_BGR1_TO_GREEN),
					              _mm_shuffle_epi8(bgr[2], K8_SHUFFLE_BGR2_TO_BLUE)));

		green = _mm_or_si128(_mm_shuffle_epi8(bgr[0], K8_SHUFFLE_BGR0_TO_GREEN),
			         _mm_or_si128(_mm_shuffle_epi8(bgr[1], K8_SHUFFLE_BGR1_TO_GREEN),
			                      _mm_shuffle_epi8(bgr[2], K8_SHUFFLE_BGR2_TO_GREEN)));

		red   =  _mm_or_si128(_mm_shuffle_epi8(bgr[0], K8_SHUFFLE_BGR0_TO_RED),
			         _mm_or_si128(_mm_shuffle_epi8(bgr[1], K8_SHUFFLE_BGR1_TO_RED),
			                      _mm_shuffle_epi8(bgr[2], K8_SHUFFLE_BGR2_TO_RED)));
	}

	__forceinline __m128i Average16(const __m128i & s0, const __m128i & s1)
	{
		return _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(s0, _mm_set1_epi8(0x01)), _mm_maddubs_epi16(s1, _mm_set1_epi8(0x01))), _mm_set1_epi16(0x0002)), 2);
	}

	__forceinline __m128i BgrToY8(__m128i b8, __m128i g8, __m128i r8)
	{
		
		//pack the low and high 8 bit into 32bit
		__m128i b03_r03 = _mm_unpacklo_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b0-b3,r0-r3
		__m128i g03_01  = _mm_unpacklo_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g0-g3
		__m128i b47_r47 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b4-b7,r4-r7
		__m128i g47_01 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g4-g7
		__m128i b811_r811 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b8-b11,r8-r11
		__m128i g811_01 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g8-g11
		__m128i b1215_r1215 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b11-b15,r11-r15
		__m128i g1215_01 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g11-g15

		//turn 32bit RGB data into Y data 
		__m128i y_03 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b03_r03, _mm_set1_epi32(BLUE_TO_Y_WEIGHT | (RED_TO_Y_WEIGHT << 16))),
			                                        _mm_madd_epi16(g03_01, _mm_set1_epi32(GREEN_TO_Y_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
									    BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i y_47 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b47_r47, _mm_set1_epi32(BLUE_TO_Y_WEIGHT | (RED_TO_Y_WEIGHT << 16))),
			                                        _mm_madd_epi16(g47_01, _mm_set1_epi32(GREEN_TO_Y_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                            BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i y_811 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b811_r811, _mm_set1_epi32(BLUE_TO_Y_WEIGHT | (RED_TO_Y_WEIGHT << 16))),
			                                        _mm_madd_epi16(g811_01, _mm_set1_epi32(GREEN_TO_Y_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                           BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i y_1215 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b1215_r1215, _mm_set1_epi32(BLUE_TO_Y_WEIGHT | (RED_TO_Y_WEIGHT << 16))),
			                                          _mm_madd_epi16(g1215_01, _mm_set1_epi32(GREEN_TO_Y_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                           BGR_TO_YUV_AVERAGING_SHIFT);

		__m128i result_lo = _mm_add_epi16(_mm_packs_epi32(y_03, y_47), _mm_set1_epi16(Y_ADJUST));
		__m128i result_hi = _mm_add_epi16(_mm_packs_epi32(y_811, y_1215), _mm_set1_epi16(Y_ADJUST));

		return _mm_packus_epi16(_mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_lo, _mm_set1_epi16(0x00))),
			                    _mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_hi, _mm_set1_epi16(0x00))));
	}               

	__forceinline __m128i BgrToU8(__m128i b8, __m128i g8, __m128i r8)
	{
		//pack the low and high 8 bit into 32bit
		__m128i b03_r03 = _mm_unpacklo_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b0-b3,r0-r3
		__m128i g03_01 = _mm_unpacklo_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g0-g3
		__m128i b47_r47 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b4-b7,r4-r7
		__m128i g47_01 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g4-g7
		__m128i b811_r811 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b8-b11,r8-r11
		__m128i g811_01 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g8-g11
		__m128i b1215_r1215 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b11-b15,r11-r15
		__m128i g1215_01 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g11-g15

		//turn 32bit RGB data into Y data 
		__m128i u_03 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b03_r03, _mm_set1_epi32(BLUE_TO_U_WEIGHT | (RED_TO_U_WEIGHT << 16))),
			                                        _mm_madd_epi16(g03_01, _mm_set1_epi32(GREEN_TO_U_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                          BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i u_47 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b47_r47, _mm_set1_epi32(BLUE_TO_U_WEIGHT | (RED_TO_U_WEIGHT << 16))),
			                                        _mm_madd_epi16(g47_01, _mm_set1_epi32(GREEN_TO_U_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                          BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i u_811 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b811_r811, _mm_set1_epi32(BLUE_TO_U_WEIGHT | (RED_TO_U_WEIGHT << 16))),
			                                         _mm_madd_epi16(g811_01, _mm_set1_epi32(GREEN_TO_U_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                          BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i u_1215 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b1215_r1215, _mm_set1_epi32(BLUE_TO_U_WEIGHT | (RED_TO_U_WEIGHT << 16))),
			                                          _mm_madd_epi16(g1215_01, _mm_set1_epi32(GREEN_TO_U_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			                           BGR_TO_YUV_AVERAGING_SHIFT);

		__m128i result_lo = _mm_add_epi16(_mm_packs_epi32(u_03, u_47), _mm_set1_epi16(UV_ADJUST));
		__m128i result_hi = _mm_add_epi16(_mm_packs_epi32(u_811, u_1215), _mm_set1_epi16(UV_ADJUST));

		return _mm_packus_epi16(_mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_lo, _mm_set1_epi16(0x00))),
			_mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_hi, _mm_set1_epi16(0x00))));
	}

	__forceinline __m128i BgrToV8(__m128i b8, __m128i g8, __m128i r8)
	{
		//pack the low and high 8 bit into 32bit
		__m128i b03_r03 = _mm_unpacklo_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b0-b3,r0-r3
		__m128i g03_01 = _mm_unpacklo_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g0-g3
		__m128i b47_r47 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpacklo_epi8(r8, _mm_set1_epi8(0x00)));//b4-b7,r4-r7
		__m128i g47_01 = _mm_unpackhi_epi16(_mm_unpacklo_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g4-g7
		__m128i b811_r811 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b8-b11,r8-r11
		__m128i g811_01 = _mm_unpacklo_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g8-g11
		__m128i b1215_r1215 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(b8, _mm_set1_epi8(0x00)), _mm_unpackhi_epi8(r8, _mm_set1_epi8(0x00)));//b11-b15,r11-r15
		__m128i g1215_01 = _mm_unpackhi_epi16(_mm_unpackhi_epi8(g8, _mm_set1_epi8(0x00)), _mm_set1_epi16(0x0001));//g11-g15

		//turn 32bit RGB data into Y data 
		__m128i v_03 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b03_r03, _mm_set1_epi32(BLUE_TO_V_WEIGHT | (RED_TO_V_WEIGHT << 16))),
			_mm_madd_epi16(g03_01, _mm_set1_epi32(GREEN_TO_V_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i v_47 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b47_r47, _mm_set1_epi32(BLUE_TO_V_WEIGHT | (RED_TO_V_WEIGHT << 16))),
			_mm_madd_epi16(g47_01, _mm_set1_epi32(GREEN_TO_V_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i v_811 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b811_r811, _mm_set1_epi32(BLUE_TO_V_WEIGHT | (RED_TO_V_WEIGHT << 16))),
			_mm_madd_epi16(g811_01, _mm_set1_epi32(GREEN_TO_V_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			BGR_TO_YUV_AVERAGING_SHIFT);
		__m128i v_1215 = _mm_srai_epi32(_mm_add_epi32(_mm_madd_epi16(b1215_r1215, _mm_set1_epi32(BLUE_TO_V_WEIGHT | (RED_TO_V_WEIGHT << 16))),
			_mm_madd_epi16(g1215_01, _mm_set1_epi32(GREEN_TO_V_WEIGHT | (BGR_TO_YUV_ROUND_TERM << 16)))),
			BGR_TO_YUV_AVERAGING_SHIFT);

		__m128i result_lo = _mm_add_epi16(_mm_packs_epi32(v_03, v_47), _mm_set1_epi16(UV_ADJUST));
		__m128i result_hi = _mm_add_epi16(_mm_packs_epi32(v_811, v_1215), _mm_set1_epi16(UV_ADJUST));

		return _mm_packus_epi16(_mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_lo, _mm_set1_epi16(0x00))),
			_mm_min_epi16(_mm_set1_epi16(0x00FF), _mm_max_epi16(result_hi, _mm_set1_epi16(0x00))));
	}

	__forceinline void BgrToYuv420p(const uint8_t * bgr0, size_t bgrStride, uint8_t * y0, size_t yStride, uint8_t * u, uint8_t * v)
	{   
		const uint8_t * bgr1 = bgr0 + bgrStride;
		uint8_t * y1 = y0 + yStride;

		__m128i blue[2][2], green[2][2], red[2][2];

		LoadBgr((__m128i*)bgr0 + 0, blue[0][0], green[0][0], red[0][0]);
		_mm_store_si128((__m128i*)y0 + 0, BgrToY8(blue[0][0], green[0][0], red[0][0]));

		LoadBgr((__m128i*)bgr0 + 3, blue[0][1], green[0][1], red[0][1]);
		_mm_store_si128((__m128i*)y0 + 1, BgrToY8(blue[0][1], green[0][1], red[0][1]));

		LoadBgr((__m128i*)bgr1 + 0, blue[1][0], green[1][0], red[1][0]);
		_mm_store_si128((__m128i*)y1 + 0, BgrToY8(blue[1][0], green[1][0], red[1][0]));

		LoadBgr((__m128i*)bgr1 + 3, blue[1][1], green[1][1], red[1][1]);
		_mm_store_si128((__m128i*)y1 + 1, BgrToY8(blue[1][1], green[1][1], red[1][1]));

		blue[0][0] = Average16(blue[0][0], blue[1][0]);
		blue[0][1] = Average16(blue[0][1], blue[1][1]);
		green[0][0] = Average16(green[0][0], green[1][0]);
		green[0][1] = Average16(green[0][1], green[1][1]);
		red[0][0] = Average16(red[0][0], red[1][0]);
		red[0][1] = Average16(red[0][1], red[1][1]);

		blue[0][0] = Average16(blue[0][0], blue[0][1]);
		green[0][0] = Average16(green[0][0], green[0][1]);
		red[0][0] = Average16(red[0][0], red[0][1]);
		
		_mm_store_si128((__m128i*)u, BgrToU8(blue[0][0], green[0][0], red[0][0]));
		_mm_store_si128((__m128i*)u, BgrToV8(blue[0][0], green[0][0], red[0][0]));
	}

	void BgrToYuv420p(const uint8_t * bgr, size_t width, size_t height, size_t bgrStride, uint8_t * y, size_t yStride,
		uint8_t * u, size_t uStride, uint8_t * v, size_t vStride)
	{
		assert((width % 2 == 0) && (height % 2 == 0) && (width >= 2 *sizeof(__m128i)) && (height >= 2));
		size_t alignedWidth = width & ~(2 * sizeof(__m128i) - 1);
		const size_t A6 = sizeof(__m128i) * 6;
		for (size_t row = 0; row < height; row += 2)
		{
			for (size_t colUV = 0, colY = 0, colBgr = 0; colY < alignedWidth; colY += 2*sizeof(__m128i), colUV += sizeof(__m128i), colBgr += A6)
				BgrToYuv420p(bgr + colBgr, bgrStride, y + colY, yStride, u + colUV, v + colUV);
			y += 2 * yStride;
			u += uStride;
			v += vStride;
			bgr += 2 * bgrStride;
		}
	}
}