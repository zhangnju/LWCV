#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include "LWCV.h"

using namespace std;

uint8_t g_rand[UINT16_MAX];
bool InitRand()
{
	for (size_t i = 0, n = UINT16_MAX; i < n; ++i)
		g_rand[i] = ::rand();
	return true;
}
bool g_inited = false;

void FillRandom(uint8_t * data, size_t width, size_t height,uint8_t lo, uint8_t hi)
{
	assert(data!=NULL || width >0 || height >0);

	if (!g_inited)
		g_inited = InitRand();

	bool fast = (lo == 0) && (hi == 255);
	for (size_t row = 0; row < height; ++row)
	{
		ptrdiff_t offset = row*width;
		if (fast)
		{
			for (size_t col = 0; col < width; col += INT16_MAX)
				memcpy(data + offset + col, g_rand + (::rand()&INT16_MAX), std::min<size_t>(INT16_MAX, width - col));
		}
		else
		{
			for (size_t col = 0; col < width; ++col, ++offset)
				data[offset] = lo + ((::rand()&INT16_MAX)*(hi - lo + 1)) / INT16_MAX;
		}
	}
}
extern void TestMath();
int main()
{
	//LwcvVersion();
	TestMath();

	size_t width = 1920;
	size_t height = 1080;
	void *data= LwcvAllocate(width*height*3, 16);
	void *result = LwcvAllocate(width*height * 3, 16);
	FillRandom((uint8_t*)data, width*3, height, 0, 255);
	LwcvBgrToHsv((uint8_t*)data, width, height, width, (uint8_t*)result, width);
	LwcvFree(data);
	LwcvFree(result);
	return 0;
}