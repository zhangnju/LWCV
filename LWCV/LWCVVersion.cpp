#include "stdafx.h"
#include "LWCVVersion.h"
#include <iostream>
#include <string>

#include <stdint.h>
#include <intrin.h>

#include <Windows.h>

using namespace std;

struct cpu_x86{
	//  Vendor
	bool Vendor_AMD;
	bool Vendor_Intel;

	//  OS Features
	bool OS_x64;
	bool OS_AVX;
	bool OS_AVX512;

	//  Misc.
	bool HW_MMX;
	bool HW_x64;
	bool HW_ABM;
	bool HW_RDRAND;
	bool HW_BMI1;
	bool HW_BMI2;
	bool HW_ADX;
	bool HW_PREFETCHWT1;
	bool HW_MPX;

	//  SIMD: 128-bit
	bool HW_SSE;
	bool HW_SSE2;
	bool HW_SSE3;
	bool HW_SSSE3;
	bool HW_SSE41;
	bool HW_SSE42;
	bool HW_SSE4a;
	bool HW_AES;
	bool HW_SHA;

	//  SIMD: 256-bit
	bool HW_AVX;
	bool HW_XOP;
	bool HW_FMA3;
	bool HW_FMA4;
	bool HW_AVX2;

	//  SIMD: 512-bit
	bool HW_AVX512_F;
	bool HW_AVX512_PF;
	bool HW_AVX512_ER;
	bool HW_AVX512_CD;
	bool HW_AVX512_VL;
	bool HW_AVX512_BW;
	bool HW_AVX512_DQ;
	bool HW_AVX512_IFMA;
	bool HW_AVX512_VBMI;

public:
	cpu_x86();
	void detect_host();

	void print() const;
	static void print_host();

	static void cpuid(int32_t out[4], int32_t x);
	static std::string get_vendor_string();

private:
	static void print(const char* label, bool yes);

	static bool detect_OS_x64();
	static bool detect_OS_AVX();
	static bool detect_OS_AVX512();
};

cpu_x86::cpu_x86(){
	memset(this, 0, sizeof(*this));
}

void cpu_x86::cpuid(int32_t out[4], int32_t x){
	__cpuidex(out, x, 0);
}
__int64 xgetbv(unsigned int x){
	return _xgetbv(x);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Detect 64-bit - Note that this snippet of code for detecting 64-bit has been copied from MSDN.
typedef BOOL(WINAPI *LPFN_ISWOW64PROCESS) (HANDLE, PBOOL);
BOOL IsWow64()
{
	BOOL bIsWow64 = FALSE;

	LPFN_ISWOW64PROCESS fnIsWow64Process = (LPFN_ISWOW64PROCESS)GetProcAddress(
		GetModuleHandle(TEXT("kernel32")), "IsWow64Process");

	if (NULL != fnIsWow64Process)
	{
		if (!fnIsWow64Process(GetCurrentProcess(), &bIsWow64))
		{
			printf("Error Detecting Operating System.\n");
			printf("Defaulting to 32-bit OS.\n\n");
			bIsWow64 = FALSE;
		}
	}
	return bIsWow64;
}
bool cpu_x86::detect_OS_x64(){
#ifdef _M_X64
	return true;
#else
	return IsWow64() != 0;
#endif
}

bool cpu_x86::detect_OS_AVX(){
	bool avxSupported = false;

	int cpuInfo[4];
	cpuid(cpuInfo, 1);

	bool osUsesXSAVE_XRSTORE = (cpuInfo[2] & (1 << 27)) != 0;
	bool cpuAVXSuport = (cpuInfo[2] & (1 << 28)) != 0;

	if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
	{
		uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
		avxSupported = (xcrFeatureMask & 0x6) == 0x6;
	}

	return avxSupported;
}
bool cpu_x86::detect_OS_AVX512(){
	if (!detect_OS_AVX())
		return false;

	uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
	return (xcrFeatureMask & 0xe6) == 0xe6;
}
std::string cpu_x86::get_vendor_string(){
	int32_t CPUInfo[4];
	char name[13];

	cpuid(CPUInfo, 0);
	memcpy(name + 0, &CPUInfo[1], 4);
	memcpy(name + 4, &CPUInfo[3], 4);
	memcpy(name + 8, &CPUInfo[2], 4);
	name[12] = '\0';

	return name;
}

void cpu_x86::detect_host(){
	//  OS Features
	OS_x64 = detect_OS_x64();
	OS_AVX = detect_OS_AVX();
	OS_AVX512 = detect_OS_AVX512();

	//  Vendor
	std::string vendor(get_vendor_string());
	if (vendor == "GenuineIntel"){
		Vendor_Intel = true;
	}
	else if (vendor == "AuthenticAMD"){
		Vendor_AMD = true;
	}

	int info[4];
	cpuid(info, 0);
	int nIds = info[0];

	cpuid(info, 0x80000000);
	uint32_t nExIds = info[0];

	//  Detect Features
	if (nIds >= 0x00000001){
		cpuid(info, 0x00000001);
		HW_MMX = (info[3] & ((int)1 << 23)) != 0;
		HW_SSE = (info[3] & ((int)1 << 25)) != 0;
		HW_SSE2 = (info[3] & ((int)1 << 26)) != 0;
		HW_SSE3 = (info[2] & ((int)1 << 0)) != 0;

		HW_SSSE3 = (info[2] & ((int)1 << 9)) != 0;
		HW_SSE41 = (info[2] & ((int)1 << 19)) != 0;
		HW_SSE42 = (info[2] & ((int)1 << 20)) != 0;
		HW_AES = (info[2] & ((int)1 << 25)) != 0;

		HW_AVX = (info[2] & ((int)1 << 28)) != 0;
		HW_FMA3 = (info[2] & ((int)1 << 12)) != 0;

		HW_RDRAND = (info[2] & ((int)1 << 30)) != 0;
	}
	if (nIds >= 0x00000007){
		cpuid(info, 0x00000007);
		HW_AVX2 = (info[1] & ((int)1 << 5)) != 0;

		HW_BMI1 = (info[1] & ((int)1 << 3)) != 0;
		HW_BMI2 = (info[1] & ((int)1 << 8)) != 0;
		HW_ADX = (info[1] & ((int)1 << 19)) != 0;
		HW_MPX = (info[1] & ((int)1 << 14)) != 0;
		HW_SHA = (info[1] & ((int)1 << 29)) != 0;
		HW_PREFETCHWT1 = (info[2] & ((int)1 << 0)) != 0;

		HW_AVX512_F = (info[1] & ((int)1 << 16)) != 0;
		HW_AVX512_CD = (info[1] & ((int)1 << 28)) != 0;
		HW_AVX512_PF = (info[1] & ((int)1 << 26)) != 0;
		HW_AVX512_ER = (info[1] & ((int)1 << 27)) != 0;
		HW_AVX512_VL = (info[1] & ((int)1 << 31)) != 0;
		HW_AVX512_BW = (info[1] & ((int)1 << 30)) != 0;
		HW_AVX512_DQ = (info[1] & ((int)1 << 17)) != 0;
		HW_AVX512_IFMA = (info[1] & ((int)1 << 21)) != 0;
		HW_AVX512_VBMI = (info[2] & ((int)1 << 1)) != 0;
	}
	if (nExIds >= 0x80000001){
		cpuid(info, 0x80000001);
		HW_x64 = (info[3] & ((int)1 << 29)) != 0;
		HW_ABM = (info[2] & ((int)1 << 5)) != 0;
		HW_SSE4a = (info[2] & ((int)1 << 6)) != 0;
		HW_FMA4 = (info[2] & ((int)1 << 16)) != 0;
		HW_XOP = (info[2] & ((int)1 << 11)) != 0;
	}
}


CPUInfo CheckCPUInfo()
{
	CPUInfo Info;
	cpu_x86 Features;
	Features.detect_host();
	Info.bSupportSSE = Features.HW_SSE;
	Info.bSupportSSE2 = Features.HW_SSE2;
	Info.bSupportSSE3 = Features.HW_SSE3;
	Info.bSupportSSE4_1 = Features.HW_SSE41;
	Info.bSupportSSE4_2 = Features.HW_SSE42;
	Info.bSupportAVX = Features.HW_AVX;
	Info.bSupportAVX2 = Features.HW_AVX2;

	return Info;
}