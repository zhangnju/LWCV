#ifndef LWCV_VERSION_H
#define LWCV_VERSION_H

struct CPUInfo{
public:
	bool bSupportSSE;
	bool bSupportSSE2;
	bool bSupportSSE3;
	bool bSupportSSE4_1;
	bool bSupportSSE4_2;
	bool bSupportAVX;
	bool bSupportAVX2;
};

CPUInfo CheckCPUInfo();
#endif