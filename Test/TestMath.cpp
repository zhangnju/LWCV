#include "LWCV.h"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<immintrin.h>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

typedef __m128 v4sf;
typedef __m256 v8sf;

#define ALIGN16 __declspec(align(16))
typedef ALIGN16 union {
	float f[4];
	int i[4];
	v4sf  v;
}V4SF;

#define ALIGN32 __declspec(align(32))
typedef ALIGN32 union {
	float f[8];
	int i[8];
	v8sf  v;
}V8SF;

int check_sincos_precision(float xmin, float xmax) {
	unsigned nb_trials = 100000;
	cout << "checking sines on [" << xmin*M_PI << " " << xmax*M_PI <<"]"<<endl;


	float max_err_sin_ref = 0, max_err_sin_x = 0;
	float max_err_cos_ref = 0, max_err_cos_x = 0;
	float max_err_sum_sqr_test = 0;
	float max_err_sum_sqr_ref = 0;
	xmin *= M_PI; xmax *= M_PI;
	for (unsigned int i = 0; i < nb_trials; ++i) {
#if 0
		V4SF vx, sin4, cos4, sin4_2, cos4_2;
		vx.f[0] = i*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[1] = (i + .5)*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[2] = ((float)rand() / RAND_MAX)*(xmax - xmin);
		vx.f[3] = (i / 32)*M_PI / ((i % 32) + 1);
		if (vx.f[3] < xmin || vx.f[3] > xmax) vx.f[3] = ((float)rand() / RAND_MAX)*(xmax - xmin);

		
		LwcvSinPS(vx.f,4,sin4.f);
		LwcvCosPS(vx.f,4,cos4.f);
		LwcvSinCosPS(vx.f, 4,sin4_2.f, cos4_2.f);
		
		for (unsigned int j = 0; j < 4; ++j) {
#else
		V8SF vx, sin8, cos8, sin8_2, cos8_2;
		vx.f[0] = i*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[1] = (i + .5)*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[2] = ((float)rand() / RAND_MAX)*(xmax - xmin);
		vx.f[3] = (i / 32)*M_PI / ((i % 32) + 1);
		vx.f[4] = i*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[5] = (i + .5)*(xmax - xmin) / (nb_trials - 1) + xmin;
		vx.f[6] = ((float)rand() / RAND_MAX)*(xmax - xmin);
		vx.f[7] = (i / 32)*M_PI / ((i % 32) + 1);
		if (vx.f[3] < xmin || vx.f[3] > xmax) vx.f[3] = ((float)rand() / RAND_MAX)*(xmax - xmin);
		if (vx.f[7] < xmin || vx.f[7] > xmax) vx.f[3] = ((float)rand() / RAND_MAX)*(xmax - xmin);


		LwcvSinPS(vx.f, 8, sin8.f);
		LwcvCosPS(vx.f, 8, cos8.f);
		LwcvSinCosPS(vx.f, 8, sin8_2.f, cos8_2.f);

		for (unsigned int j = 0; j < 8; ++j) {
#endif
			float x = vx.f[j];
#if 0
			float sin_test = sin4.f[j];
			float cos_test = cos4.f[j];
			if (sin_test != sin4_2.f[j]) {
				cout<<"sin / sincos mismatch at x="<<x<<endl;
				return 1;
			}
			if (cos_test != cos4_2.f[j]) {
				cout<<"cos / sincos mismatch at x="<<x<<endl;
				return 1;
			}
#else
			float sin_test = sin8.f[j];
			float cos_test = cos8.f[j];
			if (sin_test != sin8_2.f[j]) {
				cout << "sin / sincos mismatch at x=" << x << endl;
				return 1;
			}
			if (cos_test != cos8_2.f[j]) {
				cout << "cos / sincos mismatch at x=" << x << endl;
				return 1;
			}
#endif
			float sin_ref = sinf(x);
			float err_sin_ref = fabs(sin_ref - sin_test);
			if (err_sin_ref > max_err_sin_ref) {
				max_err_sin_ref = err_sin_ref;
				max_err_sin_x = x;
			}
			float cos_ref = cosf(x);
			float err_cos_ref = fabs(cos_ref - cos_test);
			if (err_cos_ref > max_err_cos_ref) {
				max_err_cos_ref = err_cos_ref;
				max_err_cos_x = x;
			}
			float err_sum_sqr_test = fabs(1 - cos_test*cos_test - sin_test*sin_test);
			float err_sum_sqr_ref = fabs(1 - cos_ref*cos_ref - sin_ref*sin_ref);
			max_err_sum_sqr_ref = max<float>(max_err_sum_sqr_ref, err_sum_sqr_ref);
			max_err_sum_sqr_test = max<float>(max_err_sum_sqr_test, err_sum_sqr_test);
		}
	}
	
	cout << "deviation of sin(x)^2+cos(x)^2-1: " << max_err_sum_sqr_test << " (ref deviation is " << max_err_sum_sqr_ref << ")" << endl;

	if (max_err_sum_sqr_ref < 2e-7 && max_err_sin_ref < 2e-7 && max_err_cos_ref < 2e-7) {
		cout<<"   ->> precision OK for the sin_ps / cos_ps / sincos_ps <<-"<<endl;
		return 0;
	}
	else {
		cout<<"   WRONG PRECISION !! there is a problem"<<endl;
		return 1;
	}
}

int check_explog_precision(float xmin, float xmax) {
	unsigned nb_trials = 100000;
	cout<<"checking exp/log ["<<xmin<<","<<xmax<<"]"<<endl;

	float max_err_exp_ref = 0, max_err_exp_x = 0;
	float max_err_log_ref = 0, max_err_log_x = 0;
	float max_err_logexp_test = 0;
	float max_err_logexp_ref = 0;
	for (unsigned i = 0; i < nb_trials; ++i) {
#if 0
		V4SF vx, exp4, log4;
		vx.f[0] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[1] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[2] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[3] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		LwcvExpPS(vx.f,4,exp4.f);
		LwcvLogPS(exp4.f,4,log4.f);
		for (unsigned int j = 0; j < 4; ++j) {
			float x = vx.f[j];
			float exp_test = exp4.f[j];
			float log_test = log4.f[j];
#else
		V8SF vx, exp8, log8;
		vx.f[0] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[1] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[2] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[3] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[4] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[5] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[6] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		vx.f[7] = ((float)rand() / RAND_MAX)*(xmax - xmin) + xmin;
		LwcvExpPS(vx.f, 8, exp8.f);
		LwcvLogPS(exp8.f, 8, log8.f);
		for (unsigned int j = 0; j < 8; ++j) {
			float x = vx.f[j];
			float exp_test = exp8.f[j];
			float log_test = log8.f[j];
#endif
			float exp_ref = expf(x);
			float err_exp_ref = fabs(exp_ref - exp_test) / exp_ref;
			if (err_exp_ref > max_err_exp_ref) {
				max_err_exp_ref = err_exp_ref;
				max_err_exp_x = x;
			}

			float log_ref = logf(exp_test);
			float err_log_ref = fabs(log_ref - log_test);
			if (err_log_ref > max_err_log_ref) {
				max_err_log_ref = err_log_ref;
				max_err_log_x = x;
			}
			
			float err_logexp_test = fabs(x - log_test);
			float err_logexp_ref = fabs(x - logf(expf(x)));
			max_err_logexp_ref = max<float>(max_err_logexp_ref, err_logexp_ref);
			max_err_logexp_test = max<float>(max_err_logexp_test, err_logexp_test);
		}
	}
	
	cout << "deviation of x - log(exp(x)): " << max_err_logexp_test << "(ref deviation is " << max_err_logexp_ref << ")" << endl;

	if (max_err_logexp_test < 2e-7 && max_err_exp_ref < 2e-7 && max_err_log_ref < 2e-7) {
		cout<<"   ->> precision OK for the exp_ps / log_ps <<-"<<endl;
		return 0;
	}
	else {
		cout<<"   WRONG PRECISION !! there is a problem"<<endl;
		return 1;
	}
}

void check_sincos_performance()
{
	int niter = 1000000;     
	double t0, t1;

	float xf = 0.75, bmin_f = 0.5, bmax_f = 1.0;
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter*4; ++i) {
		xf = sinf(xf); 
		if (xf < bmin_f)
			xf = bmin_f;
		if (xf > bmax_f)
			xf = bmax_f;
	} 
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "origin sin operation is " << t1 - t0 << "seconds" << endl;

	V4SF bmin, bmax,x;
	bmin.v = _mm_set_ps1(0.5), bmax.v = _mm_set_ps1(1.0);
	x.v = _mm_set_ps1(0.75); 
	t0= (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter; ++i) {
			LwcvSinPS(x.f,4,x.f); x.v = _mm_min_ps(x.v, bmax.v); x.v = _mm_max_ps(x.v, bmin.v);      
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "SSE sin operation is " << t1 - t0 << "seconds" << endl;

	V8SF bmin8, bmax8, x8;
	bmin8.v = _mm256_set1_ps(0.5), bmax8.v = _mm256_set1_ps(1.0);
	x8.v = _mm256_set1_ps(0.75);
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter/2; ++i) {
		LwcvSinPS(x8.f, 8, x8.f); x8.v = _mm256_min_ps(x8.v, bmax8.v); x8.v = _mm256_max_ps(x8.v, bmin8.v);
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "AVX sin operation is " << t1 - t0 << "seconds" << endl;

	xf = 0.75;
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter * 4; ++i) {
		xf = cosf(xf);
		if (xf < bmin_f)
			xf = bmin_f;
		if (xf > bmax_f)
			xf = bmax_f;
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "origin cos operation is " << t1 - t0 << "seconds" << endl;

	x.v = _mm_set_ps1(0.75);
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter; ++i) {
		LwcvCosPS(x.f, 4, x.f); x.v = _mm_min_ps(x.v, bmax.v); x.v = _mm_max_ps(x.v, bmin.v);
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "SSE cos operation is " << t1 - t0 << "seconds" << endl;

	
	x8.v = _mm256_set1_ps(0.75);
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter / 2; ++i) {
		LwcvCosPS(x8.f, 8, x8.f); x8.v = _mm256_min_ps(x8.v, bmax8.v); x8.v = _mm256_max_ps(x8.v, bmin8.v);
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "AVX cos operation is " << t1 - t0 << "seconds" << endl;

}

void check_explog_performance()
{
	int niter = 1000000;
	double t0, t1;

	float xf = 0.75, bmin_f = 0.5, bmax_f = 1.0;
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter * 4; ++i) {
		xf = expf(xf);
		xf = logf(xf);
		if (xf < bmin_f)
			xf = bmin_f;
		if (xf > bmax_f)
			xf = bmax_f;
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "origin exp/log operation is " << t1 - t0 << "seconds" << endl;

	V4SF bmin, bmax, x;
	bmin.v = _mm_set_ps1(0.5), bmax.v = _mm_set_ps1(1.0);
	x.v = _mm_set_ps1(0.75);
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter; ++i) {
		LwcvExpPS(x.f, 4, x.f); LwcvLogPS(x.f, 4, x.f);
		x.v = _mm_min_ps(x.v, bmax.v); x.v = _mm_max_ps(x.v, bmin.v);
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "SSE exp/log operation is " << t1 - t0 << "seconds" << endl;

	V8SF bmin8, bmax8, x8;
	bmin8.v = _mm256_set1_ps(0.5), bmax8.v = _mm256_set1_ps(1.0);
	x8.v = _mm256_set1_ps(0.75);
	t0 = (double)clock() / (double)CLOCKS_PER_SEC;
	for (int i = 0; i < niter / 2; ++i) {
		LwcvExpPS(x8.f, 8, x8.f); LwcvLogPS(x8.f, 8, x8.f);
		x8.v = _mm256_min_ps(x8.v, bmax8.v); x8.v = _mm256_max_ps(x8.v, bmin8.v);
	}
	t1 = (double)clock() / (double)CLOCKS_PER_SEC;
	cout << "AVX exp/log operation is " << t1 - t0 << "seconds" << endl;
}
void TestMath()
{
	int err = 0;
	err += check_sincos_precision(0., 1.0);
	err += check_sincos_precision(-1000, 1000);
	err += check_explog_precision(-60, 60);

	check_sincos_performance();
    check_explog_performance();

	return;
}
