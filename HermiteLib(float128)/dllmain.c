//#include "hermite.h"
//#include "pch.h"

#include <windows.h>
#include <stdint.h>
#include <stdlib.h>
#include <quadmath.h>
#include "hermite.h"
#include "functions.h"


herm_params __declspec(dllexport) _HermGen(int8_t funcnum, int8_t linknum, double a, double b, double nu)
{
	herm_params res = { 0, (int)0, (int)0, (double*)0, (double*)0,  (__float128*)0,  (__float128*)0};
	res.type = linknum;
	res.param_count = 5 - linknum % 2;
	function f[] = { funcs[funcnum], dfuncs[funcnum] };
	HermGen(f, &res, (__float128)a, (__float128)b, (__float128)nu);


	return res;
}

void  __declspec(dllexport) _free(herm_params hp)
{
	free(hp.A);
	free(hp.A128);
	free(hp.X);
	free(hp.X128);
}

double __declspec(dllexport) _HermiteSpline(const herm_params hp, const double x, int8_t der)
{
	return (double)HermiteSpline(hp, (__float128)x, der);
}

double __declspec(dllexport) _Func(const int8_t funcnum, const double x, int8_t der)
{
	return (double)(der ? dfuncs[funcnum]((__float128)x) : funcs[funcnum]((__float128)x));
}

double __declspec(dllexport) _MaxError(const herm_params hp, const int8_t funcnum, const double from, const double to)
{
	__float128 res = 0.0;
	for (int i = 0; i < 0xFFFF; i++) {
		const __float128 step = ((__float128)to - (__float128)from) / (0xFFFF);
		__float128 err = fabsq(HermiteSpline(hp, (__float128)from + step * i, 0) - funcs[funcnum]((__float128)from + step * i));
		if (res < err)
			res = err;
	}
	return (double)res;
}








