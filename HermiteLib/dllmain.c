//#include "hermite.h"
//#include "pch.h"

#include <windows.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "hermite.h"
#include "functions.h"


#ifdef __cplusplus
extern "C" 
{
#endif


	herm_params __declspec(dllexport) _HermGen(int8_t funcnum, int8_t linknum, double a, double b, double nu)
	{
		herm_params res = { 0, (int)0, (int)0, (double*)0, (double*)0 };
		res.type = linknum;
		res.param_count = 5 - linknum % 2;
		function f[] = { funcs[funcnum], dfuncs[funcnum] };
		HermGen(f, &res, a, b, nu);


		return res;
	}

	void  __declspec(dllexport) _free(herm_params hp)
	{
		free(hp.A);
		free(hp.X);
	}

	double __declspec(dllexport) _HermiteSpline(const herm_params hp, const double x, int8_t der)
	{
		return HermiteSpline(hp, x, der);
	}

	double __declspec(dllexport) _Func(const int8_t funcnum, const double x, int8_t der)
	{
		return der ? dfuncs[funcnum](x) : funcs[funcnum](x);
	}

	double __declspec(dllexport) _MaxError(const herm_params hp, const int8_t funcnum, const double from, const double to)
	{
		double res = 0.0;
		for (int i = 0; i < 0xFFFF; i++) {
			const double step = (to - from) / (0xFFFF);
			double err = fabs(HermiteSpline(hp, from + step * i, 0) - funcs[funcnum](from + step * i));
			if (res < err)
				res = err;
		}
		return res;
	}
	/*
	double __declspec(dllexport) shit(double shits[], const int32_t count)
	{
		double res = 0.0;
		for (int i = 0; i < count; i++) {
			res += shits[i];
		}
		return res;
	}
	
	int32_t __declspec(dllexport) alloc()
	{
		return (int32_t)malloc(0xFFF);
	}

	int32_t __declspec(dllexport) fre(int32_t ptr)
	{
		free((void*)ptr);
		return 0;
	}
	*/
	BOOL APIENTRY DllMain(HMODULE hModule,
		DWORD ul_reason_for_call,
		LPVOID lpReserved) {
		switch (ul_reason_for_call) {
		case DLL_PROCESS_ATTACH:
			break;
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
		}
		return TRUE;
	}



#ifdef __cplusplus
}
#endif









