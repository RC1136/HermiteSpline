#include <quadmath.h>
#include "functions.h"
#include "hermite.h"

__float128 testfunc0(const __float128 _x)
{
	//a = 1, b = -3, c = 4, d = 7
	//return 0.0;
	return PE4(_x, (__float128[]) { 1, -3, 4, 7 });
}

__float128 testdf0(const __float128 _x)
{
	//return 0.0;
	return dPE4(_x, (__float128[]) { 1, -3, 4, 7 });
}

__float128 testfunc1(const __float128 _x)
{
	//a = 1, b = -2, c = 5, d = 3, v = 2
	//return 0.0;
	return PE5(_x, (__float128[]) { 1, -2, 5, 3, 2 });
}

__float128 testdf1(const __float128 _x)
{
	//return 0.0;
	return dPE5(_x, (__float128[]) { 1, -2, 5, 3, 2 });
}

__float128 testfunc2(const __float128 _x)
{
	//a = 3, b = 2, c = -4, d = 1
	//return 0.0;
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x;
}

__float128 testdf2(const __float128 _x)
{
	return 2 - 8 * _x + 3 * _x * _x;
}

__float128 testfunc3(const __float128 _x)
{
	//a = 3, b = 2, c = -4, d = 1, h = -1
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x - 1 * _x * _x * _x * _x;
}

__float128 testdf3(const __float128 _x)
{
	return 2 - 8 * _x + 3 * _x * _x - 4 * _x * _x * _x;
}

__float128 testfunc4(const __float128 _x)
{
	return expq(powq(_x - 4, 3) / 42);
}

__float128 testdf4(const __float128 _x)
{
	return expq(powq(_x - 4., 3.) / 42.) * powq(_x - 4., 2.) / 14.;
}

__float128 testfunc5(const __float128 _x)
{
	return sinq(_x);
}

__float128 testdf5(const __float128 _x)
{
	return cosq(_x);
}

__float128 testfunc6(const __float128 _x)
{
	return sqrtq(powq(_x, 3.) + 1.);
}

__float128 testdf6(const __float128 _x)
{
	return (3. / 2.) * ((_x * _x) / sqrtq(_x * _x * _x + 1));
}

__float128 testfunc7(const __float128 _x)
{
	return 1. / (powq(_x, 2.) + 0.4);
}

__float128 testdf7(const __float128 _x)
{
	return (-2 * _x) / powq(_x * _x + 0.4, 2);
}

__float128 testfunc8(const __float128 _x)
{
	return tanq(_x) / (_x + 2.) + 2;
}

__float128 testdf8(const __float128 _x)
{
	return (1 + powq(tanq(_x), 2)) / (_x + 2) - (tanq(_x)) / powq(_x + 2, 2);
}

__float128 testfunc9(const __float128 _x)
{
	return logq(_x);
}

__float128 testdf9(const __float128 _x)
{
	return 1 / _x;
}

__float128 testfunc10(const __float128 _x)
{
	return expq(sinq(_x) - cosq(_x));
}

__float128 testdf10(const __float128 _x)
{
	return (cosq(_x) + sinq(_x)) * testfunc10(_x);
}

__float128 testfunc11(const __float128 _x)
{
	return expq(sinq(_x)) + _x;
}

__float128 testdf11(const __float128 _x)
{
	return cosq(_x) * expq(sinq(_x)) + 1;
}

__float128 testfunc12(const __float128 _x)
{
	return tanq(_x);
}

__float128 testdf12(const __float128 _x)
{
	return 1. / powq(cosq(_x), 2);
}

__float128 testfunc13(const __float128 _x)
{
	return 1. / (_x * _x * _x + 1);
}

__float128 testdf13(const __float128 _x)
{
	return (-3 * _x * _x) / powq(_x * _x * _x + 1, 2);
}



const function funcs[] = { testfunc0, testfunc1, testfunc2, testfunc3, testfunc4, testfunc5, testfunc6, testfunc7, testfunc8, testfunc9, testfunc10, testfunc11, testfunc12, testfunc13 };
const function dfuncs[] = { testdf0, testdf1, testdf2, testdf3, testdf4, testdf5, testdf6, testdf7, testdf8, testdf9, testdf10, testdf11, testdf12, testdf13 };
