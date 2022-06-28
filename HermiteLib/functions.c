#include <math.h>
#include "functions.h"
#include "hermite.h"

double testfunc0(const double _x)
{
	//a = 1, b = -3, c = 4, d = 7
	//return 0.0;
	return PE4(_x, (double[]) { 1, -3, 4, 7 });
}

double testdf0(const double _x)
{
	//return 0.0;
	return dPE4(_x, (double[]) { 1, -3, 4, 7 });
}

double testfunc1(const double _x)
{
	//a = 1, b = -2, c = 5, d = 3, v = 2
	//return 0.0;
	return PE5(_x, (double[]) { 1, -2, 5, 3, 2 });
}

double testdf1(const double _x)
{
	//return 0.0;
	return dPE5(_x, (double[]) { 1, -2, 5, 3, 2 });
}

double testfunc2(const double _x)
{
	//a = 3, b = 2, c = -4, d = 1
	//return 0.0;
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x;
}

double testdf2(const double _x)
{
	return 2 - 8 * _x + 3 * _x * _x;
}

double testfunc3(const double _x)
{
	//a = 3, b = 2, c = -4, d = 1, h = -1
	return 3 + 2 * _x - 4 * _x * _x + 1 * _x * _x * _x - 1 * _x * _x * _x * _x;
}

double testdf3(const double _x)
{
	return 2 - 8 * _x + 3 * _x * _x - 4 * _x * _x * _x;
}

double testfunc4(const double _x)
{
	return exp(pow(_x - 4, 3) / 42);
}

double testdf4(const double _x)
{
	return exp(pow(_x - 4., 3.) / 42.) * pow(_x - 4., 2.) / 14.;
}

double testfunc5(const double _x)
{
	return sin(_x);
}

double testdf5(const double _x)
{
	return cos(_x);
}

double testfunc6(const double _x)
{
	return sqrt(pow(_x, 3.) + 1.);
}

double testdf6(const double _x)
{
	return (3. / 2.) * ((_x * _x) / sqrt(_x * _x * _x + 1));
}

double testfunc7(const double _x)
{
	return 1. / (_x * _x + 0.4);
}

double testdf7(const double _x)
{
	return (-2 * _x) / pow(_x * _x + 0.4, 2);
}

double testfunc8(const double _x)
{
	return tan(_x) / (_x + 2.) + 2;
}

double testdf8(const double _x)
{
	return (1 + pow(tan(_x), 2)) / (_x + 2) - (tan(_x)) / pow(_x + 2, 2);
}

double testfunc9(const double _x)
{
	return log(_x);
}

double testdf9(const double _x)
{
	return 1 / _x;
}

double testfunc10(const double _x)
{
	return exp(sin(_x) - cos(_x));
}

double testdf10(const double _x)
{
	return (cos(_x) + sin(_x)) * testfunc10(_x);
}

double testfunc11(const double _x)
{
	return exp(sin(_x)) + _x;
}

double testdf11(const double _x)
{
	return cos(_x) * exp(sin(_x)) + 1;
}

double testfunc12(const double _x)
{
	return tan(_x);
}

double testdf12(const double _x)
{
	return 1. / pow(cos(_x), 2);
}

double testfunc13(const double _x)
{
	return 1. / (_x * _x * _x + 1);
}

double testdf13(const double _x)
{
	return (-3 * _x * _x) / pow(_x * _x * _x + 1, 2);
}



const function funcs[] = { testfunc0, testfunc1, testfunc2, testfunc3, testfunc4, testfunc5, testfunc6, testfunc7, testfunc8, testfunc9, testfunc10, testfunc11, testfunc12, testfunc13 };
const function dfuncs[] = { testdf0, testdf1, testdf2, testdf3, testdf4, testdf5, testdf6, testdf7, testdf8, testdf9, testdf10, testdf11, testdf12, testdf13 };
