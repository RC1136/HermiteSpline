#include "hermite.h"
#include <stdlib.h>
#include <quadmath.h>
//#include <stdio.h>

#define isodd(x) (x & 0b1)

//���������-��������������� ����� � ������� �����������
__float128 PE4(const __float128 x, const __float128 a[4])
{
	return a[0] * powq((__float128)x, a[1] + a[2] * x) * expq(a[3] * x);
}

//������� ���������-��������������� ����� � ������� �����������
__float128 dPE4(const __float128 x, const __float128 a[4])
{
	return PE4(x, a) * (a[2] * (1 + logq(x)) + a[1] / x + a[3]);
}

//���������-��������������� ����� � �'����� �����������
__float128 PE5(const __float128 x, const __float128 a[5])
{
	return a[0] * powq((__float128)x, a[1] + a[2] * x + a[3] * x * x) * expq(a[4] * x);
}

//������� ���������-��������������� ����� � �'����� �����������
__float128 dPE5(const __float128 x, const __float128 a[5])
{
	return PE5(x, a) * ((a[2] + 2 * a[3] * x) * logq((__float128)x) + a[1] / x + a[2] + a[3] * x + a[4]);
}

//������ � count+1 ����������� 
__float128 Polynomial(const __float128 x, const __float128* a, const int count) 
{
	__float128 res = a[0], tmpx = 1.0;
	for (int i = 1; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx;
	}
	return res;
}

//������� �������
__float128 PolynomialDerivative(const __float128 x, const __float128* a, const int count)
{
	__float128 res = a[1], tmpx = 1.0;
	for (int i = 2; i < count; i++) {
		tmpx *= x;
		res += a[i] * tmpx * i;
	}
	return res;
}

//����������� ����� � ������� �����������
__float128 PN4(const __float128 x, const __float128 a[4])
{
	return Polynomial(x, a, 4);
}

//������� ����������� ����� � ������� �����������
__float128 dPN4(const __float128 x, const __float128 a[4])
{
	return PolynomialDerivative(x, a, 4);
}

//����������� ����� � �'����� �����������
__float128 PN5(const __float128 x, const __float128 a[5])
{
	return Polynomial(x, a, 5);
}

//������� ����������� ����� � �'����� �����������
__float128 dPN5(const __float128 x, const __float128 a[5])
{
	return PolynomialDerivative(x, a, 5);
}

//�������� ������ � ����������� hp. ���� derivative == 1, �� ������������ �������
__float128 HermiteSpline(const herm_params hp, const __float128 x, const char derivative)
{
	__float128 (*const link[2][4])(const __float128, const __float128[]) = {
		{PE4, PE5, PN4, PN5}, {dPE4, dPE5, dPN4, dPN5}
	};
	if (x < hp.X128[0] || x > hp.X128[hp.link_count])
		return nanq("");
	int i = 0;
	while (x > hp.X128[++i]);
	return link[derivative][hp.type - 1](x, hp.A128 + (i-1) * hp.param_count);
}

//����'��� ���� ������� �����
int SolveGauss(const __float128** a, const __float128* b, const int count, __float128* out)
{
	//����� ������� �������, ��� ���� ���� ������ ���� ������������
	__float128** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
	*mat = calloc(count * count, sizeof(**mat)); if (*mat == NULL) return -1;
	for (int i = 0; i < count; i++) {
		mat[i] = *mat + i * count;
		for (int j = 0; j < count; j++) {
			mat[i][j] = a[i][j];
		}
	}

	//����� ������� �������, ��� ���� ���� ������ ���� ������������
	__float128* vec = calloc(count, sizeof(*vec)); if (vec == NULL) return -1;
	for (int i = 0; i < count; i++) {
		vec[i] = b[i];
	}

	//������ ��� ������ ������
	for (int row = 0; row < count; row++) {
		//output(mat, vec, count);
		for (int i = count - 1; i > row; i--) {
			//��� � ����������� ��������, �� �������� a[row][row], a[row][j], vec[row] ������ �������� � ������� ����� 
			__float128 l = mat[i][row] / mat[row][row];
			for (int j = count - 1; j >= row; j--) {
				mat[i][j] -= l * mat[row][j];
			}
			vec[i] -= l * vec[row];
		}
	}

	//�������� ���
	for (int i = count - 1; i >= 0; i--) {
		out[i] = vec[i];

		for (int j = count - 1; j > i; j--) {
			out[i] -= mat[i][j] * out[j];	//���������� ������ �������, �� ��� ���� ������� ����� �� ��� ������
		}
		out[i] /= mat[i][i];
	}

	//�������� �� �����
	free(vec);
	free(*mat);
	free(mat);

	return 0;
}

//�������� ��������� ����� PE4
int HermGenPE4(const __float128 f[4], const __float128 x0, const __float128 x1, __float128* out)
{
	//���. notes10(3)
	out[2] = (f[1] / f[0] - (logq(f[0] / f[2]) / (x0 - x1)) - (1. / x0 - logq(x0 / x1) / (x0 - x1)) * (f[1] / f[0] - f[3] / f[2]) / (1. / x0 - 1. / x1)) /
		(1 + logq(x0) - logq(powq(x0, x0) / powq(x1, x1)) / (x0 - x1) - (logq(x0 / x1) / (1 / x0 - 1 / x1)) * (1. / x0 - logq(x0 / x1) / (x0 - x1)));
	out[1] = (f[1] / f[0] - f[3] / f[2] - out[2] * logq(x0 / x1)) / (1. / x0 - 1. / x1);
	//out[3] = (log(f[0] / f[2]) - log(pow(x0, out[2] * x0 + out[1]) / pow(x1, out[2] * x1 + out[1]))) / (x0 - x1);
	out[3] = (logq(f[0] / f[2]) - out[2] * logq(powq(x0, x0) / powq(x1, x1)) - out[1] * logq(x0 / x1)) / (x0 - x1);

	out[0] = f[0] * powq(x0, -(out[2] * x0 + out[1])) * expq(-(out[3] * x0));

/*
#ifndef _DEBUG
	for (int i = 0; i < 4; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif
*/
	return 0;
}

//�������� ��������� ����� PE5
int HermGenPE5(const __float128 f[5], const __float128 x0, const __float128 x2, __float128* out)
{
	//������� �� 8-�� �����
	const __float128 x1 = (x2 + x0) * 0.5;

	//alphas
	const __float128 vec[3] = {
		logq(f[2] / f[0]) / (x1 - x0) - logq(f[3] / f[0]) / (x2 - x0), //alpha1
		f[1] / f[0] - logq(f[2] / f[0]) / (x1 - x0),	//alpha2
		f[4] / f[3] - logq(f[2] / f[0]) / (x1 - x0)	//alpha3
	};

	__float128 matrix[3][3] = {
		{
			logq(x1 / x0) / (x1 - x0) - logq(x2 / x0) / (x2 - x0), //beta1
			logq(powq(x1,x1) / powq(x0,x0)) / (x1 - x0) - logq(powq(x2,x2) / powq(x0,x0)) / (x2 - x0), //gamma1
			logq(powq(x1,x1 * x1) / powq(x0,x0 * x0)) / (x1 - x0) - logq(powq(x2,x2 * x2) / powq(x0,x0 * x0)) / (x2 - x0) //delta1
		},
		{
			1. / x0 - logq(x1 / x0) / (x1 - x0), //beta2
			logq(x0) + 1. - logq(powq(x1,x1) / powq(x0,x0)) / (x1 - x0), //gamma2
			2. * x0 * logq(x0) + x0 - logq(powq(x1,x1 * x1) / powq(x0,x0 * x0)) / (x1 - x0) //delta2
		},
		{
			1. / x2 - logq(x1 / x0) / (x1 - x0), //beta3
			logq(x2) + 1. - logq(powq(x1,x1) / powq(x0,x0)) / (x1 - x0), //gamma3
			2. * x2 * logq(x2) + x2 - logq(powq(x1,x1 * x1) / powq(x0,x0 * x0)) / (x1 - x0) // delta3
		}
	};

	__float128* mat[3] = { &matrix[0][0], &matrix[1][0], &matrix[2][0] };
	SolveGauss((const __float128**)mat, vec, 3, out + 1); // => a[1], a[2], a[3] aka b, c, d

	/*
	out[4] = (log(f[2] / f[0]) -
		log(pow(x1, out[1] + out[2] * x1 + out[3] * x1 * x1) /
			pow(x0, out[1] + out[2] * x0 + out[3] * x0 * x0))
		) / (x1 - x0); // h
	*/
	out[4] = (logq(f[2] / f[0]) - (
			out[1] * logq(x1 / x0) + 
			out[2] * logq(powq(x1, x1) / powq(x0, x0)) + 
			out[3] * logq(powq(x1, x1 * x1) / powq(x0, x0 * x0))
		)
	) / (x1 - x0); // h

	out[0] = f[0] * powq(x0, -(out[1] + out[2] * x0 + out[3] * x0 * x0)) * expq(-(out[4] * x0)); //a
	//out[0] = log(f[0]) - out[1] * log(x0) - out[2] * x0 * log(x0) - out[3] * x0 * x0 * log(x0) - out[4] * x0;
	//out[0] = exp(out[0]);
/*
#ifndef _DEBUG
	for (int i = 0; i < 5; i++)
		if (isnan(out[i]) || isinf(out[i]))
			return -1;
#endif
*/
	return 0;
}

int HermGenPN(const __float128* f, const __float128 x0, const __float128 x2, const int count, __float128* out)
{
	const int odd = isodd(count);
	const __float128 X[2] = { x0, x2 };
	__float128** mat = calloc(count, sizeof(*mat)); if (mat == NULL) return -1;
	*mat = calloc(count * count, sizeof(**mat)); if (*mat == NULL) return -1;
	for (int i = 0; i < count; i++) mat[i] = *mat + i * count;

	const int prod[2][5] = { {1,1,1,1,1}, {0,1,2,3,4} };

	for (int k = 0; k < 2; k++) { //�������� �������� (� �������)
		for (int i = k; i < count ; i++) {	//�������� �����
			//����� ����� ������� mat - �� �������� ��� ������� ��������� a[i]
			for (int j = 0; j < 2; j++) {
				mat[k + j * (2 + odd)][i] = prod[k][i] * powq(X[j], i - k);	//��� ��������� ���������� �� ��, �� mat � �������, � �� ��������
			}
		}
	}
	if (odd) {
		const __float128 x2 = (X[0] + X[1]) * 0.5;
		for (int i = 0; i < count; i++) {
			mat[2][i] = powq(x2, i);
		}
	}
	
	__float128* vec = calloc(count, sizeof(*vec)); if (vec == NULL) return -1;
	for(int i = 0; i < count; i++){
		vec[i] = (__float128)f[i];
	}
	
	SolveGauss((const __float128**)mat, vec, count, out);
	/*
	for (int i = 0; i < 4; i++){
		if (isnan(out[i])){
			return -1;
		}
	}
	*/
	free(vec);
	free(*mat);
	free(mat);
	
	return 0;
}

//�������� ��������� ����� PN4
int HermGenPN4(const __float128 f[4], const __float128 x0, const __float128 x2, __float128* out)
{
	return HermGenPN(&f[0], x0, x2, 4, out);
}

//�������� ��������� ����� PN5
int HermGenPN5(const __float128 f[5], const __float128 x0, const __float128 x2, __float128* out)
{
	return HermGenPN(&f[0], x0, x2, 5, out);
}


__float128 finderr(__float128 (*link)(const __float128, const __float128[]), const __float128 params[], function f, const __float128 from, const __float128 to)
{
	double res = 0.0;
	for (int i = 0; i < 0x7FF; i++) {
		const __float128 step = (to - from) / (0x7FF);
		__float128 err = fabsq(link(from + step * i, &(params[0])) - f(from + step * i));
		if (res < err)
			res = err;
	}
	return res;
}


//�������� ��������� ������� � �������� nu
int HermGen(function _f[], herm_params* hp, const __float128 a, const __float128 b, const __float128 nu)
{
	const __float128 eps = nu * 1e-5;
	int (* const Gen[])(const __float128[], const __float128, const __float128, __float128*) = 
	{ HermGenPE4, HermGenPE5, HermGenPN4, HermGenPN5 };
	__float128 (* const link[4])(const __float128, const __float128[]) = 
	{ PE4, PE5, PN4, PN5 };

	const int odd = isodd(hp->param_count),	//�� �-��� ��������� �������
		linknum = hp->type - 1;				//��� �����
	__float128 x0, x2 = a, delta = b;					//��� ���������� �� ����� ��������

	int errorcode = 0;
	

	struct linkparam{
		__float128 params[5];
		__float128 xright;
		void* next;
	}*cur, *top;		//������ ��������� �������
	top = malloc(sizeof(*top)); if (top == NULL) return -1;
	cur = top;

	int count = 0;	//�-��� �����
	do{
		x0 = x2;
		x2 = fminq(b, x2 + delta * 2);
		delta = x2 - x0;
		++count;
		__float128 f[5];	//�������� �������
		f[0] = _f[0](x0), f[1] = _f[1](x0);
		f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);
		if (odd) f[2] = _f[0]((x2 + x0) * 0.5);
		errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
		if (errorcode)
			return errorcode;
		__float128 nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		if (nu1 < nu) {
			cur->xright = x2;
			break;
		}
		__float128 xl = x0, xr = x2, xprev = x2;
		while (fabsq(nu - nu1) > eps) {
			xprev = x2;
			if (nu < nu1) {
				xr = x2;
				x2 = 0.5 * (x2 + xl);
			}
			else {
				xl = x2;
				x2 = 0.5 * (x2 + xr);
			}
			if (xprev == x2)
				return -1;	
			f[2 + odd] = _f[0](x2), f[3 + odd] = _f[1](x2);
			if (odd) f[2] = _f[0]((x2 + x0) * 0.5);
			/*
			for (int i = 0; i < hp->param_count; i++)
				printf("%lf\n", f[i]);
			putchar('\n');
			*/
			errorcode = Gen[linknum](&f[0], x0, x2, &(cur->params[0]));
			if (errorcode)
				return errorcode;
			nu1 = finderr(link[linknum], cur->params, _f[0], x0, x2);
		}
		/*
		for (int i = 0; i < hp->param_count; i++) {
			printf("A[%d][%d] == %lf\n", count - 1, i, cur->params[i]);
		}
		putchar('\n');
		putchar('\n');
		*/
		cur->xright = x2;
		cur->next = malloc(sizeof(*cur));
		cur = cur->next; if (cur == NULL) return -1;
	} while (x2 < b);

	hp->link_count = count;
	hp->A = calloc(count * hp->param_count, sizeof(*(hp->A))); if (hp->A == NULL) return -1;
	hp->A128 = calloc(count * hp->param_count, sizeof(*(hp->A128))); if (hp->A128 == NULL) return -1;
	hp->X = calloc(count + 1, sizeof(*(hp->X))); if (hp->X == NULL) return -1;
	hp->X128 = calloc(count + 1, sizeof(*(hp->X128))); if (hp->X128 == NULL) return -1;
	struct linkparam* iterator = top;
	hp->X[0] = (double)a;
	hp->X128[0] = a;
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < hp->param_count; j++) {
			hp->A[i * hp->param_count + j] = (double)(iterator->params[j]);
			hp->A128[i * hp->param_count + j] = iterator->params[j];
			//printf("A[%d][%d] == %lf\n", i, j, iterator->params[j]);
		}
		hp->X[i + 1] = (double)iterator->xright;
		hp->X128[i + 1] = iterator->xright;
		void* tmp = iterator;
		iterator = iterator->next;
		free(tmp);
	}
	//putchar('\n');


	return 0;
}
