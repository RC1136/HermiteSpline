

//Вигляд ланки
enum linktype {
	powexp4 = 1,
	powexp5 = 2,
	poly4 = 3,
	poly5 = 4
};

typedef struct {
	enum linktype type;
	int param_count; //к-сть параметрів у одній ланці
	int link_count;	 //кількість ланок
	double* A;		 //параметри сплайна (link_conut*param_count елементів)
	double* X;		 //точки наближення (link_count+1 елементів)
	__float128* A128;
	__float128* X128;
} herm_params;

typedef __float128 (*function)(const __float128);

__float128 PE4(const __float128 x, const __float128 a[4]);

__float128 dPE4(const __float128 x, const __float128 a[4]);

__float128 PE5(const __float128 x, const __float128 a[5]);

__float128 dPE5(const __float128 x, const __float128 a[5]);

__float128 PN4(const __float128 x, const __float128 a[4]);

__float128 dPN4(const __float128 x, const __float128 a[4]);

__float128 PN5(const __float128 x, const __float128 a[5]);

__float128 dPN5(const __float128 x, const __float128 a[5]);

__float128 HermiteSpline(const herm_params hp, const __float128 x, const char derivative);

int HermGen(function _f[], herm_params* hp, const __float128 a, const __float128 b, const __float128 nu);

__float128 finderr(__float128 (*link)(const __float128, const __float128[]), const __float128 params[], function f, const __float128 from, const __float128 to);
