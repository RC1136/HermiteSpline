
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
	void* A128;
	void* X128;
} herm_params;

typedef double (*function)(const double);

double PE4(const double x, const double a[4]);

double dPE4(const double x, const double a[4]);

double PE5(const double x, const double a[5]);

double dPE5(const double x, const double a[5]);

double PN4(const double x, const double a[4]);

double dPN4(const double x, const double a[4]);

double PN5(const double x, const double a[5]);

double dPN5(const double x, const double a[5]);

double HermiteSpline(const herm_params hp, const double x, const char derivative);

int HermGen(function _f[], herm_params* hp, const double a, const double b, const double nu);

double finderr(double (*link)(const double, const double[]), const double params[], function f, const double from, const double to);
