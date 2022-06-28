#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h>




int main(void)
{
	char buf[128];
	__float128 shit1 = 542e-210;
	__float128 shit2 = 312e-209;
	
	quadmath_snprintf(buf, sizeof(buf), "%Qe", powq(9999, 999));
	printf(buf);
	
	return 0;
}