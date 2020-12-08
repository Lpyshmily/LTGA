#include <iostream>
#include "localConst.h"
#include "ga_nlopt.h"

celestial_body Earth(-1, 54000, 0.999988049532578, 1.671681163160e-2, 0.8854353079654e-3*D2R, 175.40647696473*D2R, 287.61577546182*D2R, 257.60683707535*D2R, 0, 1);
celestial_body Jupiter(-1, 56000, 5.202646676075939E+00, 4.894528634473490E-02,	1.304070596206167E+00*D2R, 1.004989825909724E+02*D2R, 2.738038582221466E+02*D2R, 3.324137067661252E+01*D2R, 0, 1);

int main()
{
	printf("Begin!\n");
	clock_t start, stop;
	start = clock();

	GA_PSO();

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);
	return 0;
}