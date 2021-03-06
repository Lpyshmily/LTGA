#include <iostream>
#include "localConst.h"
#include "ga_nlopt.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

int main()
{
	printf("Begin!\n");
	clock_t start, stop;
	start = clock();

	// test_GA_obj_nlopt();
	// test_GA_obj_PSO();

	GA_PSO();
	// GA_nlopt();

	// output_ode();
	// output_u();

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);
	return 0;
}