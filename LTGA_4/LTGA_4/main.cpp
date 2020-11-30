#include <iostream>
#include "localConst.h"
#include "ga_nlopt.h"
#include "ga_nlopt_ee.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

int main()
{
	printf("Hello!\n");
	// test_GA_obj_nlopt();
	// GA_nlopt();

	// test_GA_obj_nlopt_6();
	// GA_nlopt_6();

	// double x[6] = {0.999709099304766, 0.004722508503925, 0.470617494430747, 0.239538394402327, 0.275113612765333, 0.666440534595599};
	// double result = GA_obj_PSO_6(x, NULL);
	// printf("剩余质量为%.15fkg\n", result);
	// GA_PSO_6();

	double x[6] = {0.999709099304766/2, 0.004722508503925, 0.470617494430747, 0.239538394402327, 0.275113612765333, 0.666440534595599};
	double result = GA_obj_PSO_6_new(x, NULL);
	printf("剩余质量为%.15fkg\n", result);
	// GA_PSO_6_new();

	// double x[6] = {0.999709099304766/2, 0.004722508503925, 0.470617494430747, 0.239538394402327, 0.275113612765333, 0.666440534595599};
	// double result = GA_obj_nlopt_6_ee(6, x, NULL, NULL);
	// printf("剩余质量为%.15fkg\n", result);
	return 0;
}