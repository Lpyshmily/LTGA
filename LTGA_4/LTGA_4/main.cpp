#include <iostream>
#include "ga_nlopt.h"

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

	// double x[6] = {0.999709099304766/2, 0.004722508503925, 0.470617494430747, 0.239538394402327, 0.275113612765333, 0.666440534595599};
	// double result = GA_obj_PSO_6_new(x, NULL);
	// printf("剩余质量为%.15fkg\n", result);
	GA_PSO_6_new();
	return 0;
}