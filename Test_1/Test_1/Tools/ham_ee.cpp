#include "ham_ee.h"
#include "dynamics_ee.h"
#include "vector_operation.h"

double ham_ee_top(const double* x, double lam0)
{
	double dx[14] = {0.0};
	int flag = dynamics_ee_top(0, x, dx, NULL);

	double H = V_Dot(&x[7], dx, 7);
	H += lam0;
	return H;
}