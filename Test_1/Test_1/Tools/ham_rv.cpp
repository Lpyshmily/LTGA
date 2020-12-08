#include "ham_rv.h"
#include "dynamics_rv.h"
#include "vector_operation.h"
#include "constants.h"

double ham_rv_top(const double* x, double lam0)
{
	double dx[14] = {0.0};
	int flag = dynamics_rv_top(0, x, dx, NULL);

	double H = V_Dot(&x[7], dx, 7);
	H += lam0;
	return H;
}

double ham_rv_fop(const double* x, double epsi, double lam0)
{
	// 计算推力大小
	double r = V_Norm2(x, 3);
	double m = x[6];
	double costate[7] = {0.0};
	V_Copy(costate, &x[7], 7);
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/(2*epsi);
	// 计算状态变量的导数
	double dx[14] = {0.0};
	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;
    int flag = dynamics_rv_fop(0, x, dx, dfpara);

    double H = V_Dot(&x[7], dx, 7);
	H += lam0*TmaxNU/Ispg0NU*(u - epsi*u*(1-u));
    return H;
}