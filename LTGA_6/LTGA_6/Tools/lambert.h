#ifndef _LAMBERT_
#define _LAMBERT_

#include "vector_operation.h"
#include "constants.h"
#include "base_functions.h"
#include <math.h>

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2, 
				double tf, const double* unith, int& flag, double mu, int way=0, int N=0, int branch=0, int Maxiter=60, double tol=1.0e-12);

void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
				const double* rv0, const double* rv1, double t, double GM);

class lambert_solver
{
public:
	lambert_solver(const double* rv1,const double* rv2,const double tof,const double mu,const int way=0);
	void solve_single(const int N=0,const int branch=0); // 求解给定圈数和分支时的速度增量
	int solve_multi();
	void get_dv1(double* dv1);
	void get_dv2(double* dv2);
	double get_Mdv1();
	double get_Mdv2();
	void get_originalRVT(double* rv);

	double m_rv1[6],m_rv2[6],m_tof,m_mu;
	int m_way;
	int m_flag;
	double m_dv1[3],m_dv2[3],m_Mdv1,m_Mdv2; // 两个速度增量及其大小
	int m_N,m_branch;

	double m_original_rvt[6];
};
#endif