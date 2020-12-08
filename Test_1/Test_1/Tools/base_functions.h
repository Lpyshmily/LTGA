#ifndef _BASE_FUNCTION_
#define _BASE_FUNCTION_

#include "constants.h"

double e2f(int& flag, double E, double e);
// E2M 根据偏近点角和偏心率求平近点角
double e2m(int& flag, double E, double e);
// f2E 根据真近点角和偏心率求偏近点角
double f2e(int& flag, double f, double e);
// M2E 根据平近点角和偏心率求偏近点角
double m2e(int& flag, double M, double e, int MaxIter=100, double epsilon=1.0e-14);

double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter=100, double epsilon=1.0e-14);
	
void coe2rv(int& flag, double* rv, const double* coe, const double mu);

void rv2coe(int& flag, double* coe, const double* RV, const double mu);

void coe2ee(int&flag, double* ee, const double* coe, const double mu);

void ee2coe(int&flag, double* coe, const double* ee, const double mu);

void ee2rv(int&flag, double* rv, const double* ee, const double mu);

void rv2ee(int&flag, double* ee, const double* RV, const double mu);

//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int&flag, double* coe1, const double* coe0, double dt, double mu);
//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int&flag, double* rv1, const double* rv0, double dt, double mu);
//根据初始时刻状态ee0求末端时刻dt的状态ee1，按二体推进。若计算成功,flag返回1
void ee02eef(int&flag, double* ee1, const double* ee0, double dt, double mu);

/**************************************************************************************************************************************/
/*********************************************************考虑j2摄动的轨道根数转换*****************************************************/
/**************************************************************************************************************************************/
void coe2soe(const double* coe, double* soe);
void soe2coe(const double* soe, double* coe);
//输入为平均经典轨道根数，输出为瞬时轨道根数
void M2O(const double* soe_m, double* soe_o, double j2=J2[Ncenter]);
void O2M(const double* soe_o, double* soe_m, double j2=J2[Ncenter]);
int myfun(int n, const double* x, double* fvec, int iflag, const double* para);
/**************************************************************************************************************************************/
/************************************************************考虑j2摄动的轨道递推******************************************************/
/**************************************************************************************************************************************/
void j2mcoe02mcoef(const double* me0, const double dt, double* mef);
void j2ocoe02ocoef(double* ocoef, const double* ocoe0, const double dt);
void j2rv02rvf(const double* rv0, const double dt, double* rvf);
	
inline int dyn_j2(double t, const double* x, double* dx, const double* hpara)
{
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	double r3 = r*r*r;
	double frz = 1+1.5*J2[Ncenter]*Ra[Ncenter]*Ra[Ncenter]/r/r*(1-5*x[2]*x[2]/r/r);
	double fr = 3*J2[Ncenter]*Ra[Ncenter]*Ra[Ncenter]/r/r;

	dx[0] = x[3];
	dx[1] = x[4];
	dx[2] = x[5];
	dx[3] = -GMMu[Ncenter]*x[0]/r3*frz;
	dx[4] = -GMMu[Ncenter]*x[1]/r3*frz;
	dx[5] = -GMMu[Ncenter]*x[2]/r3*(frz+fr);

	return 1;
}

inline int dyn(double t, const double* x, double* dx, const double* hpara)
{
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	double r3 = r*r*r;

	dx[0] = x[3];
	dx[1] = x[4];
	dx[2] = x[5];
	dx[3] = -GMMu[Ncenter]*x[0]/r3;
	dx[4] = -GMMu[Ncenter]*x[1]/r3;
	dx[5] = -GMMu[Ncenter]*x[2]/r3;

	return 1;
}


#endif