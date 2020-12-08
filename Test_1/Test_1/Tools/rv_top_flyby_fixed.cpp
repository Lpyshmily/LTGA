#include "rv_top_flyby_fixed.h"
#include "vector_operation.h"
#include "ODE45.h"
#include "dynamics_rv.h"
#include "ham_rv.h"
#include "constants.h"
#include "MinPack\cminpack.h"

// n 打靶变量个数 9
/* x[9] 打靶变量
[0] lam0
[1-7] 协态变量
[8] tf */
/* fvec[9] 打靶偏差 此函数的输出量
[0-2] 末端位置约束
[3-5] 末端速度协态为零
[6] 末端质量协态为零
[7] 协态归一化
[8] 静态条件*/
// iflag为0时表示上方调用函数告知不必算残差
/* sfpara[14]
[0-5] rv0
[6] m0
[7-12] rv1 目标天体的位置速度
[13] 是否输出的标志 */
// 计算成功时，返回0或大于0的值，不成功时返回小于0的值
int fvec_rv_top_flyby_fixed(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	int i;
	
	double x0[14] = {0.0};
	double rvf[6] = {0.0};

	V_Copy(x0, sfpara, 7); // 状态初值
	V_Copy(rvf, &sfpara[7], 6);
	int outflag = NINT(sfpara[13]);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7); // 协态初值
	double tf = x[8];

	if(iflag==0)
	{
		return 0;
	}

	double AbsTol[14]={0.0};
	for(i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;
	double work[140]={0.0};
	FILE *fid=NULL;//fopen("temp.txt","w");//如果设定有效文件路径，最后需要关闭文件
	int flag, NumPoint;
	flag = ode45(dynamics_rv_top, x0, NULL, 0.0, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);

	V_Minus(fvec, x0, rvf, 3);
	V_Copy(&fvec[3], &x0[10], 3); // 速度协态
	fvec[6] = x0[13];
	fvec[7] = V_Norm2(x, 8) - 1.0;
	fvec[8] = ham_rv_top(x0, lam0); // 静态条件 H(tf)=0

	if (outflag>0)
	{
		fvec[0] = x0[6]; // 输出剩余质量
	}
	// fclose(fid);
	return 0;
}
// Out[10] [0]剩余质量 [1-9]9个打靶变量
int solve_rv_top_flyby_fixed(double* Out, const double* rv0, const double* rvf, double m0, int MaxGuessNum)
{
	int n, j, info, flag=0;
	
	double sfpara[14] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rvf, 6);
	sfpara[13] = 0.0;

	n = 9; // 打靶变量个数
	double x[9] = {0.0};
	double fvec[9] = {0.0};

	Out[0] = 0.0;
	int num = 1;
	double wa[200]={0.0};
	double xtol=1.0e-8;
	double angle[7];
	while(num<=MaxGuessNum)
	{
		for(j=0;j<7;j++)
				angle[j] = (double)rand()/RAND_MAX;
		angle[0] = M_PI2*angle[0];
		angle[1] = M_PI2*angle[1];
		angle[2] = M_PI2*angle[2];
		angle[3] = M_PI*(angle[3] - 0.5);
		angle[4] = M_PI*(angle[4] - 0.5);
		angle[5] = M_2PI*angle[5];
		angle[6] = M_2PI*angle[6];

		x[0] = sin(angle[0]);
		x[1] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*cos(angle[5]);
		x[2] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*sin(angle[5]);
		x[3] = cos(angle[0])*cos(angle[1])*cos(angle[2])*sin(angle[3]);
		x[4] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*cos(angle[6]);
		x[5] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*sin(angle[6]);
		x[6] = cos(angle[0])*cos(angle[1])*sin(angle[2])*sin(angle[4]);
		x[7] = cos(angle[0])*sin(angle[1]);

		x[8] = (double)rand()/RAND_MAX*5*M_2PI;
		
		info = hybrd1(fvec_rv_top_flyby_fixed, n, x, fvec, sfpara, wa, xtol, 20, 2000);
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[13]=1.0;
			j=fvec_rv_top_flyby_fixed(n, x, fvec, 1, sfpara);
			if(fvec[0]>0.0 && x[8]>0.0) // 剩余质量和转移时间必须为正
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 9);
				break;
			}
			sfpara[13]=0.0;
		}
		num++;
	}
	printf("随机猜测次数：%d\n", num);
	return flag;
}