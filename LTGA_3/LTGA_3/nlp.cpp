#include "nlp.h"
#include "localConst.h"
#include<iostream>
#include <stdlib.h>
#include"cnpsol.h"
#include"Use.h"

#pragma comment(lib,"cnpsol.lib")//必须包含的静态库，通过此静态库调用动态库cnpsol.dll

using namespace std;

extern celestial_body Mars;

//下面实现了两个例子，分别为matlab的fmincon中的例子和NPSOL中的例子，执行哪个例子可将另一个例子相关的代码注释掉

//用户根据实际问题编写的指标函数
bool fnobj(int n, const double* x, double& objf)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.05 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	printf("**********\n");
	printf("factor=%f\n", factor);
	printf("rp=%f*rmin\n", x[1] + 1.0);
	printf("angle=%f*pi\n", x[2]);

	int i, j, flag;
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0}, temprv[6] = {0.0};
	double Out1[10] = {0.0}; // 第一段时间最优飞越输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out2[15] = {0.0}; // 第一段燃料最优飞越输出结果，[0-6]末端状态，[7-14]8个需要打靶的协态初值
	double Out3[10] = {0.0}; // 第二段时间最优交会输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out4[9] = {0.0}; // 第二段燃料最优交会输出结果，[0]剩余质量，[1-8]8个打靶变量
	double m0, tempm, tf, epsi, t1, t2, shortest1, shortest2;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	m0 = 20000.0/MUnit;
	const int RepeatTime = 10; // 时间最优重复求解次数

	tempm = m0;
	shortest1 = MaxNum;// 首先设置成一个很大的值
	shortest2 = MaxNum;// 首先设置成一个很大的值
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);

	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段
	// 时间最优飞越
	for (j=0;j<RepeatTime;++j)
	{
		printf("第%d次求解第一段时间最优飞越问题\n", j+1);
		flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
		if (!flag)
		{
			objf = MaxNum;
			return true;
		}
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out1[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out1[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out1[i]);
		if (Out1[9] < shortest1)
			shortest1 = Out1[9];
	}
	printf("最短转移时间为:%.3f天\n", shortest1*TUnit/86400);
	// 判断能否完成第一段转移
	if (shortest1 > t1)
	{
		printf("无法完成第一段轨迹转移\n");
		// fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
		objf = MaxNum;
		return true;
	}
	// 燃料最优飞越
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
	{
		objf = MaxNum;
		return true;
	}
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out2[6]*MUnit);
	printf("末端状态量为:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("打靶变量值为:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);


	// 引力辅助
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, unit1[3] = {0.0}, unit2[3] = {0.0}, unit3[3] = {0.0}, tempVec[3] = {0.0};
	double delta, norm_vin, temp;
	V_Copy(vm, &rvm[3], 3); // 火星速度
	V_Minus(vin, &Out2[3], vm, 3); // 引力辅助前的相对速度

	norm_vin = V_Norm2(vin, 3);
	V_Divid(unit1, vin, norm_vin, 3); // 沿vin方向的单位矢量
	V_Cross(tempVec, vin, vm);
	temp = V_Norm2(tempVec, 3);
	V_Divid(unit3, tempVec, temp, 3); // 垂直vin和vm平面的单位向量
	V_Cross(unit2, unit3, unit1); // j单位矢量

	delta = 2*asin(muNU_MARS/(muNU_MARS + rp*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] + sin(delta)*sin(phi)*unit2[i] + sin(delta)*cos(phi)*unit3[i]);
	}


	// 第二段 用temprv表示引力辅助后的位置速度
	V_Copy(temprv, rvm, 3); // 引力辅助后的位置
	V_Add(&temprv[3], vm, vout, 3); // 引力辅助后的速度
	tempm = Out2[6];// 引力辅助后的质量
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		printf("第%d次求解第二段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, temprv, rv1, tempm, MaxGuessNum);
		if (!flag)
		{
			objf = MaxNum;
			return true;
		}
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out3[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out3[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out3[i]);
		if (Out3[9] < shortest2)
			shortest2 = Out3[9];
	}
	printf("最短转移时间为:%.3f天\n", shortest2*TUnit/86400);
	// 判断能否完成第二段转移
	if (shortest2 > t2)
	{
		printf("剩余时间无法完成第二段轨迹转移\n");
		objf = MaxNum;
		return true;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out4, temprv, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
	{
		objf = MaxNum;
		return true;
	}
	printf("求解成功%d\n",flag);
	printf("**********\n");
	printf("x[0]=%f\n", x[0]);
	printf("x[1]=%f\n", x[1]);
	printf("x[2]=%f\n", x[2]);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);

	objf = -Out4[0]*MUnit;
	return true;
}

//用户根据实际问题编写的非线性约束函数
bool fncon(int n, int ncnln, const double* x, double* nlc)
{
	return true;
}

void test_nlp()
{
	int i, j;
	int n=3, nclin=0, ncnln=0;
	int nctot=n+nclin+ncnln, inform=0;//nctot为包括自变量、线性约束、非线性约束在内总约束的个数，也是bu、bl的维数
	//开辟一个2维数组存储线性约束系数矩阵
	double** AM=NULL;//new double*[nclin];
	/*for(i=0;i<nclin;i++)
		AM[i]=new double[n];*/
	////开辟1维数组存储约束上、下界
	double* bl=new double[nctot];
	double* bu=new double[nctot];
	double objf=0.0;

	double bigbnd=1.0E21;//如果无约束，当作约束处理，但取个较大的值，例如对于x1>0，取相应的bl=0,bu=bigbnd，成为0<=x1<=bigbnd
	//for(i=0;i<nclin;i++)
	//	for(j=0;j<n;j++)
	//		AM[i][j]=0.0;
	
	//先全部置为正、负无穷，相当于无约束，再对约束逐个设置
	for(j=0;j<nctot;j++)
	{
		bl[j]=-bigbnd;
		bu[j]=bigbnd;
	}
	//matlab的fmincon中的例子
	bl[0]=0.0;
	bu[0]=1.0;
	bl[1]=0.0;
	bu[1]=1.0;
	bl[2] = 0.0;
	bu[2] = 1.0;

	//初值
	double x[3]={0.87704,0.000017,0.469045};
	inform=NPSol(x, objf, n, nclin, ncnln, AM, bl, bu);
	//输出最优结果和指标及计算成功与否姿态标识
	for(i=0;i<n;i++) printf("x[%d]=%.12e,\n",i,x[i]);
	printf("obj=%.12e,\n",objf);
	printf("inform=%d,\n",inform);	
	
	//释放动态内在，必须有
	delete[] bu;
	delete[] bl;
	/*for(i=0;i<nclin;i++)
		delete[] AM[i];
	delete[] AM;*/
	
	printf("Hello Down!\n");
}