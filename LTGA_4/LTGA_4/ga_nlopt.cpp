#include "ga_nlopt.h"
#include "localConst.h"
#include "nlopt.h"

#pragma comment(lib,"nlopt.lib")//必须包含的静态库，通过此静态库调用动态库nlopt.dll

int global_count = 0;

// x[3] 0-引力辅助时间 1-引力辅助半径 2-vout方向角
// 三个分量的变化范围均为0-1，对应以下不同含义
// 引力辅助时间指第引力辅助前第一段轨迹的转移时间，变化范围为[0.3,0.35]tf
// 引力辅助半径变化范围为[1,2]rmin
// vout方向角指vout和vin、vp组成的平面的夹角，变化范围为[0,pi]
// para NULL
double GA_obj_nlopt(unsigned n, const double* x, double* grad, void* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.05 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	
	int i, j, flag;
	printf("**********\n");
	for (i=0;i<3;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");
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
			return MaxNum;
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
		return MaxNum;
	}
	// 燃料最优飞越
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
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
			return MaxNum;
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
		return MaxNum;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out4, temprv, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	printf("**********\n");
	printf("迭代次数 i=%d\n", ++global_count);
	for (i=0;i<3;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");

	return -Out4[0]*MUnit;
}

void test_GA_obj_nlopt()
{
	double x_test[3] = {0.876959627917576,0.0,0.471123454181230};
	printf("obj返回值：%f\n", GA_obj_nlopt(3, x_test, NULL, NULL));
}

void GA_nlopt()
{
	double f_min = 0;
	double tol = 1e-5;
	// double x[3] = {0.87704,0.0,0.469045};
	double x[3] = {0.5, 0.5, 0.5};
	double lb[3] = {0.0, 0.0, 0.0};
	double rb[3] = {1.0, 1.0, 1.0};

	nlopt_opt opter = nlopt_create(NLOPT_LN_COBYLA, 3);
	nlopt_set_min_objective(opter, GA_obj_nlopt, NULL);
	nlopt_set_lower_bounds(opter, lb);
	nlopt_set_upper_bounds(opter, rb);
	nlopt_set_xtol_rel(opter, tol);

	nlopt_result res = nlopt_optimize(opter, x, &f_min);
	for (int i=0;i<3;++i)
		printf("%.15f,\n", x[i]);
	printf("最小值为%.15f\n", f_min);
}

// x[6] 0-引力辅助时间 1-引力辅助半径 2-vout方向角 3-5 引力辅助前相对速度vin的3个分量
// 求解两段燃料最优交会问题
// 三个分量的变化范围均为0-1，对应以下不同含义
// 引力辅助时间指第引力辅助前第一段轨迹的转移时间，变化范围为[0.3,0.35]tf
// 引力辅助半径变化范围为[1,2]rmin
// vout方向角指vout和vin、vp组成的平面的夹角，变化范围为[0,pi]
// vin3个分量的变化范围均为[-0.1,0.1]
// para NULL
double GA_obj_nlopt_6(unsigned n, const double* x, double* grad, void* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.05 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	double vin[3] = {x[3]*0.2-0.1, x[4]*0.2-0.1, x[5]*0.2-0.1}; // 引力辅助前的相对速度
	
	int i, j, flag;
	printf("**********\n");
	for (i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0}, rv_before[6] = {0.0}, rv_after[6] = {0.0}; // 火星位置速度、引力辅助前的位置速度、引力辅助后的位置速度
	double Out1[10] = {0.0}; // 第一段时间最优交会输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out2[9] = {0.0}; // 第一段燃料最优交会输出结果，[0]剩余质量，[1-8]8个需要打靶的协态初值
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

	V_Copy(rv_before, rvm, 3);
	for (i=0;i<3;++i)
		rv_before[i+3] = rvm[i+3] + vin[i];


	// 求解算法的一些参数设置
	int MaxGuessNum = 1000;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		printf("第%d次求解第一段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out1, rv0, rv_before, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
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
		return MaxNum;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out2, rv0, rv_before, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
	{
		printf("第一段燃料最优交会问题不收敛\n");
		return MaxNum;
	}
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out2[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out2[i]);
	

	// 引力辅助
	double vout[3] = {0.0}, vm[3] = {0.0}, unit1[3] = {0.0}, unit2[3] = {0.0}, unit3[3] = {0.0}, tempVec[3] = {0.0};
	double delta, norm_vin, temp;
	V_Copy(vm, &rvm[3], 3); // 火星速度
	// 引力辅助前的相对速度vin作为优化变量，在函数一开始定义

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
	V_Copy(rv_after, rvm, 3); // 引力辅助后的位置
	V_Add(&rv_after[3], vm, vout, 3); // 引力辅助后的速度
	tempm = Out2[0];// 引力辅助后的质量
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		printf("第%d次求解第二段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, rv_after, rv1, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
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
		return MaxNum;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out4, rv_after, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
	{
		printf("第二段燃料最优交会问题不收敛\n");
		return MaxNum;
	}
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);

	printf("**********\n");
	printf("迭代次数 i=%d\n", ++global_count);
	for (i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");

	return -Out4[0]*MUnit;
}

void test_GA_obj_nlopt_6()
{
	// double x_test[6] = {0.876959627917576,0.0,0.471123454181230, 0.201, 0.3573, 0.6245};
	double x_test[6] = {0.8754359,0.0,0.4370229,0.2028572463,0.3,0.623};
	printf("obj返回值：%f\n", GA_obj_nlopt_6(6, x_test, NULL, NULL));
}

void GA_nlopt_6()
{
	double f_min = 0;
	double tol = 1e-5;
	double x[6] = {0.876959627917576,0.0,0.471123454181230, 0.201, 0.3573, 0.6245};
	double lb[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double rb[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	nlopt_opt opter = nlopt_create(NLOPT_LN_COBYLA, 6);
	nlopt_set_min_objective(opter, GA_obj_nlopt_6, NULL);
	nlopt_set_lower_bounds(opter, lb);
	nlopt_set_upper_bounds(opter, rb);
	nlopt_set_xtol_rel(opter, tol);

	nlopt_result res = nlopt_optimize(opter, x, &f_min);
	for (int i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("最小值为%.15f\n", f_min);
}


double GA_obj_PSO_6(const double* x, const double* para)
{
	return GA_obj_nlopt_6(6, x, NULL, NULL);
}

void GA_PSO_6()
{
	double xbest[6] = {0.0};
	double fbest;
	int D, Np;
	D = 6;
	Np = 20;
	double wa[400];
	PSO(GA_obj_PSO_6, xbest, fbest, NULL, D, Np, wa);
	for (int i=0;i<D;++i)
		printf("xbest[%d]=%.15f,\n", i, xbest[i]);
	printf("fbest=%.15f\n", fbest);
}

double GA_obj_PSO_6_new(const double* x, const double* para)
{
	double x_nlopt[6] = {0.0};
	V_Copy(x_nlopt, x, 6);
	x_nlopt[0] = x[0]*2;
	
	return GA_obj_nlopt_6(6, x_nlopt, NULL, NULL);
}

void GA_PSO_6_new()
{
	double xbest[6] = {0.0};
	double fbest;
	int D, Np;
	D = 6;
	Np = 20;
	double wa[400];
	PSO(GA_obj_PSO_6_new, xbest, fbest, NULL, D, Np, wa);
	for (int i=0;i<D;++i)
		printf("xbest[%d]=%.15f,\n", i, xbest[i]);
	printf("fbest=%.15f\n", fbest);
}
