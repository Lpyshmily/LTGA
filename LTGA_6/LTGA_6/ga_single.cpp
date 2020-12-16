#include "ga_single.h"
#include "nlopt.h"
#include "localConst.h"

#pragma comment(lib,"nlopt.lib")//必须包含的静态库，通过此静态库调用动态库nlopt.dll

double ga_single1_obj(double x)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000
	
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
	double rvm[6] = {0.0}, rv_after[6] = {0.0};
	double Out1[10] = {0.0}; // 第一段时间最优飞越输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out2[15] = {0.0}; // 第一段燃料最优飞越输出结果，[0-6]末端状态，[7-14]8个需要打靶的协态初值
	double Out3[10] = {0.0}; // 第二段时间最优交会输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out4[9] = {0.0}; // 第二段燃料最优交会输出结果，[0]剩余质量，[1-8]8个打靶变量
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	m0 = 20000.0/MUnit;
	tempm = m0;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	factor = x;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	shortest1 = MaxNum;// 首先设置成一个很大的值
	shortest2 = MaxNum;// 首先设置成一个很大的值

	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段
	// 时间最优飞越
	for (j=0;j<RepeatTime;++j)
	{
		// printf("第%d次求解第一段时间最优飞越问题\n", j+1);
		flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out1[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out1[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out1[i]);
		*/
		if (Out1[9] < shortest1)
			shortest1 = Out1[9];
	}
	// printf("最短转移时间为:%.3f天\n", shortest1*TUnit/86400);
	// 判断能否完成第一段转移
	if (shortest1 > t1)
	{
		// printf("无法完成第一段轨迹转移\n");
		// fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
		return MaxNum;
	}
	// 燃料最优飞越
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out2[6]*MUnit);
	printf("末端状态量为:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("打靶变量值为:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);
	*/
	
	// 已知引力辅助半径时的引力辅助
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

	delta = 2*asin(muNU_MARS/(muNU_MARS + rminU_MARS*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] + sin(delta)*unit2[i] + 0.0*unit3[i]);
	}
	


	// 第二段 用temprv表示引力辅助后的位置速度
	V_Copy(rv_after, rvm, 3); // 引力辅助后的位置
	V_Add(&rv_after[3], vm, vout, 3); // 引力辅助后的速度
	tempm = Out2[6];// 引力辅助后的质量
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		// printf("第%d次求解第二段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, rv_after, rv1, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out3[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out3[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out3[i]);
		*/
		if (Out3[9] < shortest2)
			shortest2 = Out3[9];
	}
	// printf("最短转移时间为:%.3f天\n", shortest2*TUnit/86400);
	// 判断能否完成第二段转移
	if (shortest2 > t2)
	{
		// printf("剩余时间无法完成第二段轨迹转移\n");
		return MaxNum;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out4, rv_after, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	*/

	return -Out4[0]*MUnit;
}

void ga_single1_obj_test()
{
	double x = 0.33;
	printf("%f对应的目标函数值为%f\n", x, ga_single1_obj(x));
}

void ga_single1_list()
{
	FILE *fid = fopen("ga_single1_obj_list.txt", "w");
	double x0 = 0.0, xf = 1.0, dx = 0.01;
	int segNum = 100;
	for (int i=0;i<=segNum;++i)
	{
		fprintf(fid, "%f\t%f\n", x0 + i*dx, ga_single1_obj(x0 + i*dx));
	}
	fclose(fid);
}

double ga_single2_obj(double x)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000
	
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
	double rvm[6] = {0.0}, rv_after[6] = {0.0};
	double Out1[10] = {0.0}; // 第一段时间最优飞越输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out2[15] = {0.0}; // 第一段燃料最优飞越输出结果，[0-6]末端状态，[7-14]8个需要打靶的协态初值
	double Out3[10] = {0.0}; // 第二段时间最优交会输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out4[9] = {0.0}; // 第二段燃料最优交会输出结果，[0]剩余质量，[1-8]8个打靶变量
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	m0 = 20000.0/MUnit;
	tempm = m0;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	factor = x;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	shortest1 = MaxNum;// 首先设置成一个很大的值
	shortest2 = MaxNum;// 首先设置成一个很大的值

	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段
	// 时间最优飞越
	for (j=0;j<RepeatTime;++j)
	{
		// printf("第%d次求解第一段时间最优飞越问题\n", j+1);
		flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out1[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out1[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out1[i]);
		*/
		if (Out1[9] < shortest1)
			shortest1 = Out1[9];
	}
	// printf("最短转移时间为:%.3f天\n", shortest1*TUnit/86400);
	// 判断能否完成第一段转移
	if (shortest1 > t1)
	{
		// printf("无法完成第一段轨迹转移\n");
		// fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
		return MaxNum;
	}
	// 燃料最优飞越
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out2[6]*MUnit);
	printf("末端状态量为:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("打靶变量值为:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);
	*/
	
	// 已知引力辅助半径时的引力辅助
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

	delta = 2*asin(muNU_MARS/(muNU_MARS + rminU_MARS*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] - sin(delta)*unit2[i] + 0.0*unit3[i]);
	}
	


	// 第二段 用temprv表示引力辅助后的位置速度
	V_Copy(rv_after, rvm, 3); // 引力辅助后的位置
	V_Add(&rv_after[3], vm, vout, 3); // 引力辅助后的速度
	tempm = Out2[6];// 引力辅助后的质量
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		// printf("第%d次求解第二段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, rv_after, rv1, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out3[0]*MUnit);
		printf("转移时间为:%.3f天\n", Out3[9]*TUnit/86400);
		printf("打靶变量值为:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out3[i]);
		*/
		if (Out3[9] < shortest2)
			shortest2 = Out3[9];
	}
	// printf("最短转移时间为:%.3f天\n", shortest2*TUnit/86400);
	// 判断能否完成第二段转移
	if (shortest2 > t2)
	{
		// printf("剩余时间无法完成第二段轨迹转移\n");
		return MaxNum;
	}
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out4, rv_after, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	*/

	return -Out4[0]*MUnit;
}

void ga_single2_list()
{
	FILE *fid = fopen("ga_single2_obj_list.txt", "w");
	double x0 = 0.0, xf = 1.0, dx = 0.01;
	int segNum = 100;
	for (int i=0;i<=segNum;++i)
	{
		fprintf(fid, "%f\t%f\n", x0 + i*dx, ga_single2_obj(x0 + i*dx));
	}
	fclose(fid);
}

double ga_single1_obj_nlopt(unsigned n, const double* x, double* grad, void* para)
{
	return ga_single1_obj(*x);
}

void ga_single1_nlopt()
{
	double f_min = 0;
	double tol = 1e-5;
	double x = 0.34;
	double lb = 0.0;
	double rb = 1.0;

	nlopt_opt opter = nlopt_create(NLOPT_LN_COBYLA, 1);
	nlopt_set_min_objective(opter, ga_single1_obj_nlopt, NULL);
	nlopt_set_lower_bounds(opter, &lb);
	nlopt_set_upper_bounds(opter, &rb);
	nlopt_set_xtol_rel(opter, tol);

	nlopt_result res = nlopt_optimize(opter, &x, &f_min);
	printf("%.15f,\n", x);
	printf("最小值为%.15f\n", f_min);
}
