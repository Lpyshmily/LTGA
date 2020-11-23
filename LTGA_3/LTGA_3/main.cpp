#include <iostream>
#include <ctime>
#include "localConst.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

// x[1] 引力辅助时间
// 引力辅助半径限定为最小引力辅助半径
// para NULL
double GA_obj_planar(const double* x, const double* para)
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
	double rvm[6] = {0.0}, temprv[6] = {0.0};
	double Out1[10] = {0.0}; // 第一段时间最优飞越输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out2[15] = {0.0}; // 第一段燃料最优飞越输出结果，[0-6]末端状态，[7-14]8个需要打靶的协态初值
	double Out3[10] = {0.0}; // 第二段时间最优交会输出结果，[0]剩余质量，[1-9]9个打靶变量
	double Out4[9] = {0.0}; // 第二段燃料最优交会输出结果，[0]剩余质量，[1-8]8个打靶变量
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	factor = *x;
	m0 = 20000.0/MUnit;
	const int RepeatTime = 10; // 时间最优重复求解次数
	
	printf("**********\nfractor=%f\n", factor);
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
	// 错误模型
	/*
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, tempV1[3] = {0.0}, tempV2[3] = {0.0};
	double temp, norm_vm;

	V_Copy(vm, &rvm[3], 3); // 火星速度
	norm_vm = V_Norm2(vm, 3);
	V_Minus(vin, &Out2[3], vm, 3); // 引力辅助前的相对速度
	temp = V_Dot(vin, vm, 3)/(norm_vm*norm_vm);
	V_Multi(tempV1, vm, temp, 3); // 相对速度沿vm的分量
	V_Minus(tempV2, vin, tempV1, 3); // 相对速度垂直于vm的分量
	V_Minus(vout, tempV2, tempV1, 3); // 引力辅助后的相对速度
	// 计算错误模型下的引力辅助半径
	// 结果显示引力辅助半径确实小于所允许的最小引力辅助半径，模型不合适
	double norm_vin = V_Norm2(vin, 3);
	double delta = acos(V_Dot(vin, vout, 3)/(norm_vin*norm_vin));
	double rp = muNU_MARS*(1.0/sin(delta/2) - 1.0)/(norm_vin*norm_vin);
	*/
	
	// 正确模型
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

	return -Out4[0]*MUnit;
}

// 给定一系列x的值，依次计算
void GA_planar_list()
{
	time_t now = time(NULL);
	char filename[30];
	sprintf(filename, "info_%d.txt", now);
	FILE *fid = fopen(filename, "w");

	double factor = 0.33, mf;
	while (factor<=0.35+1e-3)
	{
		printf("**********\nfractor=%f\n", factor);
		mf = GA_obj_planar(&factor, NULL);
		fprintf(fid, "%f\t%f\n", factor, mf);
		factor += 0.001;
	}

	fclose(fid);
}

// x[3] 0-引力辅助时间 1-引力辅助半径 2-vout方向角
// 引力辅助时间指第引力辅助前第一段轨迹的转移时间，变化范围为[0,tf]
// 引力辅助半径变化范围为[rmin,2rmin]
// vout方向角指vout和vin、vp组成的平面的夹角，变化范围为[0,2pi]
// para NULL
double GA_obj_triple(const double* x, const double* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0];
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_2PI;
	
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
	
	printf("**********\nfractor=%f\n", factor);
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

	return -Out4[0]*MUnit;
}

// 给定一系列x的值，依次计算
void GA_triple_list()
{
	time_t now = time(NULL);
	char filename[30];
	sprintf(filename, "info_%d.txt", now);
	FILE *fid = fopen(filename, "w");

	double x[3] = {0.0};
	double mf = 0.0;
	int i, j, k;
	for (i=0;i<=20;++i)
	{
		for (j=0;j<=10;++j)
		{
			for (k=0;k<=10;++k)
			{
				x[0] = 0.25 + i*0.01;
				x[1] = j*0.1;
				x[2] = k*0.1;
				printf("**********\n");
				printf("factor=%f\n", x[0]);
				printf("rp=%f\n", x[1]);
				printf("angle=%f\n", x[2]);
				mf = GA_obj_triple(x, NULL);
				fprintf(fid, "%f\t%f\t%f\t%f\n", x[0], x[1], x[2], mf);
			}
			fprintf(fid, "\n");
		}
	}

	fclose(fid);
}

// x[3] 0-引力辅助时间 1-引力辅助半径 2-vout方向角
// 引力辅助时间指第引力辅助前第一段轨迹的转移时间，变化范围为[0.3,0.35]tf
// 引力辅助半径变化范围为[1,2]rmin
// vout方向角指vout和vin、vp组成的平面的夹角，变化范围为[0,pi]
// para NULL
double GA_obj_PSO(const double* x, const double* para)
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

	return -Out4[0]*MUnit;
}

void GA_PSO()
{
	double xbest[3] = {0.0}, fbest;
	int D, Np;
	D = 3; // 变量个数
	Np = 20; // 粒子数量
	double wa[150]; // Np+D+3*Np*D
	PSO(GA_obj_PSO, xbest, fbest, NULL, D, Np, wa);
	printf("最优的比例系数为%f\n", xbest[0]*0.05 + 0.3);
	printf("rp=%f*rmin\n", xbest[1] + 1.0);
	printf("angle=%f*pi\n", xbest[2]);
	printf("剩余质量为%f\n", fbest);
}

int main()
{
	printf("Hello!\n");

	clock_t start, stop;
	start = clock();

	// double factor = 0.344;
	// printf("obj返回值：%f\n", GA_obj_planar(&factor, NULL));

	// GA_planar_list();

	// double x_test[3] = {0.344, 0.0, 0.25};
	// printf("obj返回值：%f\n", GA_obj_triple(x_test, NULL));

	// GA_triple_list();

	// double x_test[3] = {0.8, 0.0, 0.6}; // 对应引力辅助时间0.34，引力辅助半径rmin，旋转角0.6pi
	// printf("obj返回值：%f\n", GA_obj_PSO(x_test, NULL));

	GA_PSO();

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}