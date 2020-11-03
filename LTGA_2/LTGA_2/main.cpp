#include <iostream>
#include <ctime>
#include "localConst.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

void test_ee_fop_rend()
{
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // 飞行时间为1000天，进行归一化
	// 初始位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// 末端位置和速度,，单位分别为AU和AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// 至此，所有初始条件归一化完毕

	
	// 先用春分点轨道根数进行求解
	int flag;
	// 将位置速度转化为春分点轨道根数
	double ee0[6]={0.0}, ee1[6]={0.0};
	rv2ee(flag, ee0, rv0, muNU);
	rv2ee(flag, ee1, rv1, muNU);
	ee1[5] = ee1[5] + 3*M_2PI;
	// 求解算法的一些参数设置
	double epsi = 1.0;//取定一个较小的同伦参数直接求解近似邦邦控制的结果
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[9] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = solve_ee_fop_rend(Out, ee0, ee1, m0, tof, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0为:%.6e\n", Out[1]);
	printf("7个初始协态变量值为:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_fop_rend()
{
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // 飞行时间为1000天，进行归一化
	// 初始位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// 末端位置和速度,，单位分别为AU和AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// 至此，所有初始条件归一化完毕

	
	// 先用春分点轨道根数进行求解
	int flag;
	// 求解算法的一些参数设置
	double epsi = 1.0;//取定一个较小的同伦参数直接求解近似邦邦控制的结果
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[9] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = solve_rv_fop_rend(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0为:%.6e\n", Out[1]);
	printf("7个初始协态变量值为:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_top_rend_fixed()
{
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	// 初始位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// 末端位置和速度,，单位分别为AU和AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// 至此，所有初始条件归一化完毕

	
	// 先用春分点轨道根数进行求解
	int flag;
	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[10] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = solve_rv_top_rend_fixed(Out, rv0, rv1, m0, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("转移时间为:%.3f天\n", Out[9]*TUnit/86400);
	printf("7个初始协态变量值为:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_fop_flyby()
{
	int i;
	
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // 飞行时间为1000天，进行归一化
	// 初始位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// 末端位置和速度,，单位分别为AU和AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// 至此，所有初始条件归一化完毕

	
	// 先用春分点轨道根数进行求解
	int flag;
	// 求解算法的一些参数设置
	double epsi = 1.0;//取定一个较小的同伦参数直接求解近似邦邦控制的结果
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[15] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = solve_rv_fop_flyby(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[6]*MUnit);
	printf("末端状态量为:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out[i]);
	printf("打靶变量值为:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out[i]);
}

void test_GA()
{
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=6000 Tmax=2.26 m0=20000
	
	int i, flag;
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0};
	double Out1[15] = {0.0}; // 输出计算结果，[0-6]-末端状态，[7-14]-8个需要打靶的协态初值
	double Out2[9] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	double Out3[10] = {0.0};
	double m0, tf, epsi, factor, t1, t2;
	tf = 2201*86400/TUnit;
	epsi = 1.0;
	m0 = 20000.0/MUnit;

	// 首先随机给定一个参数，代表第一段轨迹转移时间在总时间中的占比
	factor = 0.37;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段 燃料最优飞越
	flag = solve_rv_fop_flyby(Out1, rv0, rvm, m0, t1, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out1[6]*MUnit);
	printf("末端状态量为:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out1[i]);
	printf("打靶变量值为:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out1[i]);
	

	// 引力辅助
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, tempV1[3] = {0.0}, tempV2[3] = {0.0};
	double temp, norm_vm;

	V_Copy(vm, &rvm[3], 3); // 火星速度
	norm_vm = V_Norm2(vm, 3);
	V_Minus(vin, &Out1[3], vm, 3); // 引力辅助前的相对速度
	temp = V_Dot(vin, vm, 3)/(norm_vm*norm_vm);
	V_Multi(tempV1, vm, temp, 3); // 相对速度沿vm的分量
	V_Minus(tempV2, vin, tempV1, 3); // 相对速度垂直于vm的分量
	V_Minus(vout, tempV2, tempV1, 3); // 引力辅助后的相对速度

	// 第二段 用rv0表示引力辅助后的位置速度
	V_Copy(rv0, rvm, 3); // 引力辅助后的位置
	V_Add(&rv0[3], vm, vout, 3); // 引力辅助后的速度
	m0 = Out1[6];// 引力辅助后的质量
	// 时间最优交会
	
	flag = solve_rv_top_rend_fixed(Out3, rv0, rv1, m0, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out3[0]*MUnit);
	printf("转移时间为:%.3f天\n", Out3[9]*TUnit/86400);
	printf("打靶变量值为:\n");
	for (i=1; i<10; i++)
		printf("%.15e,\n", Out3[i]);

	// 判断能否完成第二段转移
	if (Out3[9] > t2)
	{
		printf("剩余时间无法完成第二段轨迹转移\n");
		return;
	}

	// 燃料最优交会
	flag = solve_rv_fop_rend(Out2, rv0, rv1, m0, t2, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out2[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out2[i]);
}

void test_GA_factor()
{
	// 需要对Tools.h/constants.h中的参数进行修改
	// Isp=6000 Tmax=2.26 m0=20000

	FILE *fid = fopen("info_5.txt", "w");
	
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
	m0 = 20000.0/MUnit;
	const int RepeatTime = 10; // 时间最优重复求解次数

	// 首先随机给定一个参数，代表第一段轨迹转移时间在总时间中的占比
	factor = 0.34;
	while (factor<=0.36)
	{
		printf("**********\nfractor=%f\n", factor);
		tempm = m0;
		shortest1 = 1e7;// 首先设置成一个很大的值
		shortest2 = 1e7;// 首先设置成一个很大的值
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
			printf("求解成功%d\n",flag);
			printf("剩余质量为:%.3fkg\n", Out1[0]*MUnit);
			printf("转移时间为:%.3f天\n", Out1[9]*TUnit/86400);
			printf("打靶变量值为:\n");
			for (i=1; i<10; i++)
				printf("%.15e,\n", Out1[i]);
			if (flag && (Out1[9] < shortest1))
				shortest1 = Out1[9];
		}
		printf("最短转移时间为:%.3f天\n", shortest1*TUnit/86400);
		// 判断能否完成第一段转移
		if (shortest1 > t1)
		{
			printf("无法完成第一段轨迹转移\n");
			fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
			factor += 0.001;
			continue;
		}
		// 燃料最优飞越
		flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out2[6]*MUnit);
		printf("末端状态量为:\n");
		for (i=0;i<7;++i)
			printf("%.15e,\n", Out2[i]);
		printf("打靶变量值为:\n");
		for (i=7; i<15; i++)
			printf("%.15e,\n", Out2[i]);
	

		// 引力辅助
		double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, tempV1[3] = {0.0}, tempV2[3] = {0.0};
		double temp, norm_vm;

		V_Copy(vm, &rvm[3], 3); // 火星速度
		norm_vm = V_Norm2(vm, 3);
		V_Minus(vin, &Out2[3], vm, 3); // 引力辅助前的相对速度
		temp = V_Dot(vin, vm, 3)/(norm_vm*norm_vm);
		V_Multi(tempV1, vm, temp, 3); // 相对速度沿vm的分量
		V_Minus(tempV2, vin, tempV1, 3); // 相对速度垂直于vm的分量
		V_Minus(vout, tempV2, tempV1, 3); // 引力辅助后的相对速度

		// 第二段 用temprv表示引力辅助后的位置速度
		V_Copy(temprv, rvm, 3); // 引力辅助后的位置
		V_Add(&temprv[3], vm, vout, 3); // 引力辅助后的速度
		tempm = Out2[6];// 引力辅助后的质量
		// 时间最优交会
		for (j=0;j<RepeatTime;++j)
		{
			printf("第%d次求解第二段时间最优交会问题\n", j+1);
			flag = solve_rv_top_rend_fixed(Out3, temprv, rv1, tempm, MaxGuessNum);
			printf("求解成功%d\n",flag);
			printf("剩余质量为:%.3fkg\n", Out3[0]*MUnit);
			printf("转移时间为:%.3f天\n", Out3[9]*TUnit/86400);
			printf("打靶变量值为:\n");
			for (i=1; i<10; i++)
				printf("%.15e,\n", Out3[i]);
			if (flag && (Out3[9] < shortest2))
				shortest2 = Out3[9];
		}
		printf("最短转移时间为:%.3f天\n", shortest2*TUnit/86400);
		// 判断能否完成第二段转移
		if (shortest2 > t2)
		{
			printf("剩余时间无法完成第二段轨迹转移\n");
			fprintf(fid, "%f\t%f\t%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400, t2*TUnit/86400, shortest2*TUnit/86400);
			factor += 0.001;
			continue;
		}
		// 燃料最优交会
		flag = solve_rv_fop_rend(Out4, temprv, rv1, tempm, t2, epsi, MaxGuessNum);
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
		printf("打靶变量值为:\n");
		for (i=1; i<9; i++)
			printf("%.15e,\n", Out4[i]);

		fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400, t2*TUnit/86400, shortest2*TUnit/86400, Out4[0]*MUnit);

		factor += 0.001;
	}
	fclose(fid);
}

int main()
{
	printf("Hello!\n");

	clock_t start, stop;
	start = clock();

	/*
	test_ee_fop_rend();
	test_rv_fop_rend();
	test_rv_fop_flyby();
	test_rv_top_rend_fixed();
	*/

	// test_GA();
	test_GA_factor();

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}