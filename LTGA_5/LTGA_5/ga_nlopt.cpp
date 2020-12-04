#include "ga_nlopt.h"
#include "nlopt.h"
#include "localConst.h"

#pragma comment(lib,"nlopt.lib")//必须包含的静态库，通过此静态库调用动态库nlopt.dll

int global_count = 0;

// x[6] 0-引力辅助时间 1-引力辅助半径 2-vout方向角 3-5 引力辅助前相对速度vin的3个分量
// 求解两段燃料最优交会问题
// 三个分量的变化范围均为0-1，对应以下不同含义
// 引力辅助时间指第引力辅助前第一段轨迹的转移时间，变化范围为[0.3,0.4]tf
// 引力辅助半径变化范围为[1,2]rmin
// vout方向角指vout和vin、vp组成的平面的夹角，变化范围为[0,pi]
// vin3个分量的变化范围均为[-0.2,0.2]
// para NULL
double GA_obj_nlopt(unsigned n, const double* x, double* grad, void* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.1 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	double vin[3] = {x[3]*0.4-0.2, x[4]*0.4-0.2, x[5]*0.4-0.2}; // 引力辅助前的相对速度
	
	int i, j, flag;
	
	
	// nlopt时输出，PSO时不输出
	printf("**********\n");
	printf("迭代次数 i=%d\n", ++global_count);
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
	int MaxGuessNum = 500;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数

	// 求解
	// 第一段
	// 时间最优交会
	for (j=0;j<RepeatTime;++j)
	{
		// printf("第%d次求解第一段时间最优交会问题\n", j+1);
		flag = solve_rv_top_rend_fixed(Out1, rv0, rv_before, tempm, MaxGuessNum);
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
	// 燃料最优交会
	flag = solve_rv_fop_rend(Out2, rv0, rv_before, tempm, t1, epsi, MaxGuessNum);
	// flag = homotopy_rv_fop_rend(Out2, rv0, rv_before, tempm, t1, MaxGuessNum);
	if (!flag)
	{
		// printf("第一段燃料最优交会问题不收敛\n");
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
	// 测试已知最优解时需要
	/*
	vout[0] = 0.11925935378439018; vout[1] = -0.016224539415498812; vout[2] = 0.0028024583987117900;
	double temp1 = V_Dot(vout, unit1, 3)/norm_vin/cos(delta);
	double temp2 = V_Dot(vout, unit2, 3)/norm_vin/sin(delta);
	double temp3 = V_Dot(vout, unit3, 3)/norm_vin/sin(delta);
	*/
	


	// 第二段 用temprv表示引力辅助后的位置速度
	V_Copy(rv_after, rvm, 3); // 引力辅助后的位置
	V_Add(&rv_after[3], vm, vout, 3); // 引力辅助后的速度
	tempm = Out2[0];// 引力辅助后的质量
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
	{
		// printf("第二段燃料最优交会问题不收敛\n");
		return MaxNum;
	}
	
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	

	// nlopt时输出，PSO时不输出
	printf("剩余质量为:%.3fkg\n", Out4[0]*MUnit);

	// output需要的数据
	printf("rv0:\n");
	for (i=0;i<6;++i)
		printf("%.15f\n", rv0[i]);
	printf("rv_before:\n");
	for (i=0;i<6;++i)
		printf("%.15f\n", rv_before[i]);
	printf("rv_after:\n");
	for (i=0;i<6;++i)
		printf("%.15f\n", rv_after[i]);
	printf("rv1:\n");
	for (i=0;i<6;++i)
		printf("%.15f\n", rv1[i]);
	printf("t1:%.15f\n", t1);
	printf("t2:%.15f\n", t2);
	printf("m1:%.15f\n", m0);
	printf("m2:%.15f\n", tempm);

	return -Out4[0]*MUnit;
}

double GA_obj_PSO(const double* x, const double* para)
{
	return GA_obj_nlopt(6, x, NULL, NULL);
}

void test_GA_obj_nlopt()
{
	double x[6] = {0.863655768384999,0.000017319770127,0.499034175810391,0.621295996122853,0.228521973614765,0.515503534531755}; // 16021.090445199869
	printf("GA_obj_nlopt返回值：%f\n", GA_obj_nlopt(6, x, NULL, NULL));
}

void test_GA_obj_PSO()
{
	double x[6] = {0.863655768384999,0.000017319770127,0.499034175810391,0.621295996122853,0.228521973614765,0.515503534531755}; // 16021.090445199869
	printf("GA_obj_PSO返回值：%f\n", GA_obj_PSO(x, NULL));
}

void GA_PSO()
{
	double xbest[6] = {0.0};
	double fbest;
	int D, Np;
	D = 6;
	Np = 40;
	PSO(GA_obj_PSO, xbest, fbest, NULL, D, Np, 10, 1);
	for (int i=0;i<D;++i)
		printf("xbest[%d]=%.15f,\n", i, xbest[i]);
	printf("fbest=%.15f\n", fbest);
}

void GA_nlopt()
{
	double f_min = 0;
	double tol = 1e-5;
	double x[6] = {0.863655768384999,0.000017319770127,0.499034175810391,0.621295996122853,0.228521973614765,0.515503534531755};
	double lb[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double rb[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	nlopt_opt opter = nlopt_create(NLOPT_LN_COBYLA, 6);
	nlopt_set_min_objective(opter, GA_obj_nlopt, NULL);
	nlopt_set_lower_bounds(opter, lb);
	nlopt_set_upper_bounds(opter, rb);
	nlopt_set_xtol_rel(opter, tol);

	nlopt_result res = nlopt_optimize(opter, x, &f_min);
	for (int i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("最小值为%.15f\n", f_min);
}

void output_ode()
{
	// 需要提前对rv_fop_rend.cpp文件名称进行改动
	// 需要确定打靶变量
	printf("输出积分过程信息\n");
	double rv0[6] = {0.587642000000000,0.795462700000000,-0.000038452030000,-0.820580823694579,0.590131095461690,-0.000050802341948};
	double rv_before[6] = {0.785122098448566,-1.157915016177971,-0.043524651865926,0.758861072171299,0.419239168514373,-0.001528616508470};
	double rv_after[6] = {0.785122098448566,-1.157915016177971,-0.043524651865926,0.823166898923659,0.509786445169425,-0.003389557797740};
	double rv1[6] = {-5.204974000000000,1.495369000000000,0.110244400000000,-0.126321626888245,-0.401565532241840,0.004493379047356};
	double t1 = 14.682743622728498;
	double t2 = 23.179076166223364;
	double m1 = 1.000000000000000;
	double m2 = 0.913380468124813;
	double epsi = 1.0e-5;
	int flag, MaxGuessNum = 500;
	double Out[9] = {0.0}; // 燃料最优交会输出结果，[0]剩余质量，[1-8]8个需要打靶的协态初值

	// 第一段
	// flag = solve_rv_fop_rend(Out, rv0, rv_before, m1, t1, epsi, MaxGuessNum);
	// 第二段
	flag = solve_rv_fop_rend(Out, rv_after, rv1, m2, t2, epsi, MaxGuessNum);
}

void output_u()
{
	double x1[14] = {0.587642000000000,0.795462700000000,-0.000038452030000,-0.820580823694579,0.590131095461690,-0.000050802341948,1.000000000000000,-2.443289425469893e-001,-3.388805570248314e-001,-3.685229465217167e-002,3.223612567425680e-001,-2.402788863810115e-001,-9.699340788212746e-002,7.501358724377988e-002};
	double x2[14] = {0.785122098448566,-1.157915016177971,-0.043524651865926,0.823166898923659,0.509786445169425,-0.003389557797740,0.913380468124813,-1.273152789553519e-001,2.120126113944062e-001,-6.022980525892011e-003,-3.820774629894042e-001,-2.637371289846143e-001,1.664044729452408e-003,1.161888933616963e-001};
	double x1_original[14], x2_original[14];
	V_Copy(x1_original, x1, 14);
	V_Copy(x2_original, x2, 14);
	double epsi = 1.0e-5, lam01 = 8.046239878045387e-001, lam02 = 8.424738020525979e-001;
	double t1 = 14.682743622728498;
	double t2 = 23.179076166223364;
	int flag, NumPoint;
	double work[140]={0.0};

	// 第一段
	/*
	time_t now = time(NULL);
	char filename[30];
	sprintf(filename, "u_%d.txt", now);
	FILE *ufid = fopen(filename, "w");
	FILE *fid = NULL;

	int TimeNodeNum = 1000;
	double* TimeNode = new double[TimeNodeNum];
	double seg = t1/TimeNodeNum;
	for (int i=0;i<TimeNodeNum;++i)
		TimeNode[i] = seg*i;
	TimeNode[TimeNodeNum-1] = t1;

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam01;

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	double norm_lambdaV, rou, u;

	for (int i=0;i<TimeNodeNum;++i)
	{
		V_Copy(x1, x1_original, 14);
		// 积分到第i个时间节点
		flag = ode45(dynamics_rv_fop, x1, dfpara, 0.0, TimeNode[i], 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
		norm_lambdaV = V_Norm2(&x1[10], 3); // 速度协态的范数
		rou=1.0-(Ispg0NU*norm_lambdaV/x1[6] + x1[13])/lam01; // 开关函数
		if (rou > epsi)
			u = 0.0;
		else if (rou < -epsi)
			u = 1.0;
		else
			u = 0.5 - rou/(2*epsi);
		fprintf(ufid, "%.6f, %.6f\n", TimeNode[i], u);
	}
	fclose(ufid);
	delete[] TimeNode;
	*/

	// 第二段
	time_t now = time(NULL);
	char filename[30];
	sprintf(filename, "u_%d.txt", now);
	FILE *ufid = fopen(filename, "w");
	FILE *fid = NULL;

	int TimeNodeNum = 1000;
	double* TimeNode = new double[TimeNodeNum];
	double seg = t2/TimeNodeNum;
	for (int i=0;i<TimeNodeNum;++i)
		TimeNode[i] = seg*i;
	TimeNode[TimeNodeNum-1] = t2;

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam02;

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	double norm_lambdaV, rou, u;

	for (int i=0;i<TimeNodeNum;++i)
	{
		V_Copy(x2, x2_original, 14);
		// 积分到第i个时间节点
		flag = ode45(dynamics_rv_fop, x2, dfpara, 0.0, TimeNode[i], 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
		norm_lambdaV = V_Norm2(&x2[10], 3); // 速度协态的范数
		rou=1.0-(Ispg0NU*norm_lambdaV/x2[6] + x2[13])/lam02; // 开关函数
		if (rou > epsi)
			u = 0.0;
		else if (rou < -epsi)
			u = 1.0;
		else
			u = 0.5 - rou/(2*epsi);
		fprintf(ufid, "%.6f, %.6f\n", TimeNode[i], u);
	}
	fclose(ufid);

	delete[] TimeNode;
}