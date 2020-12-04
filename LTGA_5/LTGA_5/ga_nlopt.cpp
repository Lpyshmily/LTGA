#include "ga_nlopt.h"
#include "nlopt.h"
#include "localConst.h"

#pragma comment(lib,"nlopt.lib")//��������ľ�̬�⣬ͨ���˾�̬����ö�̬��nlopt.dll

int global_count = 0;

// x[6] 0-��������ʱ�� 1-���������뾶 2-vout����� 3-5 ��������ǰ����ٶ�vin��3������
// �������ȼ�����Ž�������
// ���������ı仯��Χ��Ϊ0-1����Ӧ���²�ͬ����
// ��������ʱ��ָ����������ǰ��һ�ι켣��ת��ʱ�䣬�仯��ΧΪ[0.3,0.4]tf
// ���������뾶�仯��ΧΪ[1,2]rmin
// vout�����ָvout��vin��vp��ɵ�ƽ��ļнǣ��仯��ΧΪ[0,pi]
// vin3�������ı仯��Χ��Ϊ[-0.2,0.2]
// para NULL
double GA_obj_nlopt(unsigned n, const double* x, double* grad, void* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.1 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	double vin[3] = {x[3]*0.4-0.2, x[4]*0.4-0.2, x[5]*0.4-0.2}; // ��������ǰ������ٶ�
	
	int i, j, flag;
	
	
	// nloptʱ�����PSOʱ�����
	printf("**********\n");
	printf("�������� i=%d\n", ++global_count);
	for (i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");
	
	
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0}, rv_before[6] = {0.0}, rv_after[6] = {0.0}; // ����λ���ٶȡ���������ǰ��λ���ٶȡ������������λ���ٶ�
	double Out1[10] = {0.0}; // ��һ��ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out2[9] = {0.0}; // ��һ��ȼ�����Ž�����������[0]ʣ��������[1-8]8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0}; // �ڶ���ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out4[9] = {0.0}; // �ڶ���ȼ�����Ž�����������[0]ʣ��������[1-8]8����б���
	double m0, tempm, tf, epsi, t1, t2, shortest1, shortest2;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	m0 = 20000.0/MUnit;
	const int RepeatTime = 10; // ʱ�������ظ�������
	
	tempm = m0;
	shortest1 = MaxNum;// �������ó�һ���ܴ��ֵ
	shortest2 = MaxNum;// �������ó�һ���ܴ��ֵ
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);

	V_Copy(rv_before, rvm, 3);
	for (i=0;i<3;++i)
		rv_before[i+3] = rvm[i+3] + vin[i];


	// ����㷨��һЩ��������
	int MaxGuessNum = 500;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

	// ���
	// ��һ��
	// ʱ�����Ž���
	for (j=0;j<RepeatTime;++j)
	{
		// printf("��%d������һ��ʱ�����Ž�������\n", j+1);
		flag = solve_rv_top_rend_fixed(Out1, rv0, rv_before, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out1[0]*MUnit);
		printf("ת��ʱ��Ϊ:%.3f��\n", Out1[9]*TUnit/86400);
		printf("��б���ֵΪ:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out1[i]);
		*/
		if (Out1[9] < shortest1)
			shortest1 = Out1[9];
	}
	// printf("���ת��ʱ��Ϊ:%.3f��\n", shortest1*TUnit/86400);
	// �ж��ܷ���ɵ�һ��ת��
	if (shortest1 > t1)
	{
		// printf("�޷���ɵ�һ�ι켣ת��\n");
		// fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
		return MaxNum;
	}
	// ȼ�����Ž���
	flag = solve_rv_fop_rend(Out2, rv0, rv_before, tempm, t1, epsi, MaxGuessNum);
	// flag = homotopy_rv_fop_rend(Out2, rv0, rv_before, tempm, t1, MaxGuessNum);
	if (!flag)
	{
		// printf("��һ��ȼ�����Ž������ⲻ����\n");
		return MaxNum;
	}
	
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out2[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out2[i]);
	
	

	// ��������
	double vout[3] = {0.0}, vm[3] = {0.0}, unit1[3] = {0.0}, unit2[3] = {0.0}, unit3[3] = {0.0}, tempVec[3] = {0.0};
	double delta, norm_vin, temp;
	V_Copy(vm, &rvm[3], 3); // �����ٶ�
	// ��������ǰ������ٶ�vin��Ϊ�Ż��������ں���һ��ʼ����

	norm_vin = V_Norm2(vin, 3);
	V_Divid(unit1, vin, norm_vin, 3); // ��vin����ĵ�λʸ��
	V_Cross(tempVec, vin, vm);
	temp = V_Norm2(tempVec, 3);
	V_Divid(unit3, tempVec, temp, 3); // ��ֱvin��vmƽ��ĵ�λ����
	V_Cross(unit2, unit3, unit1); // j��λʸ��

	delta = 2*asin(muNU_MARS/(muNU_MARS + rp*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] + sin(delta)*sin(phi)*unit2[i] + sin(delta)*cos(phi)*unit3[i]);
	}
	// ������֪���Ž�ʱ��Ҫ
	/*
	vout[0] = 0.11925935378439018; vout[1] = -0.016224539415498812; vout[2] = 0.0028024583987117900;
	double temp1 = V_Dot(vout, unit1, 3)/norm_vin/cos(delta);
	double temp2 = V_Dot(vout, unit2, 3)/norm_vin/sin(delta);
	double temp3 = V_Dot(vout, unit3, 3)/norm_vin/sin(delta);
	*/
	


	// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
	V_Copy(rv_after, rvm, 3); // �����������λ��
	V_Add(&rv_after[3], vm, vout, 3); // ������������ٶ�
	tempm = Out2[0];// ���������������
	// ʱ�����Ž���
	for (j=0;j<RepeatTime;++j)
	{
		// printf("��%d�����ڶ���ʱ�����Ž�������\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, rv_after, rv1, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		/*
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out3[0]*MUnit);
		printf("ת��ʱ��Ϊ:%.3f��\n", Out3[9]*TUnit/86400);
		printf("��б���ֵΪ:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out3[i]);
			*/
		if (Out3[9] < shortest2)
			shortest2 = Out3[9];
	}
	// printf("���ת��ʱ��Ϊ:%.3f��\n", shortest2*TUnit/86400);
	// �ж��ܷ���ɵڶ���ת��
	if (shortest2 > t2)
	{
		// printf("ʣ��ʱ���޷���ɵڶ��ι켣ת��\n");
		return MaxNum;
	}
	// ȼ�����Ž���
	flag = solve_rv_fop_rend(Out4, rv_after, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
	{
		// printf("�ڶ���ȼ�����Ž������ⲻ����\n");
		return MaxNum;
	}
	
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	

	// nloptʱ�����PSOʱ�����
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);

	// output��Ҫ������
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
	printf("GA_obj_nlopt����ֵ��%f\n", GA_obj_nlopt(6, x, NULL, NULL));
}

void test_GA_obj_PSO()
{
	double x[6] = {0.863655768384999,0.000017319770127,0.499034175810391,0.621295996122853,0.228521973614765,0.515503534531755}; // 16021.090445199869
	printf("GA_obj_PSO����ֵ��%f\n", GA_obj_PSO(x, NULL));
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
	printf("��СֵΪ%.15f\n", f_min);
}

void output_ode()
{
	// ��Ҫ��ǰ��rv_fop_rend.cpp�ļ����ƽ��иĶ�
	// ��Ҫȷ����б���
	printf("������ֹ�����Ϣ\n");
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
	double Out[9] = {0.0}; // ȼ�����Ž�����������[0]ʣ��������[1-8]8����Ҫ��е�Э̬��ֵ

	// ��һ��
	// flag = solve_rv_fop_rend(Out, rv0, rv_before, m1, t1, epsi, MaxGuessNum);
	// �ڶ���
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

	// ��һ��
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
		// ���ֵ���i��ʱ��ڵ�
		flag = ode45(dynamics_rv_fop, x1, dfpara, 0.0, TimeNode[i], 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
		norm_lambdaV = V_Norm2(&x1[10], 3); // �ٶ�Э̬�ķ���
		rou=1.0-(Ispg0NU*norm_lambdaV/x1[6] + x1[13])/lam01; // ���غ���
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

	// �ڶ���
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
		// ���ֵ���i��ʱ��ڵ�
		flag = ode45(dynamics_rv_fop, x2, dfpara, 0.0, TimeNode[i], 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
		norm_lambdaV = V_Norm2(&x2[10], 3); // �ٶ�Э̬�ķ���
		rou=1.0-(Ispg0NU*norm_lambdaV/x2[6] + x2[13])/lam02; // ���غ���
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