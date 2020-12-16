#include "ga_single.h"
#include "nlopt.h"
#include "localConst.h"

#pragma comment(lib,"nlopt.lib")//��������ľ�̬�⣬ͨ���˾�̬����ö�̬��nlopt.dll

double ga_single1_obj(double x)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000
	
	int i, j, flag;
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0}, rv_after[6] = {0.0};
	double Out1[10] = {0.0}; // ��һ��ʱ�����ŷ�Խ��������[0]ʣ��������[1-9]9����б���
	double Out2[15] = {0.0}; // ��һ��ȼ�����ŷ�Խ��������[0-6]ĩ��״̬��[7-14]8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0}; // �ڶ���ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out4[9] = {0.0}; // �ڶ���ȼ�����Ž�����������[0]ʣ��������[1-8]8����б���
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	m0 = 20000.0/MUnit;
	tempm = m0;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	factor = x;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	shortest1 = MaxNum;// �������ó�һ���ܴ��ֵ
	shortest2 = MaxNum;// �������ó�һ���ܴ��ֵ

	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// ����㷨��һЩ��������
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

	// ���
	// ��һ��
	// ʱ�����ŷ�Խ
	for (j=0;j<RepeatTime;++j)
	{
		// printf("��%d������һ��ʱ�����ŷ�Խ����\n", j+1);
		flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
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
	// ȼ�����ŷ�Խ
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out2[6]*MUnit);
	printf("ĩ��״̬��Ϊ:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("��б���ֵΪ:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);
	*/
	
	// ��֪���������뾶ʱ����������
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, unit1[3] = {0.0}, unit2[3] = {0.0}, unit3[3] = {0.0}, tempVec[3] = {0.0};
	double delta, norm_vin, temp;
	V_Copy(vm, &rvm[3], 3); // �����ٶ�
	V_Minus(vin, &Out2[3], vm, 3); // ��������ǰ������ٶ�

	norm_vin = V_Norm2(vin, 3);
	V_Divid(unit1, vin, norm_vin, 3); // ��vin����ĵ�λʸ��
	V_Cross(tempVec, vin, vm);
	temp = V_Norm2(tempVec, 3);
	V_Divid(unit3, tempVec, temp, 3); // ��ֱvin��vmƽ��ĵ�λ����
	V_Cross(unit2, unit3, unit1); // j��λʸ��

	delta = 2*asin(muNU_MARS/(muNU_MARS + rminU_MARS*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] + sin(delta)*unit2[i] + 0.0*unit3[i]);
	}
	


	// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
	V_Copy(rv_after, rvm, 3); // �����������λ��
	V_Add(&rv_after[3], vm, vout, 3); // ������������ٶ�
	tempm = Out2[6];// ���������������
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
		return MaxNum;
	/*
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);
	*/

	return -Out4[0]*MUnit;
}

void ga_single1_obj_test()
{
	double x = 0.33;
	printf("%f��Ӧ��Ŀ�꺯��ֵΪ%f\n", x, ga_single1_obj(x));
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
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0}, rv_after[6] = {0.0};
	double Out1[10] = {0.0}; // ��һ��ʱ�����ŷ�Խ��������[0]ʣ��������[1-9]9����б���
	double Out2[15] = {0.0}; // ��һ��ȼ�����ŷ�Խ��������[0-6]ĩ��״̬��[7-14]8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0}; // �ڶ���ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out4[9] = {0.0}; // �ڶ���ȼ�����Ž�����������[0]ʣ��������[1-8]8����б���
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	m0 = 20000.0/MUnit;
	tempm = m0;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	factor = x;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	shortest1 = MaxNum;// �������ó�һ���ܴ��ֵ
	shortest2 = MaxNum;// �������ó�һ���ܴ��ֵ

	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// ����㷨��һЩ��������
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

	// ���
	// ��һ��
	// ʱ�����ŷ�Խ
	for (j=0;j<RepeatTime;++j)
	{
		// printf("��%d������һ��ʱ�����ŷ�Խ����\n", j+1);
		flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
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
	// ȼ�����ŷ�Խ
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
		return MaxNum;
	/*
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out2[6]*MUnit);
	printf("ĩ��״̬��Ϊ:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("��б���ֵΪ:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);
	*/
	
	// ��֪���������뾶ʱ����������
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, unit1[3] = {0.0}, unit2[3] = {0.0}, unit3[3] = {0.0}, tempVec[3] = {0.0};
	double delta, norm_vin, temp;
	V_Copy(vm, &rvm[3], 3); // �����ٶ�
	V_Minus(vin, &Out2[3], vm, 3); // ��������ǰ������ٶ�

	norm_vin = V_Norm2(vin, 3);
	V_Divid(unit1, vin, norm_vin, 3); // ��vin����ĵ�λʸ��
	V_Cross(tempVec, vin, vm);
	temp = V_Norm2(tempVec, 3);
	V_Divid(unit3, tempVec, temp, 3); // ��ֱvin��vmƽ��ĵ�λ����
	V_Cross(unit2, unit3, unit1); // j��λʸ��

	delta = 2*asin(muNU_MARS/(muNU_MARS + rminU_MARS*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] - sin(delta)*unit2[i] + 0.0*unit3[i]);
	}
	


	// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
	V_Copy(rv_after, rvm, 3); // �����������λ��
	V_Add(&rv_after[3], vm, vout, 3); // ������������ٶ�
	tempm = Out2[6];// ���������������
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
		return MaxNum;
	/*
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
	printf("��б���ֵΪ:\n");
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
	printf("��СֵΪ%.15f\n", f_min);
}
