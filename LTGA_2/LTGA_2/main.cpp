#include <iostream>
#include <ctime>
#include "localConst.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // ������f

void test_ee_fop_rend()
{
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // ����ʱ��Ϊ1000�죬���й�һ��
	// ��ʼλ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// ĩ��λ�ú��ٶ�,����λ�ֱ�ΪAU��AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// ���ˣ����г�ʼ������һ�����

	
	// ���ô��ֵ��������������
	int flag;
	// ��λ���ٶ�ת��Ϊ���ֵ�������
	double ee0[6]={0.0}, ee1[6]={0.0};
	rv2ee(flag, ee0, rv0, muNU);
	rv2ee(flag, ee1, rv1, muNU);
	ee1[5] = ee1[5] + 3*M_2PI;
	// ����㷨��һЩ��������
	double epsi = 1.0;//ȡ��һ����С��ͬ�ײ���ֱ�������ư����ƵĽ��
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������
	// ���
	double Out[9] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	flag = solve_ee_fop_rend(Out, ee0, ee1, m0, tof, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0Ϊ:%.6e\n", Out[1]);
	printf("7����ʼЭ̬����ֵΪ:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_fop_rend()
{
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // ����ʱ��Ϊ1000�죬���й�һ��
	// ��ʼλ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// ĩ��λ�ú��ٶ�,����λ�ֱ�ΪAU��AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// ���ˣ����г�ʼ������һ�����

	
	// ���ô��ֵ��������������
	int flag;
	// ����㷨��һЩ��������
	double epsi = 1.0;//ȡ��һ����С��ͬ�ײ���ֱ�������ư����ƵĽ��
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������
	// ���
	double Out[9] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	flag = solve_rv_fop_rend(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0Ϊ:%.6e\n", Out[1]);
	printf("7����ʼЭ̬����ֵΪ:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_top_rend_fixed()
{
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	// ��ʼλ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// ĩ��λ�ú��ٶ�,����λ�ֱ�ΪAU��AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// ���ˣ����г�ʼ������һ�����

	
	// ���ô��ֵ��������������
	int flag;
	// ����㷨��һЩ��������
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������
	// ���
	double Out[10] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	flag = solve_rv_top_rend_fixed(Out, rv0, rv1, m0, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
	printf("ת��ʱ��Ϊ:%.3f��\n", Out[9]*TUnit/86400);
	printf("7����ʼЭ̬����ֵΪ:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
}

void test_rv_fop_flyby()
{
	int i;
	
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=3800 Tmax=0.33 m0=1500
	double m0 = 1500.0/MUnit;
	double tof = 1000.0*86400.0/TUnit; // ����ʱ��Ϊ1000�죬���й�һ��
	// ��ʼλ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// ĩ��λ�ú��ٶ�,����λ�ֱ�ΪAU��AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/M_2PI;
		rv1[i] = rv1[i]/M_2PI;
	}
	// ���ˣ����г�ʼ������һ�����

	
	// ���ô��ֵ��������������
	int flag;
	// ����㷨��һЩ��������
	double epsi = 1.0;//ȡ��һ����С��ͬ�ײ���ֱ�������ư����ƵĽ��
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������
	// ���
	double Out[15] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	flag = solve_rv_fop_flyby(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[6]*MUnit);
	printf("ĩ��״̬��Ϊ:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out[i]);
	printf("��б���ֵΪ:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out[i]);
}

void test_GA()
{
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=6000 Tmax=2.26 m0=20000
	
	int i, flag;
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
	}
	double rvm[6] = {0.0};
	double Out1[15] = {0.0}; // �����������[0-6]-ĩ��״̬��[7-14]-8����Ҫ��е�Э̬��ֵ
	double Out2[9] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0};
	double m0, tf, epsi, factor, t1, t2;
	tf = 2201*86400/TUnit;
	epsi = 1.0;
	m0 = 20000.0/MUnit;

	// �����������һ�������������һ�ι켣ת��ʱ������ʱ���е�ռ��
	factor = 0.37;
	t1 = factor*tf;
	t2 = (1-factor)*tf;
	Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


	// ����㷨��һЩ��������
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

	// ���
	// ��һ�� ȼ�����ŷ�Խ
	flag = solve_rv_fop_flyby(Out1, rv0, rvm, m0, t1, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out1[6]*MUnit);
	printf("ĩ��״̬��Ϊ:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out1[i]);
	printf("��б���ֵΪ:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out1[i]);
	

	// ��������
	double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, tempV1[3] = {0.0}, tempV2[3] = {0.0};
	double temp, norm_vm;

	V_Copy(vm, &rvm[3], 3); // �����ٶ�
	norm_vm = V_Norm2(vm, 3);
	V_Minus(vin, &Out1[3], vm, 3); // ��������ǰ������ٶ�
	temp = V_Dot(vin, vm, 3)/(norm_vm*norm_vm);
	V_Multi(tempV1, vm, temp, 3); // ����ٶ���vm�ķ���
	V_Minus(tempV2, vin, tempV1, 3); // ����ٶȴ�ֱ��vm�ķ���
	V_Minus(vout, tempV2, tempV1, 3); // ���������������ٶ�

	// �ڶ��� ��rv0��ʾ�����������λ���ٶ�
	V_Copy(rv0, rvm, 3); // �����������λ��
	V_Add(&rv0[3], vm, vout, 3); // ������������ٶ�
	m0 = Out1[6];// ���������������
	// ʱ�����Ž���
	
	flag = solve_rv_top_rend_fixed(Out3, rv0, rv1, m0, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out3[0]*MUnit);
	printf("ת��ʱ��Ϊ:%.3f��\n", Out3[9]*TUnit/86400);
	printf("��б���ֵΪ:\n");
	for (i=1; i<10; i++)
		printf("%.15e,\n", Out3[i]);

	// �ж��ܷ���ɵڶ���ת��
	if (Out3[9] > t2)
	{
		printf("ʣ��ʱ���޷���ɵڶ��ι켣ת��\n");
		return;
	}

	// ȼ�����Ž���
	flag = solve_rv_fop_rend(Out2, rv0, rv1, m0, t2, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out2[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out2[i]);
}

void test_GA_factor()
{
	// ��Ҫ��Tools.h/constants.h�еĲ��������޸�
	// Isp=6000 Tmax=2.26 m0=20000

	FILE *fid = fopen("info_5.txt", "w");
	
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
	double rvm[6] = {0.0}, temprv[6] = {0.0};
	double Out1[10] = {0.0}; // ��һ��ʱ�����ŷ�Խ��������[0]ʣ��������[1-9]9����б���
	double Out2[15] = {0.0}; // ��һ��ȼ�����ŷ�Խ��������[0-6]ĩ��״̬��[7-14]8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0}; // �ڶ���ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out4[9] = {0.0}; // �ڶ���ȼ�����Ž�����������[0]ʣ��������[1-8]8����б���
	double m0, tempm, tf, epsi, factor, t1, t2, shortest1, shortest2;
	tf = 2201*86400/TUnit;
	epsi = 1.0e-5;
	m0 = 20000.0/MUnit;
	const int RepeatTime = 10; // ʱ�������ظ�������

	// �����������һ�������������һ�ι켣ת��ʱ������ʱ���е�ռ��
	factor = 0.34;
	while (factor<=0.36)
	{
		printf("**********\nfractor=%f\n", factor);
		tempm = m0;
		shortest1 = 1e7;// �������ó�һ���ܴ��ֵ
		shortest2 = 1e7;// �������ó�һ���ܴ��ֵ
		t1 = factor*tf;
		t2 = (1-factor)*tf;
		Mars.GetRV(rvm, 59534.0 + t1*TUnit/86400, muNU);


		// ����㷨��һЩ��������
		int MaxGuessNum = 100;//�����������²����
		srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

		// ���
		// ��һ��
		// ʱ�����ŷ�Խ
		for (j=0;j<RepeatTime;++j)
		{
			printf("��%d������һ��ʱ�����ŷ�Խ����\n", j+1);
			flag = solve_rv_top_flyby_fixed(Out1, rv0, rvm, tempm, MaxGuessNum);
			printf("���ɹ�%d\n",flag);
			printf("ʣ������Ϊ:%.3fkg\n", Out1[0]*MUnit);
			printf("ת��ʱ��Ϊ:%.3f��\n", Out1[9]*TUnit/86400);
			printf("��б���ֵΪ:\n");
			for (i=1; i<10; i++)
				printf("%.15e,\n", Out1[i]);
			if (flag && (Out1[9] < shortest1))
				shortest1 = Out1[9];
		}
		printf("���ת��ʱ��Ϊ:%.3f��\n", shortest1*TUnit/86400);
		// �ж��ܷ���ɵ�һ��ת��
		if (shortest1 > t1)
		{
			printf("�޷���ɵ�һ�ι켣ת��\n");
			fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
			factor += 0.001;
			continue;
		}
		// ȼ�����ŷ�Խ
		flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out2[6]*MUnit);
		printf("ĩ��״̬��Ϊ:\n");
		for (i=0;i<7;++i)
			printf("%.15e,\n", Out2[i]);
		printf("��б���ֵΪ:\n");
		for (i=7; i<15; i++)
			printf("%.15e,\n", Out2[i]);
	

		// ��������
		double vin[3] = {0.0}, vout[3] = {0.0}, vm[3] = {0.0}, tempV1[3] = {0.0}, tempV2[3] = {0.0};
		double temp, norm_vm;

		V_Copy(vm, &rvm[3], 3); // �����ٶ�
		norm_vm = V_Norm2(vm, 3);
		V_Minus(vin, &Out2[3], vm, 3); // ��������ǰ������ٶ�
		temp = V_Dot(vin, vm, 3)/(norm_vm*norm_vm);
		V_Multi(tempV1, vm, temp, 3); // ����ٶ���vm�ķ���
		V_Minus(tempV2, vin, tempV1, 3); // ����ٶȴ�ֱ��vm�ķ���
		V_Minus(vout, tempV2, tempV1, 3); // ���������������ٶ�

		// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
		V_Copy(temprv, rvm, 3); // �����������λ��
		V_Add(&temprv[3], vm, vout, 3); // ������������ٶ�
		tempm = Out2[6];// ���������������
		// ʱ�����Ž���
		for (j=0;j<RepeatTime;++j)
		{
			printf("��%d�����ڶ���ʱ�����Ž�������\n", j+1);
			flag = solve_rv_top_rend_fixed(Out3, temprv, rv1, tempm, MaxGuessNum);
			printf("���ɹ�%d\n",flag);
			printf("ʣ������Ϊ:%.3fkg\n", Out3[0]*MUnit);
			printf("ת��ʱ��Ϊ:%.3f��\n", Out3[9]*TUnit/86400);
			printf("��б���ֵΪ:\n");
			for (i=1; i<10; i++)
				printf("%.15e,\n", Out3[i]);
			if (flag && (Out3[9] < shortest2))
				shortest2 = Out3[9];
		}
		printf("���ת��ʱ��Ϊ:%.3f��\n", shortest2*TUnit/86400);
		// �ж��ܷ���ɵڶ���ת��
		if (shortest2 > t2)
		{
			printf("ʣ��ʱ���޷���ɵڶ��ι켣ת��\n");
			fprintf(fid, "%f\t%f\t%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400, t2*TUnit/86400, shortest2*TUnit/86400);
			factor += 0.001;
			continue;
		}
		// ȼ�����Ž���
		flag = solve_rv_fop_rend(Out4, temprv, rv1, tempm, t2, epsi, MaxGuessNum);
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
		printf("��б���ֵΪ:\n");
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
	printf("������ʱΪ��%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}