#include "ga_nlopt_ee.h"
#include "localConst.h"

int global_count_ee = 0;

// x[6] 0-��������ʱ�� 1-���������뾶 2-vout����� 3-5 ��������ǰ����ٶ�vin��3������
// �������ȼ�����Ž�������
// ���������ı仯��Χ��Ϊ0-1����Ӧ���²�ͬ����
// ��������ʱ��ָ����������ǰ��һ�ι켣��ת��ʱ�䣬�仯��ΧΪ[0.3,0.35]tf
// ���������뾶�仯��ΧΪ[1,2]rmin
// vout�����ָvout��vin��vp��ɵ�ƽ��ļнǣ��仯��ΧΪ[0,pi]
// vin3�������ı仯��Χ��Ϊ[-0.1,0.1]
// para NULL
double GA_obj_nlopt_6_ee(unsigned n, const double* x, double* grad, void* para)
{
	// Tools.h/constants.h
	// Isp=6000 Tmax=2.26 m0=20000

	double factor = x[0]*0.1 + 0.3;
	double rp = rminU_MARS*(x[1] + 1.0);
	double phi = x[2]*M_PI;
	double vin[3] = {x[3]*0.2-0.1, x[4]*0.2-0.1, x[5]*0.2-0.1}; // ��������ǰ������ٶ�
	
	int i, j, flag;
	printf("**********\n");
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
	double ee0[6], ee1[6], ee_before[6], ee_after[6];
	double Out1[10] = {0.0}; // ��һ��ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out2[9] = {0.0}; // ��һ��ȼ�����Ž�����������[0]ʣ��������[1-8]8����Ҫ��е�Э̬��ֵ
	double Out3[10] = {0.0}; // �ڶ���ʱ�����Ž�����������[0]ʣ��������[1-9]9����б���
	double Out4[9] = {0.0}; // �ڶ���ȼ�����Ž�����������[0]ʣ��������[1-8]8����б���
	double m0, tempm, tf, epsi, t1, t2, shortest1, shortest2;
	tf = 2201*86400/TUnit;
	epsi = 1.0;
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

	rv2ee(flag, ee0, rv0, muNU);
	rv2ee(flag, ee1, rv1, muNU);
	rv2ee(flag, ee_before, rv_before, muNU);


	// ����㷨��һЩ��������
	int MaxGuessNum = 1000;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������

	// ���
	// ��һ��
	// ʱ�����Ž���
	for (j=0;j<RepeatTime;++j)
	{
		printf("��%d������һ��ʱ�����Ž�������\n", j+1);
		flag = solve_ee_top_rend_fixed(Out1, ee0, ee_before, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out1[0]*MUnit);
		printf("ת��ʱ��Ϊ:%.3f��\n", Out1[9]*TUnit/86400);
		printf("��б���ֵΪ:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out1[i]);
		if (Out1[9] < shortest1)
			shortest1 = Out1[9];
	}
	printf("���ת��ʱ��Ϊ:%.3f��\n", shortest1*TUnit/86400);
	// �ж��ܷ���ɵ�һ��ת��
	if (shortest1 > t1)
	{
		printf("�޷���ɵ�һ�ι켣ת��\n");
		// fprintf(fid, "%f\t%f\t%f\n", factor, t1*TUnit/86400, shortest1*TUnit/86400);
		return MaxNum;
	}
	// ȼ�����Ž���
	flag = solve_ee_fop_rend(Out2, ee0, ee_before, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
	{
		printf("��һ��ȼ�����Ž������ⲻ����\n");
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
	


	// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
	V_Copy(rv_after, rvm, 3); // �����������λ��
	V_Add(&rv_after[3], vm, vout, 3); // ������������ٶ�
	rv2ee(flag, ee_after, rv_after, muNU);
	tempm = Out2[0];// ���������������
	// ʱ�����Ž���
	for (j=0;j<RepeatTime;++j)
	{
		printf("��%d�����ڶ���ʱ�����Ž�������\n", j+1);
		flag = solve_ee_top_rend_fixed(Out3, ee_after, ee1, tempm, MaxGuessNum);
		if (!flag)
			return MaxNum;
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out3[0]*MUnit);
		printf("ת��ʱ��Ϊ:%.3f��\n", Out3[9]*TUnit/86400);
		printf("��б���ֵΪ:\n");
		for (i=1; i<10; i++)
			printf("%.15e,\n", Out3[i]);
		if (Out3[9] < shortest2)
			shortest2 = Out3[9];
	}
	printf("���ת��ʱ��Ϊ:%.3f��\n", shortest2*TUnit/86400);
	// �ж��ܷ���ɵڶ���ת��
	if (shortest2 > t2)
	{
		printf("ʣ��ʱ���޷���ɵڶ��ι켣ת��\n");
		return MaxNum;
	}
	// ȼ�����Ž���
	flag = solve_rv_fop_rend(Out4, ee_after, ee1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
	{
		printf("�ڶ���ȼ�����Ž������ⲻ����\n");
		return MaxNum;
	}
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);

	printf("**********\n");
	printf("�������� i=%d\n", ++global_count_ee);
	for (i=0;i<6;++i)
		printf("%.15f,\n", x[i]);
	printf("**********\n");

	return -Out4[0]*MUnit;
}