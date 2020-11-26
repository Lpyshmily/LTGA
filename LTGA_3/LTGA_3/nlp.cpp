#include "nlp.h"
#include "localConst.h"
#include<iostream>
#include <stdlib.h>
#include"cnpsol.h"
#include"Use.h"

#pragma comment(lib,"cnpsol.lib")//��������ľ�̬�⣬ͨ���˾�̬����ö�̬��cnpsol.dll

using namespace std;

extern celestial_body Mars;

//����ʵ�����������ӣ��ֱ�Ϊmatlab��fmincon�е����Ӻ�NPSOL�е����ӣ�ִ���ĸ����ӿɽ���һ��������صĴ���ע�͵�

//�û�����ʵ�������д��ָ�꺯��
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
		if (!flag)
		{
			objf = MaxNum;
			return true;
		}
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
		objf = MaxNum;
		return true;
	}
	// ȼ�����ŷ�Խ
	flag = solve_rv_fop_flyby(Out2, rv0, rvm, tempm, t1, epsi, MaxGuessNum);
	if (!flag)
	{
		objf = MaxNum;
		return true;
	}
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out2[6]*MUnit);
	printf("ĩ��״̬��Ϊ:\n");
	for (i=0;i<7;++i)
		printf("%.15e,\n", Out2[i]);
	printf("��б���ֵΪ:\n");
	for (i=7; i<15; i++)
		printf("%.15e,\n", Out2[i]);


	// ��������
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

	delta = 2*asin(muNU_MARS/(muNU_MARS + rp*norm_vin*norm_vin));
	for (i=0;i<3;++i)
	{
		vout[i] = norm_vin*(cos(delta)*unit1[i] + sin(delta)*sin(phi)*unit2[i] + sin(delta)*cos(phi)*unit3[i]);
	}


	// �ڶ��� ��temprv��ʾ�����������λ���ٶ�
	V_Copy(temprv, rvm, 3); // �����������λ��
	V_Add(&temprv[3], vm, vout, 3); // ������������ٶ�
	tempm = Out2[6];// ���������������
	// ʱ�����Ž���
	for (j=0;j<RepeatTime;++j)
	{
		printf("��%d�����ڶ���ʱ�����Ž�������\n", j+1);
		flag = solve_rv_top_rend_fixed(Out3, temprv, rv1, tempm, MaxGuessNum);
		if (!flag)
		{
			objf = MaxNum;
			return true;
		}
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
		objf = MaxNum;
		return true;
	}
	// ȼ�����Ž���
	flag = solve_rv_fop_rend(Out4, temprv, rv1, tempm, t2, epsi, MaxGuessNum);
	if (!flag)
	{
		objf = MaxNum;
		return true;
	}
	printf("���ɹ�%d\n",flag);
	printf("**********\n");
	printf("x[0]=%f\n", x[0]);
	printf("x[1]=%f\n", x[1]);
	printf("x[2]=%f\n", x[2]);
	printf("ʣ������Ϊ:%.3fkg\n", Out4[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<9; i++)
		printf("%.15e,\n", Out4[i]);

	objf = -Out4[0]*MUnit;
	return true;
}

//�û�����ʵ�������д�ķ�����Լ������
bool fncon(int n, int ncnln, const double* x, double* nlc)
{
	return true;
}

void test_nlp()
{
	int i, j;
	int n=3, nclin=0, ncnln=0;
	int nctot=n+nclin+ncnln, inform=0;//nctotΪ�����Ա���������Լ����������Լ��������Լ���ĸ�����Ҳ��bu��bl��ά��
	//����һ��2ά����洢����Լ��ϵ������
	double** AM=NULL;//new double*[nclin];
	/*for(i=0;i<nclin;i++)
		AM[i]=new double[n];*/
	////����1ά����洢Լ���ϡ��½�
	double* bl=new double[nctot];
	double* bu=new double[nctot];
	double objf=0.0;

	double bigbnd=1.0E21;//�����Լ��������Լ��������ȡ���ϴ��ֵ���������x1>0��ȡ��Ӧ��bl=0,bu=bigbnd����Ϊ0<=x1<=bigbnd
	//for(i=0;i<nclin;i++)
	//	for(j=0;j<n;j++)
	//		AM[i][j]=0.0;
	
	//��ȫ����Ϊ����������൱����Լ�����ٶ�Լ���������
	for(j=0;j<nctot;j++)
	{
		bl[j]=-bigbnd;
		bu[j]=bigbnd;
	}
	//matlab��fmincon�е�����
	bl[0]=0.0;
	bu[0]=1.0;
	bl[1]=0.0;
	bu[1]=1.0;
	bl[2] = 0.0;
	bu[2] = 1.0;

	//��ֵ
	double x[3]={0.87704,0.000017,0.469045};
	inform=NPSol(x, objf, n, nclin, ncnln, AM, bl, bu);
	//������Ž����ָ�꼰����ɹ������̬��ʶ
	for(i=0;i<n;i++) printf("x[%d]=%.12e,\n",i,x[i]);
	printf("obj=%.12e,\n",objf);
	printf("inform=%d,\n",inform);	
	
	//�ͷŶ�̬���ڣ�������
	delete[] bu;
	delete[] bl;
	/*for(i=0;i<nclin;i++)
		delete[] AM[i];
	delete[] AM;*/
	
	printf("Hello Down!\n");
}