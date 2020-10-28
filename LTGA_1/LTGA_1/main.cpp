#include <iostream>
#include <ctime>
#include "localConst.h"
#include "GA_fop.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // ������f

void test_GA_fop()
{
	int flag, i;

	// **********
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
		// rv0[i] = rv0[i]/D2PI;
		// rv1[i] = rv1[i]/D2PI;
	}
	double rvm[6] = {0.0};
	Mars.GetRV(rvm, 59534, muNU);
	double m0 = 20000/MUnit;

	// ����㷨��һЩ��������
	int MaxGuessNum = 100;//�����������²����
	srand( (unsigned)time( NULL ) );//�趨��������ӣ���û�д����ã�ÿ�β���һ���������
	// ���
	double Out[18] = {0.0}; // �����������0-ʣ��������1~8-8����Ҫ��е�Э̬��ֵ
	double tf = 2201*86400/TUnit, epsi = 1.0;
	flag = solve_GA_fop(Out, rv0, rv1, rvm, m0, tf, epsi, MaxGuessNum);
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
	printf("��б���ֵΪ:\n");
	for (i=1; i<=17; i++)
		printf("%.15e,\n", Out[i]);
	printf("��������ʱ��Ϊ%f��\n", Out[17]*TUnit/86400);
}

int main()
{
	printf("Hello!\n");

	clock_t start, stop;
	start = clock();
	
	test_GA_fop();

	stop = clock();
	printf("������ʱΪ��%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}