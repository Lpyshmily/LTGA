#include "rv_fop_rend.h"
#include "vector_operation.h"
#include "ODE45.h"
#include "dynamics_rv.h"
#include "constants.h"
#include "MinPack\cminpack.h"

// n ��б������� 8
/* x[8] ��б���
[0] lam0
[1-7] Э̬���� */
/* fvec[8] ���ƫ�� �˺����������
[0-5] ĩ��״̬Լ��
[6] ĩ������Э̬Ϊ��
[7] Э̬��һ�� */
// iflagΪ0ʱ��ʾ�Ϸ����ú�����֪������в�
/* sfpara[16]
[0-5] rv0
[6] m0
[7-12] rv1 Ŀ�������λ���ٶ�
[13] ת��ʱ��
[14] epsilon
[15] �Ƿ�����ı�־��0��ʾ�������1��ʾ��� */
// ����ɹ�ʱ������0�����0��ֵ�����ɹ�ʱ����С��0��ֵ
int fvec_rv_fop_rend(int n, const double* x, double* fvec, int iflag, const double* sfpara)
{
	int i;

	double x0[14] = {0.0};
	double rvf[6] = {0.0};

	V_Copy(x0, sfpara, 7); // ״̬��ֵ
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	double epsi = sfpara[14];
	int outflag = NINT(sfpara[15]);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7); // Э̬��ֵ

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;

	if(iflag==0)
	{
		return 0;
	}

	double AbsTol[14]={0.0};
	for(i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;
	double work[140]={0.0};
	FILE *fid=NULL;//fopen("temp.txt","w");//����趨��Ч�ļ�·���������Ҫ�ر��ļ�
	int flag, NumPoint;
	flag = ode45(dynamics_rv_fop, x0, dfpara, 0.0, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);

	V_Minus(fvec, x0, rvf, 6);
	fvec[6] = x0[13];
	fvec[7] = V_Norm2(x, 8) - 1.0;

	if (outflag>0)
	{
		fvec[0] = x0[6]; // ���ʣ������
	}
	// fclose(fid);
	return 0;
}
// Out[9] ������� [0]~ʣ������ [1-8]~8����б���
int solve_rv_fop_rend(double* Out, const double* rv0, const double* rvf, double m0, double tf, double epsi, int MaxGuessNum)
{
	int n, j, info, flag=0;
	
	double sfpara[16] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rvf, 6);
	sfpara[13] = tf;
	sfpara[14] = epsi;
	sfpara[15] = 0.0;

	n = 8; // ��б�������
	double x[8] = {0.0};
	double fvec[8] = {0.0};

	Out[0] = 0.0;
	int num = 0;
	double wa[200]={0.0};
	double xtol=1.0e-8;
	double angle[7];
	while(num<=MaxGuessNum)
	{
		for(j=0;j<7;j++)
				angle[j] = (double)rand()/RAND_MAX;
		angle[0] = M_PI2*angle[0];
		angle[1] = M_PI2*angle[1];
		angle[2] = M_PI2*angle[2];
		angle[3] = M_PI*(angle[3] - 0.5);
		angle[4] = M_PI*(angle[4] - 0.5);
		angle[5] = M_2PI*angle[5];
		angle[6] = M_2PI*angle[6];

		x[0] = sin(angle[0]);
		x[1] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*cos(angle[5]);
		x[2] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*sin(angle[5]);
		x[3] = cos(angle[0])*cos(angle[1])*cos(angle[2])*sin(angle[3]);
		x[4] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*cos(angle[6]);
		x[5] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*sin(angle[6]);
		x[6] = cos(angle[0])*cos(angle[1])*sin(angle[2])*sin(angle[4]);
		x[7] = cos(angle[0])*sin(angle[1]);
		
		info = hybrd1(fvec_rv_fop_rend, n, x, fvec, sfpara, wa, xtol, 20, 2000);
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			j=fvec_rv_fop_rend(n, x, fvec, 1, sfpara);
			if(fvec[0]>0.0) // ʣ����������Ϊ��
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 8);
				break;
			}
			sfpara[15]=0.0;
		}
		num++;
	}
	printf("����²������%d\n", num);
	return flag;
}

// �ⲿ�ֳ������ƣ���ʱ��ʹ��
/*
// ͬ�׹��̣���0.1ͬ�׵�1e-5
// ����0.1ʱ��������²����ɹ�
int homotopy_rv_fop_rend(double* Out, const double* rv0, const double* rvf, double m0, double tf, int MaxGuessNum)
{
	int n, j, info, flag=0;

	double epsi = 0.1, depsi = 0.01;
	
	double sfpara[16] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rvf, 6);
	sfpara[13] = tf;
	sfpara[14] = epsi;
	sfpara[15] = 0.0;

	n = 8; // ��б�������
	double x[8] = {0.0};
	double fvec[8] = {0.0};

	Out[0] = 0.0;
	int num = 0;
	double wa[200]={0.0};
	double xtol=1.0e-8;
	double angle[7];
	while(num<=MaxGuessNum)
	{
		for(j=0;j<7;j++)
				angle[j] = (double)rand()/RAND_MAX;
		angle[0] = M_PI2*angle[0];
		angle[1] = M_PI2*angle[1];
		angle[2] = M_PI2*angle[2];
		angle[3] = M_PI*(angle[3] - 0.5);
		angle[4] = M_PI*(angle[4] - 0.5);
		angle[5] = M_2PI*angle[5];
		angle[6] = M_2PI*angle[6];

		x[0] = sin(angle[0]);
		x[1] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*cos(angle[5]);
		x[2] = cos(angle[0])*cos(angle[1])*cos(angle[2])*cos(angle[3])*sin(angle[5]);
		x[3] = cos(angle[0])*cos(angle[1])*cos(angle[2])*sin(angle[3]);
		x[4] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*cos(angle[6]);
		x[5] = cos(angle[0])*cos(angle[1])*sin(angle[2])*cos(angle[4])*sin(angle[6]);
		x[6] = cos(angle[0])*cos(angle[1])*sin(angle[2])*sin(angle[4]);
		x[7] = cos(angle[0])*sin(angle[1]);
		
		info = hybrd1(fvec_rv_fop_rend, n, x, fvec, sfpara, wa, xtol, 20, 2000);
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			j=fvec_rv_fop_rend(n, x, fvec, 1, sfpara);
			if(fvec[0]>0.0) // ʣ����������Ϊ��
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 8);
				break;
			}
			sfpara[15]=0.0;
		}
		num++;
	}
	printf("epsi=%f���ɹ�������²������%d\n", epsi, num);

	// ͬ�׹���
	sfpara[15]=0.0;
	while (epsi > 1e-6)
	{
		flag = 0;
		// if (epsi < depsi + 1e-10)
			// depsi *= 0.1;
		if (depsi < 1e-10)
			break;
		epsi -= depsi;
		sfpara[14] = epsi;
		info = hybrd1(fvec_rv_fop_rend, n, x, fvec, sfpara, wa, xtol, 20, 2000);
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			printf("ͬ�׹��̳ɹ���epsi=%.15f\n", epsi);
			sfpara[15]=1.0;
			j=fvec_rv_fop_rend(n, x, fvec, 1, sfpara);
			if(fvec[0]>0.0) // ʣ����������Ϊ��
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 8);
			}
			sfpara[15]=0.0;
		}
		else
		{
			printf("ͬ�׹���ʧ��\n");
			epsi += depsi;
			depsi *= 0.1;
		}
	}

	return flag;
}
*/
