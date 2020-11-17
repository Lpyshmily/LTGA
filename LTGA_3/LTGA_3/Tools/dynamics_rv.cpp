#include "dynamics_rv.h"
#include "vector_operation.h"
#include "constants.h"

// �����õ�dfpara��ֵ�����ó�NULL����
int dynamics_rv_top(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);

	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double u = 1.0; // ʱ����������£�����ʼ������
	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3��λ�÷����ĵ�����3���ٶȷ���
		dx[i+3] = -muNU/(r*r*r)*rv[i] + u*TmaxNU/m*alpha[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = -u*TmaxNU/Ispg0NU;
	dx[13] = -u*TmaxNU/(m*m)*norm_lambdaV;
	return 1;
}

// dfpara[2];[0]~epsi;[1]~lam0
int dynamics_rv_fop(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);

	double epsi, lam0;
	epsi = dfpara[0];
	lam0 = dfpara[1];

	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou=1.0-(Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/(2*epsi);
	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3��λ�÷����ĵ�����3���ٶȷ���
		dx[i+3] = -muNU/(r*r*r)*rv[i] + u*TmaxNU/m*alpha[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = -u*TmaxNU/Ispg0NU;
	dx[13] = -u*TmaxNU/(m*m)*norm_lambdaV;
	return 1;
}