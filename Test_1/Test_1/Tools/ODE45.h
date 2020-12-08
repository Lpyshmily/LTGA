#ifndef _ODE45_H_
#define _ODE45_H_
#include<math.h>
#include<iostream>
#include<stdio.h>

//ode45��΢�ַ���fun(xΪ״̬��, dxΪ����, paraΪ����Ĳ�����)�䲽������
//NumPointΪ�䲽�����ֵĵ���, workΪ��Ҫ���ܵĹ������飬��СΪά��n��10����AbsTolΪ���������ֵ��ά��Ϊn��RelToΪ������
//NormControlΪ1ʱ��ʾ��Ҫ���д�����ƣ�MaxStep��ʾ�������������ƣ�ȡ��ֵʱ��ʾ�����ƣ�InitialStep��ʾ��ʼ�²ⲽ����ȡ��ʱ��ʾ�Զ��²ⲽ��
//�䲽�����ֵ��м������ļ�fid�У���fid=NULL��ʾ�������м�������ֻ�������ս����x��
int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol=1.0e-3,
		  int NormControl=0, double MaxStep=-1, double InitialStep=-1, FILE* fid=NULL);

//���ź���
template<class T> inline int sign(const T & x)
{
	if(x>0)
		return 1;
	else if(x<0)
		return -1;
	else
		return 0;
}

//�����ֵ
template <class T>
inline T max(T x, T y) 
{ 
	return (x>y)?x:y;
}


//����Сֵ
template <class T>
inline T min(T x, T y) 
{
	return (x<y)?x:y;
}

#endif