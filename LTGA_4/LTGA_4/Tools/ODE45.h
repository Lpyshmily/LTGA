#ifndef _ODE45_H_
#define _ODE45_H_
#include<math.h>
#include<iostream>
#include<stdio.h>

//ode45对微分方程fun(x为状态量, dx为导数, para为所需的参数集)变步长积分
//NumPoint为变步长积分的点数, work为需要开避的工作数组，大小为维数n的10倍，AbsTol为绝对误差数值，维数为n，RelTo为相对误差
//NormControl为1时表示需要进行此项控制，MaxStep表示最大迭代次数限制，取负值时表示不限制，InitialStep表示初始猜测步长，取负时表示自动猜测步长
//变步长积分的中间结果在文件fid中，若fid=NULL表示不保存中间结果，而只保存最终结果到x。
int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol=1.0e-3,
		  int NormControl=0, double MaxStep=-1, double InitialStep=-1, FILE* fid=NULL);

//符号函数
template<class T> inline int sign(const T & x)
{
	if(x>0)
		return 1;
	else if(x<0)
		return -1;
	else
		return 0;
}

//求最大值
template <class T>
inline T max(T x, T y) 
{ 
	return (x>y)?x:y;
}


//求最小值
template <class T>
inline T min(T x, T y) 
{
	return (x<y)?x:y;
}

#endif