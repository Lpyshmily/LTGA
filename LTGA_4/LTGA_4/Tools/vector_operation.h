#ifndef _VECTOROPERATION_
#define _VECTOROPERATION_

#include <math.h>
#include<assert.h>
#include <stdio.h>
using namespace std;



//*************************释放内存****************************
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B=NULL;
}

//***************************赋值******************************
//将向量A的值赋给B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}
//将三个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}
//将六个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//**************************基本运算*******************************
//B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}
//向量C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}
//向量C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}
//向量C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}
//向量C[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B;
}
//向量C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}
//向量C[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B*A[I_];
}
//向量C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}
//向量C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//**************************向量运算********************************
//求内积
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}
//求外积C=AXB,不能用V_Cross(B,B,A)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}
//求A[i]=|B[i]|
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			A[I_]=B[I_];
		else
			A[I_]=-B[I_];
	}
}
//求1-范数
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			result+=B[I_];
		else
			result-=B[I_];
	}
	return result;
}
//求2-范数
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}
//求无穷-范数
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}
//求单位向量
template<class T> inline void V_Vers(T* C, const T* B, int N)
{
	double a = V_Norm2(B, N);
	V_Divid(C,B,a,N);
}

//***************************其他运算********************************
//求最大元素
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]>result) result=B[I_];
	return result;
}
//求最大元素
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}
//求最小元素
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}
//求最小元素
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}
//向量B与A每个元素都相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++)	if(B[I_]!=A[I_])return false;
	return true;
}

//将九个值依次赋给B
template<class T> inline void M_Copy(T* B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	T temp[9]={a11, a12, a13, a21, a22, a23, a31, a32, a33};
	for(int I_=0;I_<9;I_++)
		B[I_]=temp[I_];
}

template<class T> inline void M_Copy(T* B, T angle, int axis)
{
	assert(axis==1||axis==2||axis==3);
	for(int I_=0;I_<9;I_++) B[I_]=0;
	if(axis==1)
	{
		B[0]=1.0;
		B[4]=B[8]=cos(angle);
		B[5]=sin(angle);
		B[7]=-B[5];
	}
	if(axis==2)
	{
		B[4]=1.0;
		B[0]=B[8]=cos(angle);
		B[6]=sin(angle);
		B[2]=-B[6];
	}
	if(axis==3)
	{
		B[8]=1.0;
		B[0]=B[4]=cos(angle);
		B[1]=sin(angle);
		B[3]=-B[1];
	}
}

template<class T> inline void M_Multi(T* C, const T* A, const T* B, int N, int M, int K)
{
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		C[s]=0;
		for(int i=0;i<K;i++) C[s]+=A[I_*K+i]*B[i*M+J_];
	}
}

//取最近接x的整数,返回double型.
inline double ANINT(double x)
{
	double left=fmod(x, 1.0);
	if(fabs(left)<0.5)
		return x-left;
	else if(left>=0.5)
		return x-left+1;
	else
		return x-left-1;
}
//取最近接x的整数,返回int型.
inline int NINT(double x)
{
	double left=ANINT(x);
	return (int)left;
}

inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

inline double acosh(double x)
{	
	assert(x>=1.0);
	return log(x+sqrt(x*x-1.0));
}


#endif