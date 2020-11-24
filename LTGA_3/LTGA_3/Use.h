#ifndef USE_H
#define USE_H
#include"cnpsol.h"

//NPSol调用相关的3个函数，函数npsolfnobj和npsolfncon的输入输出介绍见NPSOL5-0 Manual.pdf第4章
int npsolfnobj(int* mode, int& n, double* x, double* objf, double* objgrd, int& nstate);
int npsolfncon(int* mode, int &ncnln, int& n, int& ldJ, int* needc, double* x, double* nlc, double* cJac, int& nstate);

//另写的调用核心函数npsol_的函数
//输入:
//	x:初值
//	n:x的维数
//	nclin:线性约束的维数
//	ncnln:非线性约束的维数
//	AM:线性约束的矩阵系数,nclin行,n列
//	blow:约束下界数组,n+nclin+ncnln维
//	bup:约束上界数组,n+nclin+ncnln维
//输出:
//	x:极小值对应的解
//	objf:指标极小值
//	函数返回整数，其值为0时表示结果基本合理，其它值表示的意义见NPSOL5-0 Manual.pdf第7页
int NPSol(double* x, double& objf, int n, int nclin, int ncnln, double** AM, const double* blow, const double* bup);


//指标函数fnobj和非线性约束函数fncon
//输入:
//	n:x的维数
//	ncnln:非线性约束维数
//	x:自变量数组
//输出:
//	objf:指标值
//	nlc:非线性约束值（数组）
//	如果对指标和约束计算成功，返回true，否则应返回false
bool fnobj(int n, const double* x, double& objf);
bool fncon(int n, int ncnln, const double* x, double* nlc);

#endif