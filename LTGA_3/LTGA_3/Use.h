#ifndef USE_H
#define USE_H
#include"cnpsol.h"

//NPSol������ص�3������������npsolfnobj��npsolfncon������������ܼ�NPSOL5-0 Manual.pdf��4��
int npsolfnobj(int* mode, int& n, double* x, double* objf, double* objgrd, int& nstate);
int npsolfncon(int* mode, int &ncnln, int& n, int& ldJ, int* needc, double* x, double* nlc, double* cJac, int& nstate);

//��д�ĵ��ú��ĺ���npsol_�ĺ���
//����:
//	x:��ֵ
//	n:x��ά��
//	nclin:����Լ����ά��
//	ncnln:������Լ����ά��
//	AM:����Լ���ľ���ϵ��,nclin��,n��
//	blow:Լ���½�����,n+nclin+ncnlnά
//	bup:Լ���Ͻ�����,n+nclin+ncnlnά
//���:
//	x:��Сֵ��Ӧ�Ľ�
//	objf:ָ�꼫Сֵ
//	����������������ֵΪ0ʱ��ʾ���������������ֵ��ʾ�������NPSOL5-0 Manual.pdf��7ҳ
int NPSol(double* x, double& objf, int n, int nclin, int ncnln, double** AM, const double* blow, const double* bup);


//ָ�꺯��fnobj�ͷ�����Լ������fncon
//����:
//	n:x��ά��
//	ncnln:������Լ��ά��
//	x:�Ա�������
//���:
//	objf:ָ��ֵ
//	nlc:������Լ��ֵ�����飩
//	�����ָ���Լ������ɹ�������true������Ӧ����false
bool fnobj(int n, const double* x, double& objf);
bool fncon(int n, int ncnln, const double* x, double* nlc);

#endif