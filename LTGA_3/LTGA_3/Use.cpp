#include"Use.h"

//调用用户定义的指标函数fnobj，由于指标关于自变量的导数通过差分自动求解，objgrd不必计算
int npsolfnobj(int* mode, int& n, double* x, double* objf, double* objgrd, int& nstate)
{
	bool flag=false;
	flag=fnobj(n, x, *objf);
	if(!flag)	*mode=-1;//如果指标计算不成功，将mode的值置为负数，npsol_调用时会按异常处理
	// objgrd[0]=-400.0*(x[1]-x[0]*x[0])*x[0]-2.0*(1.0-x[0]);
	// objgrd[1]=200.0*(x[1]-x[0]*x[0]);	
	return 0;//npsol_调用时要求返回整数
}

//调用用户定义的指标函数fncon，由于非线性约束关于自变量的导数通过差分自动求解，cJac不必计算
int npsolfncon(int* mode, int &ncnln, int& n, int& ldJ, int* needc, double* x, double* nlc, double* cJac, int& nstate)
{
	bool flag=false;
	flag=fncon(n, ncnln, x, nlc);
	if(!flag)	*mode=-1;//如果指标计算不成功，将mode的值置为负数，npsol_调用时会按异常处理	
	// cJac[0]=2.0*(x[0]-1.0/3.0);
	// cJac[1]=x[1];
	// cJac[2]=2.0*(x[1]-1.0/3.0);
	// cJac[3]=x[0];
	return 0;//npsol_调用时要求返回整数
}

int NPSol(double* x, double& objf, int n, int nclin, int ncnln, double** AM, const double* blow, const double* bup)
{
	int nbnd, ldA, ldJ, ldR, leniw, lenw;
	int flag, i, j;
	int inform=0, iter=0;
	nbnd=n+nclin+ncnln;//约束上、下界数组的维数
	ldA=nclin;//线性约束的行数
	ldJ=ncnln;//非线性约束的行数
	ldR=n;//Hessian矩阵相关的统领数
	leniw=3*n+nclin+2*ncnln;//涉及中间计算需要的整形数组的维数
	lenw=2*n*n+n*nclin+2*n*ncnln+20*n+11*nclin+21*ncnln;//涉及中间计算需要的双精度数组的维数

	int idA, idbl, idbu, idnlc, idcJac, idclamda, idobjgrd, idR, idwork;
	int iistate, iiwork;
	//各种数组下标的起始指标
	idA=0;
	idbl=idA+nclin*n;
	idbu=idbl+nbnd;
	idnlc=idbu+nbnd;
	idcJac=idnlc+ncnln;
	idclamda=idcJac+ncnln*n;
	idobjgrd=idclamda+nbnd;
	idR=idobjgrd+n;
	idwork=idR+n*n;
	
	iistate=0;
	iiwork=iistate+nbnd;

	double* dworkarr=new double[idwork+lenw];//涉及中间计算需要的双精度数组
	int* iworkarr=new int[iiwork+leniw];//涉及中间计算需要的整形数组

	doublereal* A=(doublereal*)(&dworkarr[idA]);//线性约束矩阵系数矩阵，按一维存放，逐列放置
	doublereal* bl=(doublereal*)(&dworkarr[idbl]);//约束下界数组
	doublereal* bu=(doublereal*)(&dworkarr[idbu]);//约束上界数组
	doublereal* nlc=(doublereal*)(&dworkarr[idnlc]);//非线性约束数组
	doublereal* cJac=(doublereal*)(&dworkarr[idcJac]);//非线性约束雅可比矩阵，按一维存放，逐列放置
	doublereal* clamda=(doublereal*)(&dworkarr[idclamda]);//约束对应的拉氏乘子数组
	doublereal* objgrd=(doublereal*)(&dworkarr[idobjgrd]);//指标对自变量的梯度数组
	doublereal* R=(doublereal*)(&dworkarr[idR]);//Hessian矩阵相关的矩阵，按一维存放，逐列放置
	doublereal* work=(doublereal*)(&dworkarr[idwork]);//双精度工作数组
	integer* istate=(integer*)(&iworkarr[iistate]);//标识约束状态的整形数组，详细说明见NPSOL5-0 Manual.pdf第3章
	integer* iwork=(integer*)(&iworkarr[iiwork]);//整形工作数组
	//将输入的2维线性约束系数矩阵AM逐列放置于一维数组A中
	for(i=0;i<nclin;i++)
		for(j=0;j<n;j++)
			A[i+j*nclin]=AM[i][j];
	//将输入的约束上、下界数组bup、blow存放于bu、bl中
	for(i=0;i<nbnd;i++)
	{
		bl[i]=blow[i];
		bu[i]=bup[i];
	}
	//设计导数级数，很重要，Derivative level = 0表示指标和非线性约束的雅可比矩阵都非解析，通过数值差分求解（可部分），
	//1表示仅非线性约束的雅可比矩阵通过数值差分求解（可部分），2表示仅指标的雅可比矩阵通过数值差分求解（可部分），3表示全部解析给出
	//数字20表示Derivative level = 0包括空格在内共20个字符
	flag=npoptn_("Derivative level = 0",20);
	////输出选项，详见NPSOL5-0 Manual.pdf第9章
	flag=npoptn_("Print level = 5",15);
	flag=npoptn_("Print file = 10",15);
	//调用原程序中的非线性规划主函数
	flag=npsol_((integer*)(&n), (integer*)(&nclin), (integer*)(&ncnln), (integer*)(&ldA), (integer*)(&ldJ), (integer*)(&ldR),
		A, bl, bu, (U_fp)npsolfncon, (U_fp)npsolfnobj, (integer*)(&inform), (integer*)(&iter), istate, nlc, cJac, clamda,
		(doublereal*)(&objf), objgrd, R, (doublereal*)(x), iwork, (integer*)(&leniw), work, (integer*)(&lenw));

	delete[] dworkarr;
	delete[] iworkarr;
	return inform;
}