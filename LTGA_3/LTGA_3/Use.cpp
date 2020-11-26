#include"Use.h"

//�����û������ָ�꺯��fnobj������ָ������Ա����ĵ���ͨ������Զ���⣬objgrd���ؼ���
int npsolfnobj(int* mode, int& n, double* x, double* objf, double* objgrd, int& nstate)
{
	bool flag=false;
	flag=fnobj(n, x, *objf);
	if(!flag)	*mode=-1;//���ָ����㲻�ɹ�����mode��ֵ��Ϊ������npsol_����ʱ�ᰴ�쳣����
	// objgrd[0]=-400.0*(x[1]-x[0]*x[0])*x[0]-2.0*(1.0-x[0]);
	// objgrd[1]=200.0*(x[1]-x[0]*x[0]);	
	return 0;//npsol_����ʱҪ�󷵻�����
}

//�����û������ָ�꺯��fncon�����ڷ�����Լ�������Ա����ĵ���ͨ������Զ���⣬cJac���ؼ���
int npsolfncon(int* mode, int &ncnln, int& n, int& ldJ, int* needc, double* x, double* nlc, double* cJac, int& nstate)
{
	bool flag=false;
	flag=fncon(n, ncnln, x, nlc);
	if(!flag)	*mode=-1;//���ָ����㲻�ɹ�����mode��ֵ��Ϊ������npsol_����ʱ�ᰴ�쳣����	
	// cJac[0]=2.0*(x[0]-1.0/3.0);
	// cJac[1]=x[1];
	// cJac[2]=2.0*(x[1]-1.0/3.0);
	// cJac[3]=x[0];
	return 0;//npsol_����ʱҪ�󷵻�����
}

int NPSol(double* x, double& objf, int n, int nclin, int ncnln, double** AM, const double* blow, const double* bup)
{
	int nbnd, ldA, ldJ, ldR, leniw, lenw;
	int flag, i, j;
	int inform=0, iter=0;
	nbnd=n+nclin+ncnln;//Լ���ϡ��½������ά��
	ldA=nclin;//����Լ��������
	ldJ=ncnln;//������Լ��������
	ldR=n;//Hessian������ص�ͳ����
	leniw=3*n+nclin+2*ncnln;//�漰�м������Ҫ�����������ά��
	lenw=2*n*n+n*nclin+2*n*ncnln+20*n+11*nclin+21*ncnln;//�漰�м������Ҫ��˫���������ά��

	int idA, idbl, idbu, idnlc, idcJac, idclamda, idobjgrd, idR, idwork;
	int iistate, iiwork;
	//���������±����ʼָ��
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

	double* dworkarr=new double[idwork+lenw];//�漰�м������Ҫ��˫��������
	int* iworkarr=new int[iiwork+leniw];//�漰�м������Ҫ����������

	doublereal* A=(doublereal*)(&dworkarr[idA]);//����Լ������ϵ�����󣬰�һά��ţ����з���
	doublereal* bl=(doublereal*)(&dworkarr[idbl]);//Լ���½�����
	doublereal* bu=(doublereal*)(&dworkarr[idbu]);//Լ���Ͻ�����
	doublereal* nlc=(doublereal*)(&dworkarr[idnlc]);//������Լ������
	doublereal* cJac=(doublereal*)(&dworkarr[idcJac]);//������Լ���ſɱȾ��󣬰�һά��ţ����з���
	doublereal* clamda=(doublereal*)(&dworkarr[idclamda]);//Լ����Ӧ�����ϳ�������
	doublereal* objgrd=(doublereal*)(&dworkarr[idobjgrd]);//ָ����Ա������ݶ�����
	doublereal* R=(doublereal*)(&dworkarr[idR]);//Hessian������صľ��󣬰�һά��ţ����з���
	doublereal* work=(doublereal*)(&dworkarr[idwork]);//˫���ȹ�������
	integer* istate=(integer*)(&iworkarr[iistate]);//��ʶԼ��״̬���������飬��ϸ˵����NPSOL5-0 Manual.pdf��3��
	integer* iwork=(integer*)(&iworkarr[iiwork]);//���ι�������
	//�������2ά����Լ��ϵ������AM���з�����һά����A��
	for(i=0;i<nclin;i++)
		for(j=0;j<n;j++)
			A[i+j*nclin]=AM[i][j];
	//�������Լ���ϡ��½�����bup��blow�����bu��bl��
	for(i=0;i<nbnd;i++)
	{
		bl[i]=blow[i];
		bu[i]=bup[i];
	}
	//��Ƶ�������������Ҫ��Derivative level = 0��ʾָ��ͷ�����Լ�����ſɱȾ��󶼷ǽ�����ͨ����ֵ�����⣨�ɲ��֣���
	//1��ʾ��������Լ�����ſɱȾ���ͨ����ֵ�����⣨�ɲ��֣���2��ʾ��ָ����ſɱȾ���ͨ����ֵ�����⣨�ɲ��֣���3��ʾȫ����������
	//����20��ʾDerivative level = 0�����ո����ڹ�20���ַ�
	flag=npoptn_("Derivative level = 0",20);
	////���ѡ����NPSOL5-0 Manual.pdf��9��
	flag=npoptn_("Print level = 5",15);
	flag=npoptn_("Print file = 10",15);
	//����ԭ�����еķ����Թ滮������
	flag=npsol_((integer*)(&n), (integer*)(&nclin), (integer*)(&ncnln), (integer*)(&ldA), (integer*)(&ldJ), (integer*)(&ldR),
		A, bl, bu, (U_fp)npsolfncon, (U_fp)npsolfnobj, (integer*)(&inform), (integer*)(&iter), istate, nlc, cJac, clamda,
		(doublereal*)(&objf), objgrd, R, (doublereal*)(x), iwork, (integer*)(&leniw), work, (integer*)(&lenw));

	delete[] dworkarr;
	delete[] iworkarr;
	return inform;
}