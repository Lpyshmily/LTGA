#include "lambert.h"

double x2tof(double x, double s, double c, int lw, int N)
{
	double t=0.0;
	double am=s/2.0;
	double a=am/(1.0-x*x);
	double beta=0.0, alfa=0.0;
	double temp;
	if(x<1.0) //ELLISSE
	{
		temp=(s-c)/2.0/a;
		if(temp<0.0) temp=0.0;
		beta=2.0*asin(sqrt(temp));
		if(lw)
			beta=-beta;
		alfa=2.0*acos(x);
	}
	else   //IPERBOLE
	{
		alfa=2.0*log(x+sqrt(x*x-1.0));//acosh(x);
		temp=(s-c)/(-2.0*a);
		if(temp<0.0) temp=0.0;
		temp=sqrt(temp);
		beta=2.0*log(temp+sqrt(temp*temp+1.0));//asinh(sqrt((s-c)/(-2.0*a)));
		if(lw)
			beta=-beta;		
	}
	if(a>0.0)
		t=a*sqrt(a)*((alfa-sin(alfa))-(beta-sin(beta))+N*M_2PI);
	else
		t=-a*sqrt(-a)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
	return t;
}

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2, 
				double tf, const double* unith, int& flag, double mu, int way, int N, int branch, int Maxiter, double tol)
{
	//Originally programmed with Matlab source code by:  Dario Izzo (Advanced Concept Team) Date:
	//28/05/2004 Revision:              1 Tested by:             ----------
	//Translated into C++ source code by Jiang FH
	//Inputs:
	//           R1=Position array at departure 
	//           R2=Position array at arrival
	//           tf=Transfer time (scalar)
	//           mu=gravitational parameter (scalar, units have to be
	//           consistent with r1,t units), default 132712439940.0 km^3/s^2 
	//           lw=0 for counterclockwise transfer, 1 for clockwise transfer, default 0     
	//           N=number of revolutions, default 0
	//           branch=0 if the left branch is selected in a problem where N
	//           is not 0 (multirevolution), 1 for the right one, default 0 (好像是lw=0时，branch=1较好)
	//           Maxiter=the maximum iteration, default 60
	//           tol=tolerance, default 1.0e-11
	//			 unith,R1和R2共线时，需提供的轨道法向单位矢量
	//
	//Outputs:
	//           v1=Velocity at departure        (consistent units)
	//           v2=Velocity at arrival
	//           a=semi major axis of the solution
	//           e=eccentricity of the solution 
	flag=1;//-1, negative time;0, not to be converged; //2, R1 and R2 collinear
	int i;
	if(tf<=0.0)
	{
//		cout<<"Negative time as input"<<endl;
		flag=-1;
		return;
	}

	//double tol=1.0e-11;  //Increasing the tolerance does not bring any advantage as the 
	//precision is usually greater anyway (due to the rectification of the tof
	//graph) except near particular cases such as parabolas in which cases a
	//lower precision allow for usual convergence.


	double r1[3], r2[3];
	for(i=0;i<3;i++)
	{
		r1[i]=R1[i];
		r2[i]=R2[i];
	}
	//Non dimensional units
	double R=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
	double V=sqrt(mu/R);
	double T=R/V;

	//working with non-dimensional radii and time-of-flight
	for(i=0;i<3;i++)
	{
		r1[i]/=R;
		r2[i]/=R;
	}
	tf/=T;                     

	//Evaluation of the relevant geometry parameters in non dimensional units
	double r2mod=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);//ArrayNorm2(r2, 3);
	double theta=(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])/r2mod;
	if(theta>=1.0)
		theta=0.0;
	else if(theta<=-1.0)
		theta=M_PI;
	else
		theta=acos(theta); //the real command is useful when theta is very 
										//close to pi and the acos function could return complex numbers
	double crsprod[3];
	crsprod[0]=r1[1]*r2[2]-r1[2]*r2[1];
	crsprod[1]=r1[2]*r2[0]-r1[0]*r2[2];
	crsprod[2]=r1[0]*r2[1]-r1[1]*r2[0];	
	if(crsprod[2]<0.0)
		theta=M_2PI-theta;
	if(way==1)
		theta=M_2PI-theta;
	int lw=0;
	if(theta>=M_PI)
		lw=1;
	
	double c=sqrt(1.0+r2mod*r2mod-2.0*r2mod*cos(theta)); //non dimensional chord
	double s=(1.0+r2mod+c)/2.0;                      //non dimensional semi-perimeter
	double am=s/2.0;                               //minimum energy ellipse semi major axis
	double lambda=sqrt(r2mod)*cos(theta/2.0)/s;    //lambda parameter defined in BATTIN's book


	if(c<=0.0)
	{
		a=tf/(2.0*(1.0+N)*M_PI);
		a=pow(a*a*mu, 1.0/3.0);
		a*=R;
		e=0.0;
		v1[0]=unith[1]*R1[2]-unith[2]*R1[1];
		v1[1]=unith[2]*R1[0]-unith[0]*R1[2];
		v1[2]=unith[0]*R1[1]-unith[1]*R1[0];
		v2[0]=v1[0];
		v2[1]=v1[1];
		v2[2]=v1[2];
		return;
	}

	//We start finding the log(x+1) value of the solution conic:
	////NO MULTI REV --> (1 SOL)
	double x=0.0, inn1=0.0, inn2=0.0, x1=0.0, x2=0.0, y1=0.0, y2=0.0;
	double xnew=0.0, ynew=0.0;
	int iter=0;
	if(N==0)
	{
		inn1=-0.5233;    //first guess point
		inn2=0.5233;     //second guess point
		x1=log(1.0+inn1);
		x2=log(1.0+inn2);
		y1=log(x2tof(inn1,s,c,lw,N))-log(tf);
		y2=log(x2tof(inn2,s,c,lw,N))-log(tf);
    
		//Newton iterations
		double err=1.0;
		int i=0;
		while ((err>tol)&&y1!=y2)//(fabs(y1-y2)>1.0e-14))
		{
			i++;
			xnew=(x1*y2-y1*x2)/(y2-y1);
			ynew=log(x2tof(exp(xnew)-1.0,s,c,lw,N))-log(tf);
			x1=x2;
			y1=y2;
			x2=xnew;
			y2=ynew;
			err=fabs(x1-xnew);
		}
		iter=i;
		x=exp(xnew)-1.0;
	}    
    
		////MULTI REV --> (2 SOL) SEPARATING RIGHT AND LEFT BRANCH
	else
	{
		if (branch==0)
		{
			inn1=-0.5234;
			inn2=-0.2234;
		}
		else
		{
			inn1=0.7234;
			inn2=0.5234;
		}
		x1=tan(inn1*M_PI/2.0);
		x2=tan(inn2*M_PI/2.0);
		y1=x2tof(inn1,s,c,lw,N)-tf;    
		y2=x2tof(inn2,s,c,lw,N)-tf;
		double err=1.0;
		int i=0;
    
		//Newton Iteration
		while ((err>tol) && (i<Maxiter) && y1!=y2)//(fabs(y1-y2)>1.0e-14))
		{
			i++;
			xnew=(x1*y2-y1*x2)/(y2-y1);
			ynew=x2tof(atan(xnew)*2.0/M_PI,s,c,lw,N)-tf;
			x1=x2;
			y1=y2;
			x2=xnew;
			y2=ynew;
			err=fabs(x1-xnew);	    
		}
		x=atan(xnew)*2.0/M_PI;
		iter=i;
	}
	if(iter>=Maxiter)
	{
//		cout<<"Solution does not seem to be converging"<<endl;
		flag=0;
		return;
	}
	//The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
	//now need the conic. As for transfer angles near to pi the lagrange
	//coefficient technique goes singular (dg approaches a zero/zero that is
	//numerically bad) we here use a different technique for those cases. When
	//the transfer angle is exactly equal to pi, then the ih unit vector is not
	//determined. The remaining equations, though, are still valid.
	a=am/(1.0-x*x);
	double beta=0.0, alfa=0.0, psi=0.0, eta=0.0, eta2=0.0;//solution semimajor axis
	double temp;
	//calcolo psi
	if(x<1.0) //ellisse
	{
		temp=(s-c)/2.0/a;
		if(temp<0.0) temp=0.0;
		beta=2.0*asin(sqrt(temp));
		if(lw)
			beta=-beta;		
		alfa=2.0*acos(x);
		psi=(alfa-beta)/2.0;
		eta2=2.0*a*sin(psi)*sin(psi)/s;
		eta=sqrt(eta2);
	}
	else //iperbole
	{
		temp=(c-s)/2.0/a;
		if(temp<0.0) temp=0.0;
		temp=sqrt(temp);
		beta=2.0*log(temp+sqrt(temp*temp+1.0));//asinh(sqrt((c-s)/2.0/a));
		if(lw)
			beta=-beta;		
		alfa=2.0*log(x+sqrt(x*x-1.0));//	alfa=2.0*acosh(x);
		psi=(alfa-beta)/2.0;
		eta2=-2.0*a*sinh(psi)*sinh(psi)/s;
		eta=sqrt(eta2);
	}
//	if(eta2==0) {a=a*R;
	double p=r2mod/am/eta2*sin(theta/2)*sin(theta/2);     //parameter of the solution
	double sigma1=1.0/eta/sqrt(am)*(2.0*lambda*am-(lambda+x*eta));
	e=1.0-p/a;
	e=(e>=0.0)?e:0.0;
	e=sqrt(e);

	temp=sqrt(crsprod[0]*crsprod[0]+crsprod[1]*crsprod[1]+crsprod[2]*crsprod[2]);
//	if(temp<=0.0) {a=a*R;cout<<"Cannot determine velocity vector"<<endl; flag=2;return;}
//	for(i=0;i<3;i++)
//		crsprod[i]/=temp;
	if(temp<=0.0)
	{
//		flag=2;
		for(i=0;i<3;i++) crsprod[i]=unith[i];
	}
	else
		for(i=0;i<3;i++) crsprod[i]/=temp;
		
	
	if(lw)
		for(i=0;i<3;i++) crsprod[i]=-crsprod[i];
	double vr1 = sigma1;
	double vt1 = sqrt(p);
	double v1tan[3]={0.0}, v2tan[3]={0.0};
	v1tan[0]=crsprod[1]*r1[2]-crsprod[2]*r1[1];
	v1tan[1]=crsprod[2]*r1[0]-crsprod[0]*r1[2];
	v1tan[2]=crsprod[0]*r1[1]-crsprod[1]*r1[0];
	v2tan[0]=crsprod[1]*r2[2]-crsprod[2]*r2[1];
	v2tan[1]=crsprod[2]*r2[0]-crsprod[0]*r2[2];
	v2tan[2]=crsprod[0]*r2[1]-crsprod[1]*r2[0];
	for(i=0;i<3;i++) v2tan[i]/=r2mod;

	double vt2=vt1/r2mod;
	double vr2=-vr1+(vt1-vt2)/tan(theta/2.0);
	for(i=0;i<3;i++)
	{
		v1[i]=V*(vr1*r1[i]+vt1*v1tan[i]);
		v2[i]=V*(vr2/r2mod*r2[i]+vt2*v2tan[i]);
	}
	a=a*R;
	return;
}
// 圈数从0到最大可能数穷举，保留特征速度最小的结果
// 输入变量
// rv0 rv1 始端和终端的位置速度
// t 转移时间
// GM 引力系数
// 输出变量
// flag 求解成功为1，失败为0
// dv0 始端速度增量 大小为Mdv0
// dv1 末端速度增量 大小为Mdv1
// N 速度增量最小的圈数
// branch 速度增量最小时的分支 0或1
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
				const double* rv0, const double* rv1, double t, double GM)
{
	double vt0[6],am,vt1[3],v0[3],v1[3],a,e,unith[3],dtemp0,dtemp1,dvmin=1.0e10,dv;
	int i,Nmax,n,flag0;
	flag=0;
	V_Cross(unith,rv0,&rv0[3]);
	dtemp0=V_Norm2(unith,3);
	if(dtemp0>0.0)	for(i=0;i<3;i++) unith[i]/=dtemp0;
	else {unith[0]=0.0;unith[1]=0.0;unith[2]=1.0;}
	for(i=0;i<3;i++) vt0[i]=rv0[i]-rv1[i];
	am=(V_Norm2(rv0,3)+V_Norm2(rv1,3)+V_Norm2(vt0,3))/4.0;//最小椭圆转移轨道半长轴
	Nmax=(int)floor(t/(M_2PI*sqrt(am*am*am/GM)));//最多可能的转移圈数
	for(n=0;n<=Nmax;n++)//圈数从0到最大可能数穷举，保留特征速度最小的结果
	{
		for(int j=0;j<2;j++)//对于多圈问题，有左右分枝两个解
		{
			lambert(v0, v1, a, e, rv0, rv1, t, unith, flag0, GM, 0, n, j);
			if(flag0!=1||e>=1.0) continue;//排除双曲和抛物轨道
			for(i=0;i<3;i++) 
			{
				vt0[i]=v0[i]-rv0[3+i];
				vt1[i]=rv1[3+i]-v1[i];
			}			
			dtemp0=V_Norm2(vt0,3); // 初始时刻速度增量
			dtemp1=V_Norm2(vt1,3); // 末端速度增量
			dv=dtemp0+dtemp1; // 两个速度增量相加
			if(dv<dvmin)
			{
				flag=1;
				N=n;
				branch=j;
				dvmin=dv;
				V_Copy(dv0,vt0,3);
				V_Copy(dv1,vt1,3);
				Mdv0=dtemp0;
				Mdv1=dtemp1;
			}
			if(n==0) break;
		}
	}
}

lambert_solver::lambert_solver(const double* rv1,const double* rv2,const double tof,const double mu,const int way)
{
	V_Copy(m_rv1,rv1,6);
	V_Copy(m_rv2,rv2,6);
	m_tof=tof;
	m_mu=mu;
	m_way=way;
	m_flag=0;
}
// 求解给定圈数和分支时的速度增量
void lambert_solver::solve_single(const int N,const int branch)
{
	double v1[3],v2[3],a,e,unith[3];
	V_Cross(unith,m_rv1,&m_rv1[3]);
	if(V_Norm2(unith,3)>0.0)
		V_Vers(unith,unith,3);
	else {unith[0]=0.0;unith[1]=0.0;unith[2]=1.0;}
	lambert(v1, v2, a, e, m_rv1, m_rv2, m_tof, unith, m_flag, m_mu, m_way, N, branch);
	V_Minus(m_dv1,v1,&m_rv1[3],3);
	V_Minus(m_dv2,&m_rv2[3],v2,3);
	m_Mdv1=V_Norm2(m_dv1,3);
	m_Mdv2=V_Norm2(m_dv2,3);

	V_Copy(m_original_rvt, m_rv2, 3);
	V_Copy(&m_original_rvt[3], v2, 3);
}
// 不指定圈数，调用LamEval函数寻找速度增量最小的圈数，并获取速度增量
int lambert_solver::solve_multi()
{
	LambEval(m_flag, m_dv1, m_dv2, m_Mdv1, m_Mdv2, m_N, m_branch, m_rv1, m_rv2, m_tof, m_mu);
	return m_flag;
}

void lambert_solver::get_dv1(double* dv1)
{
	V_Copy(dv1,m_dv1,3);
}

void lambert_solver::get_dv2(double* dv2)
{
	V_Copy(dv2,m_dv2,3);
}

double lambert_solver::get_Mdv1()
{
	return m_Mdv1;
}

double lambert_solver::get_Mdv2()
{
	return m_Mdv2;
}
// 获取未经改变的末端位置速度
void lambert_solver::get_originalRVT(double* rv)
{
	V_Copy(rv, m_original_rvt, 6);
}
