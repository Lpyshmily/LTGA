#include "base_functions.h"
#include "vector_operation.h"
#include "MinPack\cminpack.h"

// 这个函数没有在任何头文件中声明过，仅在此文件中使用
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}
double e2f(int& flag, const double E, const double e)
{
	if(e<0.0) {flag=0;return E;}

	double f=0.0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double E0=fmod(E, M_2PI);
		if(E0>M_PI)
			E0-=M_2PI;
		if(E0<-M_PI)
			E0+=M_2PI;
		f=2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*E0));
		f=f+E-E0;
	}
	else if(e>1.0)//双曲线轨道
		f=2.0*atan(sqrt((e+1.0)/(e-1.0))*tanh(0.5*E));
	else//抛物线轨道
	{
		f=E;//抛物线轨道没有定义偏近点角.在此将真近点角的值置为E
	}
	flag=1;
	return f;
}

// E2M 根据偏近点角和偏心率求平近点角
//输入参数E, e:
//   E:偏近点角,单位：弧度,对于双曲轨道,指H
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   epsilon:最小精度,默认为1.0e-12
//输出参数M:
// 	平近点角(弧度),对于双曲轨道,指N
double e2m(int& flag, double E, double e)
{
	if(e<0.0){flag=0;return E;}
	double M=0.0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double E0=fmod(E, M_2PI);
		M=E0-e*sin(E0);
		M=M+E-E0;
	}
	else if(e>1.0)//双曲线轨道
		M=e*sinh(E)-E;   
	else//(abs(e-1.0)<epsilon)抛物线轨道
	{
		M=E;
//		cout<<"抛物线轨道没有定义偏近点角.在此将平近点角值置为E."<<endl;
	}
	flag=1;
	return M;
}

double f2e(int& flag, double f, double e)
{
	if(e<0.0){flag=0;return f;}

	double E=0.0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double f0=fmod(f, M_2PI);
		if(f0>M_PI)
			f0-=M_2PI;
		if(f0<-M_PI)
			f0+=M_2PI;
		E=2.0*atan(sqrt((1.0-e)/(1.0+e))*tan(0.5*f0));
		E+=f-f0;
	}
	else if(e>1.0)//双曲线轨道
	{
		if(f>M_PI-acos(1.0/e)||f<-M_PI+acos(1.0/e))
		{
//			cout<<"不可能达到的双曲轨道."<<endl;
			flag=0;
			return f;
		}
		else
			E=2.0*atanh(sqrt((e-1.0)/(1.0+e))*tan(0.5*f));
	}
	else//if(abs(e-1.0)<epsilon)抛物线轨道
	{
		E=f;
//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为f."<<endl;
	}
	flag=1;
	return E;
}

// M2E 根据平近点角和偏心率求偏近点角
//输入参数M, e, MaxIter, epsilon:
//   M:平近点角,单位：弧度,对于双曲轨道,指N
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   MaxIter:最大迭代次数,默认为60
//   epsilon:迭代绝对误差,默认为1e-14;
//输出参数E:
// 	偏近点角（弧度），对双曲轨道，指H
double m2e(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if(epsilon<=0.0||MaxIter<1||e<0.0){flag=0;return M;}

	//迭代方法见《Solar System Dynamics》Chapter2,Carl D.Murray and Stanley F.Dermott著
	double E=0.0, Minus=0.0, DeMinus=0.0, DeDeMinus=0.0, DeDeDeMinus=0.0, Delta1=0.0, Delta2=0.0, Delta3=0.0;
	int N=0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double RM=fmod(M, M_2PI);
		if(RM<0.0)
			RM+=M_2PI;
		double sinRM=sin(RM);
		E=RM+0.85*e*Sign(sinRM);
		N=0;   
		Delta3=1.0;
		while(fabs(Delta3)>=epsilon&&N<MaxIter)
		{
			Minus=E-e*sin(E)-RM;
			DeMinus=1.0-e*cos(E);
			DeDeMinus=e*sin(E);
			DeDeDeMinus=e*cos(E);
			Delta1=-Minus/DeMinus;
			Delta2=-Minus/(DeMinus+0.5*Delta1*DeDeMinus);
			Delta3=-Minus/(DeMinus+0.5*Delta2*DeDeMinus+1.0/6.0*Delta2*Delta2*DeDeDeMinus);
			E=E+Delta3;
			N=N+1;
		}    
		E=E+M-RM;
	}
	else if(e>1.0)//双曲线轨道
	{
		E=asinh(M/e);
		Delta3=1.0;
		N=0;
		while(fabs(Delta3)>=epsilon&&N<MaxIter)
		{
			Minus=e*sinh(E)-E-M;
			DeMinus=e*cosh(E)-1.0;
			DeDeMinus=e*sinh(E);
			DeDeDeMinus=e*cosh(E);
			Delta1=-Minus/DeMinus;
			Delta2=-Minus/(DeMinus+0.5*Delta1*DeDeMinus);
			Delta3=-Minus/(DeMinus+0.5*Delta2*DeDeMinus+1.0/6.0*Delta2*Delta2*DeDeDeMinus);
			E=E+Delta3;
			N=N+1;
		}  
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		E=M;
//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为M."<<endl;
	}
	if(((e>=0.0&&e<1.0)||(e>1.0))&&fabs(Delta3)>=5.0*epsilon&&N>=MaxIter)
	{
//		cout<<"迭代不收敛,请降低精度epsilon或增加迭代次数限制."<<endl;
		flag=0;
		return M;
	}
	flag=1;
	return E;
}

double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if(mu<=0.0||MaxIter<1||a<=0.0||e<0.0){flag=0;return f0;}

	double ft=0.0;
	if((e>=0.0&&e<1.0)||(e>1.0))//圆,椭圆,双曲轨道
	{
		double E=f2e(flag,f0,e);
		if(flag==0) return f0;
		double M=e2m(flag,E,e);
		if(flag==0) return f0;
		M+=sqrt(mu/(a*a*a))*dt;
		E=m2e(flag,M,e,MaxIter,epsilon);
		if(flag==0) return f0;
		ft=e2f(flag,E,e);
		if(flag==0) return f0;
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		if((f0<-M_PI)||(f0>M_PI))
		{
//			cout<<"对于抛物线轨道，初始真近点角应在-180至180度之间."<<endl;
			flag=0;return f0;
		}
		else if(f0>M_PI||f0<-M_PI)
			ft=f0;
		else
		{
			double B=0.75*sqrt(2.0*mu/(a*a*a))*dt+0.5*tan(0.5*f0)*((tan(0.5*f0))*(tan(0.5*f0))+3.0);
			double B1B=B+sqrt(1.0+B*B);
			double tanv=0.0;
			if(fabs(dt)<M_2PI*sqrt((a*a*a)/mu)/1000.0)//推进时间为小量的情况
			{
				double A=pow(B1B, 2.0/3.0);
				tanv=2.0*A*B/(1.0+(1.0+A)*A);
			}
			else//不是小量的情况
			{
				double temp=pow(B1B, 1.0/3.0);
				tanv=temp-1.0/temp;
			}
			ft=2.0*atan(tanv);
		}
	}
	flag=1;
	return ft;
}

void coe2rv(int& flag, double* rv, const double* coe, const double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<=0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>M_PI)
		return;
	if((coe[1]*cos(coe[5]))<-1.0)
		return;

	double p=coe[0]*fabs(1.0-coe[1]*coe[1]);//半通径
	if(coe[1]==1.0)//如果是抛物线轨道,区别对待.
		p=2.0*coe[0];	

	double sini, cosi, sinO, cosO, sino, coso;
	sini=sin(coe[2]);
	cosi=cos(coe[2]);
	sinO=sin(coe[3]);
	cosO=cos(coe[3]);
	sino=sin(coe[4]);
	coso=cos(coe[4]);

	//轨道平面法向单位矢量,即角动量单位矢量
	double HVector[3]={sini*sinO, -sini*cosO, cosi};
      
	//偏心率单位矢量,或叫Laplace矢量
	double PVector[3]={cosO*coso-sinO*sino*cosi, sinO*coso+cosO*sino*cosi, sino*sini};

	//半通径方向单位矢量,PVector,QVector,HVector构成右手坐标系
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	

	V_Cross(QVector, HVector, PVector);

	double r=0.0;
	if((coe[1]*cos(coe[5]))+1.0<=0.0)
	{
//		cout<<"抛物或双曲轨道到达无穷远处."<<endl;
		r=1.0e308;
	}
	else
		r=p/(1.0+coe[1]*cos(coe[5]));	

	for(int i=0;i<3;i++)
	{
		rv[i]=r*(cos(coe[5])*PVector[i]+sin(coe[5])*QVector[i]);
		rv[3+i]=sqrt(mu/p)*(-sin(coe[5])*PVector[i]+(cos(coe[5])+coe[1])*QVector[i]);
	}
	flag=1;
	return;
}

void rv2coe(int& flag, double* coe, const double* RV, const double mu)
{
	int i;
	flag=0;
	if(mu<=0.0)
		return;

	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R,3);//距离
	double velocity=V_Norm2(V,3);//速度
	if(radius<=0.0||velocity<=0.0)
		return;
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//径向单位矢量    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector,unitR,unitV);
	double h=radius*velocity*V_Norm2(hvector,3);//角动量值
	if(h<=0.0)
		return;
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector,3);//轨道面法向单位矢量
	//偏心率矢量
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];
	coe[1]=V_Norm2(evector,3);//偏心率
	double p=h*h/mu;
	if(coe[1]==1.0)
		coe[0]=0.5*p;//抛物线轨道的近星距
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//半长轴
	bool judge=(coe[1]>0.0);
	double unite[3]={0.0};
	if(judge)
		for(i=0;i<3;i++) unite[i]=evector[i]/coe[1];//偏心率单位矢量
	coe[2]=acos(unith[2]);//轨道倾角

	double unitN[3]={-unith[1], unith[0], 0.0};//节线矢量,未归一化

	double temp[3];

	if(V_Norm2(unitN,3)==0.0)
	{
		coe[3]=0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if(!judge)
		{
			coe[4]=0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5]=atan2(unitR[1]*unith[2],unitR[0]);//真近点角
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4]=atan2(unite[1]*unith[2],unite[0]); //近星点幅角       
			coe[5]=atan2(V_Dot(unith,temp,3), V_Dot(unite,unitR,3));
		}
	}
	else
	{
		V_Cross(temp, unitN, unitR);
		coe[3]=atan2(unith[0],-unith[1]);
		coe[5]=atan2(V_Dot(unith,temp,3), V_Dot(unitN,unitR,3));
		if(!judge)
		{
			coe[4]=0.0;
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4]=atan2(V_Dot(unith,temp,3), V_Dot(unite,unitN,3));
			coe[5]=coe[5]-coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3]=fmod(coe[3], M_2PI);
	if(coe[3]<0.0)
		coe[3]+=M_2PI;
	coe[4]=fmod(coe[4], M_2PI);
	if(coe[4]<0.0)
		coe[4]+=M_2PI;
	coe[5]=fmod(coe[5], M_2PI);
	if(coe[1]>=1.0)
	{
		if(coe[5]>M_PI-acos(1.0/coe[1]))
			coe[5]-=M_2PI;
		else if(coe[5]<-M_PI+acos(1.0/coe[1]))
			coe[5]+=M_2PI;
	}
	flag=1;
	return;
}

void coe2ee(int&flag, double* ee, const double* coe, double mu)
{
flag=0;
if(mu<=0.0||coe[0]<0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>M_PI)	
	return;	
if((coe[1]*cos(coe[5]))<-1.0)
//		cout<<"不可能达到的双曲轨道."<<endl;
	return;
	
ee[0]=coe[0]*fabs(1.0-coe[1]*coe[1]);//半通径p
if(coe[1]==1.0)//如果是抛物线轨道,区别对待.
	ee[0]=2.0*coe[0];
ee[1]=coe[1]*cos(coe[4]+coe[3]);//f
ee[2]=coe[1]*sin(coe[4]+coe[3]);//g
double temp=tan(coe[2]/2.0);
ee[3]=temp*cos(coe[3]);//h
ee[4]=temp*sin(coe[3]);//k
ee[5]=coe[4]+coe[3]+coe[5];
ee[5]=fmod(ee[5],M_2PI);
if(ee[5]<0.0)
	ee[5]+=M_2PI;
flag=1;
return ;
}

void ee2coe(int& flag, double* coe, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<=0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	coe[1]=sqrt(f*f+g*g);
	if(coe[1]==1.0)
		coe[0]=0.5*p;//抛物线轨道的近星距
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//半长轴
	double temp=sqrt(h*h+k*k);
	coe[2]=2.0*atan(temp);
	if(temp<=0.0)
	{
		coe[3]=0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5]=L;//真近点角
		}
		else
		{
			coe[4]=atan2(g,f); //近星点幅角       
			coe[5]=L-coe[4];
		}
	}
	else
	{
		coe[3]=atan2(k,h);
		coe[5]=L-coe[3];
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			coe[4]=atan2(g*h-f*k,f*h+g*k);
			coe[5]=coe[5]-coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3]=fmod(coe[3], M_2PI);
	if(coe[3]<0.0)
		coe[3]+=M_2PI;
	coe[4]=fmod(coe[4], M_2PI);
	if(coe[4]<0.0)
		coe[4]+=M_2PI;
	coe[5]=fmod(coe[5], M_2PI);
	if(coe[5]<0.0)
		coe[5]+=M_2PI;
	if(coe[1]>=1.0)
	{
		if(coe[5]>M_PI-acos(1.0/coe[1]))
			coe[5]-=M_2PI;
		else if(coe[5]<-M_PI+acos(1.0/coe[1]))
			coe[5]+=M_2PI;
	}
	flag=1;
	return;
}

void ee2rv(int&flag, double* rv, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<=0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	double h_=sqrt(p/mu);
	double n=h_/(1.0+f*cos(L)+g*sin(L));
	double s2=1.0+h*h+k*k;
	double hh_kk=h*h-k*k;
	double hk2=2.0*h*k;

	double r=n*h_*mu;
	rv[0]=r/s2*(cos(L)+hh_kk*cos(L)+hk2*sin(L));
	rv[1]=r/s2*(sin(L)-hh_kk*sin(L)+hk2*cos(L));
	rv[2]=2.0*r/s2*(h*sin(L)-k*cos(L));
	rv[3]=-1.0/h_/s2*(sin(L)+hh_kk*sin(L)-hk2*cos(L)+g-hk2*f+hh_kk*g);
	rv[4]=-1.0/h_/s2*(-cos(L)+hh_kk*cos(L)+hk2*sin(L)-f+hk2*g+hh_kk*f);
	rv[5]=2.0/h_/s2*(h*cos(L)+k*sin(L)+f*h+g*k);
	flag=1;
	return;
}

void rv2ee(int&flag, double* ee, const double* RV, double mu)
{
	flag=0;
	if(mu<=0.0)
		return;
	int i;
	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R, 3);//距离
	double velocity=V_Norm2(V, 3);//速度	
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//径向单位矢量    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h=radius*velocity*V_Norm2(hvector, 3);//角动量值
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector, 3);//轨道面法向单位矢量
	
	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//偏心率矢量	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];	
	//半能径方向矢量,模为e
	double qvector[3];
	V_Cross(qvector,unith,evector);
	//定义的一个矢量
	double unitA[3];
	for(i=0;i<3;i++) unitA[i]=h/mu*V[i]-qvector[i];
	//基本关系式
	// evector=e*[cos(omega)*cos(OMEGA)-sin(omega)*sin(OMEGA)*cos(i),
	//            cos(omega)*sin(OMEGA)+sin(omega)*cos(OMEGA)*cos(i),
	//            sin(omega)*sin(i)];
	// qvector=e*[-sin(omega)*cos(OMEGA)-cos(omega)*sin(OMEGA)*cos(i),
	//            -sin(omega)*sin(OMEGA)+cos(omega)*cos(OMEGA)*cos(i),
	//            cos(omega)*sin(i)];
	// unitR=[cos(OMEGA)*cos(omega+f)-cos(i)*sin(OMEGA)*sin(omega+f),
	//        sin(OMEGA)*cos(omega+f)+cos(i)*cos(OMEGA)*sin(omega+f),
	//        sin(i)*sin(omega+f)];
	// unitA=[-cos(OMEGA)*sin(omega+f)-cos(i)*sin(OMEGA)*cos(omega+f),
	//        -sin(OMEGA)*sin(omega+f)+cos(i)*cos(OMEGA)*cos(omega+f),
	//         sin(i)*cos(omega+f)];
	//推导得出
	//evector[0]+qvector[1]=(1+cos(i))*e*cos(omega+OMEGA);
	//evector[1]-qvector[0]=(1+cos(i))*e*sin(omega+OMEGA);
	//unith[0]=(1+cos(i))*tan(i/2)*sin(OMEGA);
	//unith[1]=-(1+cos(i))*tan(i/2)*cos(OMEGA);
	// unitR[0]+unitA[1]=(1+cos(i))*cos(omega+OMEGA+f);
	// unitR[1]-unitA[0]=(1+cos(i))*sin(omega+OMEGA+f);

	ee[0]=h*h/mu;//p

	if(unith[2]+1.0<=0.0)
	{
//		cout<<"轨道倾角接近180度，不适合用天球轨道根数描述. "<<endl;
		return;
	}
	double cosiadd1=1.0+unith[2];
	ee[1]=(evector[0]+qvector[1])/cosiadd1;//f
	ee[2]=(evector[1]-qvector[0])/cosiadd1;//g
	ee[3]=-unith[1]/cosiadd1;//h
	ee[4]=unith[0]/cosiadd1;//k
	ee[5]=atan2(unitR[1]-unitA[0],unitR[0]+unitA[1]);//L
	//转换到[0,2pi)中
	ee[5]=fmod(ee[5], M_2PI);
	if(ee[5]<0.0)
		ee[5]+=M_2PI;
	flag=1;
	return;
}

void coe02coef(int&flag, double* coe1, const double* coe0, double dt, double mu)
{
	V_Copy(coe1,coe0,6);
	coe1[5]=f0dt2ft(flag,coe1[5],dt,coe1[0],coe1[1],mu);
}

void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu)
{
	double coe[6];
	rv2coe(flag,coe,rv0,mu);
	if(flag==0) return;
	coe[5]=f0dt2ft(flag,coe[5],dt,coe[0],coe[1],mu);
	if(flag==0) return;
	coe[5]=fmod(coe[5],M_2PI);
	coe2rv(flag,rv1,coe,mu);
}

void ee02eef(int&flag, double* ee1, const double* ee0, double dt, double mu)
{
	double coe[6];
	ee2coe(flag,coe,ee0,mu);
	if(flag==0) return;
	coe[5]=f0dt2ft(flag,coe[5],dt,coe[0],coe[1],mu);
	if(flag==0) return;
	coe[5]=fmod(coe[5],M_2PI);
	coe2ee(flag,ee1,coe,mu);
	}

	/**************************************************************************************************************************************/
/*********************************************************考虑j2摄动的轨道根数转换*****************************************************/
/**************************************************************************************************************************************/
void coe2soe(const double* coe, double* soe){
	soe[0] = coe[0];
	soe[1] = coe[4]+coe[5];
	soe[2] = coe[2];
	soe[3] = coe[1]*cos(coe[4]);
	soe[4] = coe[1]*sin(coe[4]);
	soe[5] = coe[3];
}

void soe2coe(const double* soe, double* coe){
	coe[0] = soe[0];
	coe[1] = sqrt(soe[3]*soe[3]+soe[4]*soe[4]);
	coe[2] = soe[2];
	coe[3] = soe[5];
	coe[4] = atan2(soe[4],soe[3]);
	coe[5] = soe[1]-coe[4];
	coe[3]=fmod(coe[3], M_2PI);
	coe[4]=fmod(coe[4], M_2PI);
	coe[5]=fmod(coe[5], M_2PI);
}

//平均根数mcoe(a,e,i,Omega,omega,m)转换为瞬时根数ocoe(a,e,i,Omega,omega,f)
void M2O(const double* mcoe, double* ocoe, double j2){

	double m_coe[6],soe_m[6],soe_o[6];
	int flag;

	V_Copy(m_coe, mcoe, 6);
	m_coe[5]=e2f(flag,m2e(flag,mcoe[5],mcoe[1]),mcoe[1]);
	coe2soe(m_coe,soe_m);

	double a     = soe_m[0];
	double theta = soe_m[1];
	double i     = soe_m[2];
	double q1    = soe_m[3];
	double q2    = soe_m[4];
	double Omega = soe_m[5];

	double cosi = cos(i);
	double sini = sin(i);
	double cosi2 = cosi*cosi;
	double sini2 = 1.0-cosi2;
	double cosi4 = cosi2*cosi2;
	double cosi6 = cosi4*cosi2;

	double cost = cos(theta);
	double sint = sin(theta);
	double cos2t = cost*cost-sint*sint;
	double sin2t = 2.0*sint*cost;
	double cos3t = cost*cos2t-sint*sin2t;
	double sin3t = sin2t*cost+cos2t*sint;
	double cos4t = cos2t*cos2t-sin2t*sin2t;
	double sin4t = 2.0*sin2t*cos2t;
	double cos5t = cos2t*cos3t-sin2t*sin3t;
	double sin5t = sin3t*cos2t+cos3t*sin2t;
	double eta2 = 1.0-q1*q1-q2*q2;
	double eta = sqrt(eta2);
	double eta4 = eta2*eta2;
	double eta6 = eta2*eta4;

	double eps1 = sqrt(1.0-eta2);
	double eps2 = q1*cost+q2*sint;
	double eps3 = q1*sint-q2*cost;
	double THETA = 1.0/(1.0-5.0*cosi2);
	double p = a*eta2;
	double r = p/(1.0+eps2);
	double beta = 1.0/(eta2+eta);
	double F = atan2(r*(1.0+beta*q1*q1)*sint-r*beta*q1*q2*cost+a*q2,r*(1.0+beta*q2*q2)*cost-r*beta*q1*q2*sint+a*q1);
	double theta_lam = theta-(F-q1*sin(F)+q2*cos(F));
	if(abs(theta_lam)>M_2PI)
		theta_lam = fmod(theta_lam, M_2PI);
	while(abs(theta_lam)>M_PI)
		theta_lam=theta_lam-theta_lam/abs(theta_lam)*M_2PI;
	//while(abs(theta_lam)>DPI)
	//	theta_lam=theta_lam-theta_lam/abs(theta_lam)*2.0*DPI;

	double alp = 0.0;
	double C1 = sini/(8.0*a*a*eta2)*(1.0-10.0*THETA*cosi2);
	double C2 = (q1*q2/(16.0*a*a*eta4))*(3.0-55.0*cosi2-280.0*THETA*cosi4-400.0*THETA*THETA*cosi6);
	double lamlp = q1*q2*sini/(1.0+eta)*C1+C2;
	double thetalp = lamlp-0.5*sini/eta2*C1*(q1*q2*(3.0+2.0*eta2/(1.0+eta))+2.0*(q1*sint+q2*cost)+0.5*eps1*sin2t);
	double ilp = 0.5*cosi/eta2*C1*(q1*q1-q2*q2);
	double q1lp = -0.5*q1*sini*C1-q2*C2;
	double q2lp = 0.5*q2*sini*C1+q1*C2;
	double Omegalp = q1*q2*cosi/(8.0*a*a*eta4)*(11.0+80.0*THETA*cosi2+200.0*THETA*THETA*cosi4);

	C1 = 0.25*(1.0-3.0*cosi2)/(a*a*eta4);
	C2 = (theta_lam+eps3)/(a*a*eta4);
	double asp1 = 2.0*a/eta2*C1*((1.0+eps2)*(1.0+eps2)*(1.0+eps2)-eta2*eta);
	double lamsp1 = eps3/(1.0+eta)*C1*((1.0+eps2)*(2.0+eps2)+eta2)+0.75*(1.0-5.0*cosi2)*C2;
	double thetasp1 = lamsp1-eps3/(1.0+eta)*C1*((1.0+eps2)*(1.0+eps2)+eta2+eta);
	double isp1 = 0.0;
	double q1sp1 = C1/(1.0+eta)*(((1.0+eps2)*(1.0+eps2)+eta2)*(q1+(1.0+eta)*cost)+(1.0+eps2)*((1.0+eta)*cost+q1*(eta-eps2)))-0.75*q2*(1.0-5.0*cosi2)*C2;
	double q2sp1 = C1/(1.0+eta)*(((1.0+eps2)*(1.0+eps2)+eta2)*(q2+(1.0+eta)*sint)+(1.0+eps2)*((1.0+eta)*sint+q2*(eta-eps2)))+0.75*q1*(1.0-5.0*cosi2)*C2;
	double Omegasp1 = 1.5*cosi*C2;

	double asp2 = -1.5*sini2/(a*eta6)*(1.0+eps2)*(1.0+eps2)*(1.0+eps2)*cos2t;
	C1 = 0.25*(3.0*(q1*sint+q2*cost)+3.0*sin2t+q1*sin3t-q2*cos3t)/(a*a*eta4);
	double lamsp2 = -0.75*eps3*sini2*cos2t/(a*a*eta4*(1.0+eta))*(1.0+eps2)*(2.0+eps2)-0.125*sini2/(a*a*eta2*(1.0+eta))*(3.0*(q1*sint+q2*cost)+q1*sin3t-q2*cos3t)-0.5*(3.0-5.0*cosi2)*C1;
	double thetasp2 = lamsp2-sini2/(32.0*a*a*eta4*(1.0+eta))*(36.0*q1*q2-4.0*(3.0*eta2+5.0*eta-1.0)*(q1*sint+q2*cost)+12.0*eps2*q1*q2-32.0*(1.0+eta)*sin2t-(eta2+12.0*eta+39.0)*(q1*sin3t-q2*cos3t)+36.0*q1*q2*cos4t-18.0*(q1*q1-q2*q2)*sin4t-3.0*(q1*q1-q2*q2)*q1*sin5t+3.0*(3.0*q1*q1-q2*q2)*q2*cos5t);
	double isp2 = -0.25*sini*cosi/(a*a*eta4)*(3.0*(q1*cost-q2*sint)+3.0*cos2t+q1*cos3t+q2*sin3t);
	double q1sp2 = 0.5*q2*(3.0-5.0*cosi2)*C1+0.125*sini2/(a*a*eta4)*(3.0*(eta2-q1*q1)*cost+3.0*q1*q2*sint-(eta2+3.0*q1*q1)*cos3t-3.0*q1*q2*sin3t)-3.0*sini2*cos2t/(16.0*a*a*eta4)*(10.0*q1+(8.0+3.0*q1*q1+q2*q2)*cost+2.0*q1*q2*sint+6.0*(q1*cos2t+q2*sin2t)+(q1*q1-q2*q2)*cos3t+2.0*q1*q2*sin3t);
	double q2sp2 = -0.5*q1*(3.0-5.0*cosi2)*C1-0.125*sini2/(a*a*eta4)*(3.0*(eta2-q2*q2)*sint+3.0*q1*q2*cost+(eta2+3.0*q2*q2)*sin3t+3.0*q1*q2*cos3t)-3.0*sini2*cos2t/(16.0*a*a*eta4)*(10.0*q2+(8.0+3.0*q1*q1+q2*q2)*sint+2.0*q1*q2*cost+6.0*(q1*sin2t-q2*cos2t)+(q1*q1-q2*q2)*sin3t-2.0*q1*q2*cos3t);
	double Omegasp2 = -cosi*C1;

	double soelp[6] = {alp, thetalp, ilp, q1lp, q2lp, Omegalp};
	double soesp1[6] = {asp1, thetasp1, isp1, q1sp1, q2sp1, Omegasp1};
	double soesp2[6] = {asp2, thetasp2, isp2, q1sp2, q2sp2, Omegasp2};
	for(int i=0;i<6;i++)
		soe_o[i] = soe_m[i]-j2*Ra[Ncenter]*Ra[Ncenter]*(soelp[i]+soesp1[i]+soesp2[i]);

	soe2coe(soe_o,ocoe);
}

//平均根数msoe转换为瞬时根数osoe
void M2O2(const double* soe_m, double* soe_o, double j2){
	double a     = soe_m[0];
	double theta = soe_m[1];
	double i     = soe_m[2];
	double q1    = soe_m[3];
	double q2    = soe_m[4];
	double Omega = soe_m[5];

	double cosi = cos(i);
	double sini = sin(i);
	double cosi2 = cosi*cosi;
	double sini2 = 1.0-cosi2;
	double cosi4 = cosi2*cosi2;
	double cosi6 = cosi4*cosi2;

	double cost = cos(theta);
	double sint = sin(theta);
	double cos2t = cost*cost-sint*sint;
	double sin2t = 2.0*sint*cost;
	double cos3t = cost*cos2t-sint*sin2t;
	double sin3t = sin2t*cost+cos2t*sint;
	double cos4t = cos2t*cos2t-sin2t*sin2t;
	double sin4t = 2.0*sin2t*cos2t;
	double cos5t = cos2t*cos3t-sin2t*sin3t;
	double sin5t = sin3t*cos2t+cos3t*sin2t;
	double eta2 = 1.0-q1*q1-q2*q2;
	double eta = sqrt(eta2);
	double eta4 = eta2*eta2;
	double eta6 = eta2*eta4;

	double eps1 = sqrt(1.0-eta2);
	double eps2 = q1*cost+q2*sint;
	double eps3 = q1*sint-q2*cost;
	double THETA = 1.0/(1.0-5.0*cosi2);
	double p = a*eta2;
	double r = p/(1.0+eps2);
	double beta = 1.0/(eta2+eta);
	double F = atan2(r*(1.0+beta*q1*q1)*sint-r*beta*q1*q2*cost+a*q2,r*(1.0+beta*q2*q2)*cost-r*beta*q1*q2*sint+a*q1);
	double theta_lam = theta-(F-q1*sin(F)+q2*cos(F));
	if(abs(theta_lam)>M_2PI)
		theta_lam = fmod(theta_lam, M_2PI);
	while(abs(theta_lam)>M_PI)
		theta_lam=theta_lam-theta_lam/abs(theta_lam)*M_2PI;
	//while(abs(theta_lam)>DPI)
	//	theta_lam=theta_lam-theta_lam/abs(theta_lam)*2.0*DPI;

	double alp = 0.0;
	double C1 = sini/(8.0*a*a*eta2)*(1.0-10.0*THETA*cosi2);
	double C2 = (q1*q2/(16.0*a*a*eta4))*(3.0-55.0*cosi2-280.0*THETA*cosi4-400.0*THETA*THETA*cosi6);
	double lamlp = q1*q2*sini/(1.0+eta)*C1+C2;
	double thetalp = lamlp-0.5*sini/eta2*C1*(q1*q2*(3.0+2.0*eta2/(1.0+eta))+2.0*(q1*sint+q2*cost)+0.5*eps1*sin2t);
	double ilp = 0.5*cosi/eta2*C1*(q1*q1-q2*q2);
	double q1lp = -0.5*q1*sini*C1-q2*C2;
	double q2lp = 0.5*q2*sini*C1+q1*C2;
	double Omegalp = q1*q2*cosi/(8.0*a*a*eta4)*(11.0+80.0*THETA*cosi2+200.0*THETA*THETA*cosi4);

	C1 = 0.25*(1.0-3.0*cosi2)/(a*a*eta4);
	C2 = (theta_lam+eps3)/(a*a*eta4);
	double asp1 = 2.0*a/eta2*C1*((1.0+eps2)*(1.0+eps2)*(1.0+eps2)-eta2*eta);
	double lamsp1 = eps3/(1.0+eta)*C1*((1.0+eps2)*(2.0+eps2)+eta2)+0.75*(1.0-5.0*cosi2)*C2;
	double thetasp1 = lamsp1-eps3/(1.0+eta)*C1*((1.0+eps2)*(1.0+eps2)+eta2+eta);
	double isp1 = 0.0;
	double q1sp1 = C1/(1.0+eta)*(((1.0+eps2)*(1.0+eps2)+eta2)*(q1+(1.0+eta)*cost)+(1.0+eps2)*((1.0+eta)*cost+q1*(eta-eps2)))-0.75*q2*(1.0-5.0*cosi2)*C2;
	double q2sp1 = C1/(1.0+eta)*(((1.0+eps2)*(1.0+eps2)+eta2)*(q2+(1.0+eta)*sint)+(1.0+eps2)*((1.0+eta)*sint+q2*(eta-eps2)))+0.75*q1*(1.0-5.0*cosi2)*C2;
	double Omegasp1 = 1.5*cosi*C2;

	double asp2 = -1.5*sini2/(a*eta6)*(1.0+eps2)*(1.0+eps2)*(1.0+eps2)*cos2t;
	C1 = 0.25*(3.0*(q1*sint+q2*cost)+3.0*sin2t+q1*sin3t-q2*cos3t)/(a*a*eta4);
	double lamsp2 = -0.75*eps3*sini2*cos2t/(a*a*eta4*(1.0+eta))*(1.0+eps2)*(2.0+eps2)-0.125*sini2/(a*a*eta2*(1.0+eta))*(3.0*(q1*sint+q2*cost)+q1*sin3t-q2*cos3t)-0.5*(3.0-5.0*cosi2)*C1;
	double thetasp2 = lamsp2-sini2/(32.0*a*a*eta4*(1.0+eta))*(36.0*q1*q2-4.0*(3.0*eta2+5.0*eta-1.0)*(q1*sint+q2*cost)+12.0*eps2*q1*q2-32.0*(1.0+eta)*sin2t-(eta2+12.0*eta+39.0)*(q1*sin3t-q2*cos3t)+36.0*q1*q2*cos4t-18.0*(q1*q1-q2*q2)*sin4t-3.0*(q1*q1-q2*q2)*q1*sin5t+3.0*(3.0*q1*q1-q2*q2)*q2*cos5t);
	double isp2 = -0.25*sini*cosi/(a*a*eta4)*(3.0*(q1*cost-q2*sint)+3.0*cos2t+q1*cos3t+q2*sin3t);
	double q1sp2 = 0.5*q2*(3.0-5.0*cosi2)*C1+0.125*sini2/(a*a*eta4)*(3.0*(eta2-q1*q1)*cost+3.0*q1*q2*sint-(eta2+3.0*q1*q1)*cos3t-3.0*q1*q2*sin3t)-3.0*sini2*cos2t/(16.0*a*a*eta4)*(10.0*q1+(8.0+3.0*q1*q1+q2*q2)*cost+2.0*q1*q2*sint+6.0*(q1*cos2t+q2*sin2t)+(q1*q1-q2*q2)*cos3t+2.0*q1*q2*sin3t);
	double q2sp2 = -0.5*q1*(3.0-5.0*cosi2)*C1-0.125*sini2/(a*a*eta4)*(3.0*(eta2-q2*q2)*sint+3.0*q1*q2*cost+(eta2+3.0*q2*q2)*sin3t+3.0*q1*q2*cos3t)-3.0*sini2*cos2t/(16.0*a*a*eta4)*(10.0*q2+(8.0+3.0*q1*q1+q2*q2)*sint+2.0*q1*q2*cost+6.0*(q1*sin2t-q2*cos2t)+(q1*q1-q2*q2)*sin3t-2.0*q1*q2*cos3t);
	double Omegasp2 = -cosi*C1;

	double soelp[6] = {alp, thetalp, ilp, q1lp, q2lp, Omegalp};
	double soesp1[6] = {asp1, thetasp1, isp1, q1sp1, q2sp1, Omegasp1};
	double soesp2[6] = {asp2, thetasp2, isp2, q1sp2, q2sp2, Omegasp2};
	for(int i=0;i<6;i++)
		soe_o[i] = soe_m[i]-j2*Ra[Ncenter]*Ra[Ncenter]*(soelp[i]+soesp1[i]+soesp2[i]);
}

//瞬时根数ocoe(a,e,i,Omega,omega,f)转换为平均根数mcoe(a,e,i,Omega,omega,m)
void O2M(const double* ocoe, double* mcoe, double j2){
	double soe_o[6], soe_m[6];

	coe2soe(ocoe,soe_o);

	const int n = 6;
	double fvec[6] = {0, 0, 0, 0, 0, 0};
	double para[7], wa[int((n*(3*n+13))/2)], xtol = 1.0e-6;
	for(int i=0;i<n;i++){
		para[i] = soe_o[i];
		soe_m[i] = soe_o[i];
	}
	para[6]=j2;
	int info = hybrd1(myfun, n, soe_m, fvec, para, wa, xtol, 0, 200);

	soe2coe(soe_m,mcoe);

	int flag;
	mcoe[5]=e2m(flag,f2e(flag,mcoe[5],mcoe[1]),mcoe[1]);
}

int myfun(int n, const double* x, double* fvec, int iflag, const double* para){
//para:瞬时根数	
	if(iflag==0)
		return 0;
	double soe_o[6];
	M2O2(x, soe_o, para[6]);
	for(int i=0;i<n;i++)
		fvec[i] = para[i]-soe_o[i];
	return 1;
}

/**************************************************************************************************************************************/
/************************************************************考虑j2摄动的轨道递推******************************************************/
/**************************************************************************************************************************************/
//a e i Omega omega m
void j2mcoe02mcoef(const double* me0, const double dt, double* mef)
{
	double a2=me0[0]*me0[0];
	double a3=a2*me0[0];
	double e2=me0[1]*me0[1];
	double n=sqrt(GMMu[Ncenter]/a3);
	double p=me0[0]*(1-e2);
	double p2=p*p;
	double JRPn=J2[Ncenter]*Ra[Ncenter]*Ra[Ncenter]/p2*n;
	double ci=cos(me0[2]);
	double ci2=ci*ci;
	double dOmega=-1.5*JRPn*ci;
	double domega=0.75*JRPn*(5*ci2-1);
	double dm=n+0.75*JRPn*sqrt(1-e2)*(3*ci2-1);
	mef[0]=me0[0];
	mef[1]=me0[1];
	mef[2]=me0[2];
	mef[3]=me0[3]+dOmega*dt;
	mef[4]=me0[4]+domega*dt;
	mef[5]=me0[5]+dm*dt;

	mef[3]=fmod(mef[3], M_2PI);
	mef[4]=fmod(mef[4], M_2PI);
	mef[5]=fmod(mef[5], M_2PI);
}

void j2ocoe02ocoef(double* ocoef, const double* ocoe0, const double dt)
{
	double mcoe0[6],mcoef[6];

	O2M(ocoe0,mcoe0);

	j2mcoe02mcoef(mcoe0, dt, mcoef);

	M2O(mcoef, ocoef);
}

void j2rv02rvf(const double* rv0, const double dt, double* rvf)
{
	int flag;
	double ocoe0[6],mcoe0[6],mcoef[6],ocoef[6];

	rv2coe(flag,ocoe0,rv0,GMMu[Ncenter]);
	O2M(ocoe0,mcoe0);

	j2mcoe02mcoef(mcoe0, dt, mcoef);

	M2O(mcoef, ocoef);
	coe2rv(flag, rvf, ocoef,GMMu[Ncenter]);
}
