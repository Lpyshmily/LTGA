#include"PSO.h"
// 通过PSO使目标函数值最小
// ObjFun,目标函数
// xbest,fbest,优化结果
// D,变量维数
// Np,种群大小
void PSO(double (*ObjFun)(const double* X, const double* para), double* xbest, double& fbest, const double* para, int D, int Np, int ItMax, int ItOut, 
		 double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max, double Vmax)
{
	int ct;
	double Omega, C1, C2, val, GBVAL, rand0;	
	double* PBVAL=new double[Np];
	double* GBPOS=new double[D];
	double** popvel=NULL;
	popvel=new double*[Np];
	for(int i=0;i<Np;i++)
		popvel[i]=new double[D];
		
	double** pop=new double*[Np];
	for(int i=0;i<Np;i++)
		pop[i]=new double[D];
	double** PBPOS=new double*[Np];
	for(int i=0;i<Np;i++)
		PBPOS[i]=new double[D];	
	
	for(int j=0;j<Np;j++)
		PBVAL[j]= 1.0e10;
	GBVAL = 1.0e10;
	srand( (unsigned)time( NULL ) );
	for(int i=0;i<Np;i++)
	{   
		for(int j=0;j<D;j++)
		{
			rand0=(double)rand()/RAND_MAX;//第一个丢掉						
			pop[i][j] = (double)rand()/RAND_MAX;
			popvel[i][j] = 0.0;
		}	
	}
	
	ct=1;
	while(ct<=ItMax)
	{	
		#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<Np;i++)
		{
			val=ObjFun(pop[i], para);
			#pragma omp critical
			{
				if (val< PBVAL[i])
				{
					PBVAL[i]=val;
					for(int j=0;j<D;j++)
						PBPOS[i][j]=pop[i][j];				
				}
				if (val< GBVAL)
				{
					GBVAL=val;
					for(int j=0;j<D;j++)
						GBPOS[j]=pop[i][j];
				}
			}
		}
		Omega=OmegaMax-(OmegaMax-OmegaMin)/ItMax*ct;
		C1=-(C1Max-C1Min)*ct/ItMax+C1Max;
		C2=(C2Max-C2Min)*ct/ItMax+C2Min;
		for(int i=0;i<Np;i++)
		{			      	
			for(int j=0;j<D;j++)
			{
				rand0=(double)rand()/RAND_MAX;				
				popvel[i][j]=Omega*popvel[i][j]+C1*(PBPOS[i][j]-pop[i][j])*(double)rand()/RAND_MAX
					+C2*(GBPOS[j]-pop[i][j])*(double)rand()/RAND_MAX;
				if(popvel[i][j]>Vmax) popvel[i][j] =Vmax;//
				else if(popvel[i][j]<-Vmax) popvel[i][j] =-Vmax;//
				pop[i][j]+=popvel[i][j];
				if ((pop[i][j] > 1.0)||(pop[i][j] < 0.0))
				{
					rand0=(double)rand()/RAND_MAX;				
					pop[i][j] = (double)rand()/RAND_MAX;
					popvel[i][j]=0.0;
				}
			}
		}
		for(int j=0;j<D;j++)
			xbest[j]=GBPOS[j];
			fbest=GBVAL;
		if (ct%ItOut==0)
		{  		   
			 cout<<"No. of iteration="<<ct<<endl;	
 			 for(int i=0;i<D;i++)
				 cout<<"xbest("<<i+1<<")="<<setprecision(15)<<GBPOS[i]<<endl;			
			 cout<<"fbest="<<setprecision(15)<<GBVAL<<endl<<endl;			 
		}
		ct=ct+1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
	for(int i=0;i<Np;i++)
		delete[] popvel[i];
	delete[] popvel;
	for(int i=0;i<Np;i++)
		delete[] pop[i];
	delete[] pop;
	for(int i=0;i<Np;i++)
		delete[] PBPOS[i];		
	delete[] PBPOS;		
}