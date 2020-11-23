#include"PSO.h"
// 通过PSO使目标函数值最小
// ObjFun,目标函数
// xbest,fbest,优化结果
// D,变量维数
// Np,种群大小
// wa，工作变量，大小至少为3*Np*D+Np+D
// 依次存储PBVAL(Np),GBPOS(D),popvel(Np*D),pop(Np*D),PBPOS(Np*D)
// 分别表示每个粒子的最优目标函数值,全局最优位置,每个粒子的速度,每个粒子的位置,每个粒子的最优位置
void PSO(double (*ObjFun)(const double* X, const double* para), double* xbest, double& fbest,const double* para,
		 int D, int Np, double* wa, int ItMax, int ItOut, double OmegaMin, double OmegaMax, double C1Min, double C1Max,
		 double C2Min, double C2Max, double Vmax)
{
	int ct, i, j, index;
	int s[5]={0,Np,Np+D,Np*D+Np+D,2*Np*D+Np+D};//PBVAL,GBPOS,popvel,pop,PBPOS的起始索引
	double Omega, C1, C2, val, GBVAL, rand0;	
	
	for(int j=0;j<Np;j++)
		wa[s[0]+j]= 1.0e10;//PBVAL[j]= 1.0e10;
	GBVAL = 1.0e10;

	// **********
	/*
	// 设置xbest,fbest为初始近优解
	for(int j=0;j<D;j++)
		wa[s[1]+j]= xbest[j];//PBVAL[j]= xbest[j];
	fbest = ObjFun(xbest, para);
	GBVAL = fbest;
	*/
	// **********

	for(i=0;i<Np;i++)
	{   
		for(j=0;j<D;j++)
		{
			index=i*D+j;
			rand0=(double)rand()/RAND_MAX;//第一个丢掉						
			wa[s[3]+index] = (double)rand()/RAND_MAX;//pop[i][j] = (double)rand()/RAND_MAX;
			wa[s[2]+index] = 0.0;//popvel[i][j] = 0.0;
		}	
	}
	
	ct=1;
	while(ct<=ItMax)
	{	
		#pragma omp parallel for schedule(dynamic)
		for(i=0;i<Np;i++)
		{
			index=i*D;
			val=ObjFun(&wa[s[3]+index],para);//val=ObjFun(pop[i],para);
			#pragma omp critical
			{
				if (val< wa[s[0]+i])//if (val< PBVAL[i])
				{
					wa[s[0]+i]=val;//PBVAL[i]=val;
					for(j=0;j<D;j++)
						wa[s[4]+index+j]=wa[s[3]+index+j];//PBPOS[i][j]=pop[i][j];
				
				}
				if (val< GBVAL)
				{
					GBVAL=val;
					for(j=0;j<D;j++)
						wa[s[1]+j]=wa[s[3]+index+j];//GBPOS[j]=pop[i][j];
				}
			}
		}
		Omega=OmegaMax-(OmegaMax-OmegaMin)/ItMax*ct;
		C1=-(C1Max-C1Min)*ct/ItMax+C1Max;
		C2=(C2Max-C2Min)*ct/ItMax+C2Min;
		for(i=0;i<Np;i++)
		{			      	
			for(j=0;j<D;j++)
			{
				index=i*D+j;
				rand0=(double)rand()/RAND_MAX;
				wa[s[2]+index]=Omega*wa[s[2]+index]+C1*(wa[s[4]+index]-wa[s[3]+index])*(double)rand()/RAND_MAX
					+C2*(wa[s[1]+j]-wa[s[3]+index])*(double)rand()/RAND_MAX;
				if(wa[s[2]+index]>Vmax) wa[s[2]+index] =Vmax;//
				else if(wa[s[2]+index]<-Vmax) wa[s[2]+index] =-Vmax;//
				wa[s[3]+index]+=wa[s[2]+index];
				if ((wa[s[3]+index] > 1.0)||(wa[s[3]+index] < 0.0))
				{
					rand0=(double)rand()/RAND_MAX;				
					wa[s[3]+index] = (double)rand()/RAND_MAX;
					wa[s[2]+index]=0.0;
				}
//				popvel[i][j]=Omega*popvel[i][j]+C1*(PBPOS[i][j]-pop[i][j])*(double)rand()/RAND_MAX
//					+C2*(GBPOS[j]-pop[i][j])*(double)rand()/RAND_MAX;
//				if(popvel[i][j]>Vmax) popvel[i][j] =Vmax;//
//				else if(popvel[i][j]<-Vmax) popvel[i][j] =-Vmax;//
//				pop[i][j]+=popvel[i][j];
//				if ((pop[i][j] > 1.0)||(pop[i][j] < 0.0))
//				{
//					rand0=(double)rand()/RAND_MAX;				
//					pop[i][j] = (double)rand()/RAND_MAX;
//					popvel[i][j]=0.0;
//				}
			}
		}
		for(j=0;j<D;j++)
			xbest[j]=wa[s[1]+j];//xbest[j]=GBPOS[j];
		fbest=GBVAL;
		if (ItOut>0&&ct%ItOut==0)
		{  		   
			 cout<<"No. of iteration="<<ct<<endl;	
 			 for(i=0;i<D;i++)
	//			 cout<<"xbest("<<i+1<<")="<<setprecision(15)<<GBPOS[i]<<endl;
				cout<<setprecision(15)<<wa[s[1]+i]<<","<<endl;//cout<<setprecision(15)<<GBPOS[i]<<","<<endl;
			 cout<<"fbest="<<setprecision(15)<<GBVAL<<endl<<endl;			 
		}
		ct=ct+1;
	}
/*	ofstream myf("results.txt");
	if(myf.fail())
	{
		cout<<"文件打开失败"<<endl;
		exit(1);
	}
	for(i=0;i<D;i++)
		myf<<setprecision(15)<<GBPOS[i]<<endl;
	myf.close();*/

}

