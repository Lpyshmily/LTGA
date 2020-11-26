#include "dynamics_ee.h"
#include "constants.h"
#include "vector_operation.h"

//M:6*3, dM/dx:6*3*6, dD6/dx:6, ddM/dx:6*3*6*6, ddD6/dx:6*6
//step:1,只算M, D6; 2, 算到dM, dD6;3, 算到ddM, ddD6
void MatMD(double * M, double& D6, double * dM, double* dD6, double* ddM, double* ddD6, const double* ee, int step)
{
	int i=0, j=0;

	double p=ee[0];
	double f=ee[1];
	double g=ee[2];
	double h=ee[3];
	double k=ee[4];
	double L=ee[5];

	double H=sqrt(p/muNU);
	double cosL=cos(L);
	double sinL=sin(L);
	double W=1.0+f*cosL+g*sinL;
	double S=1.0+h*h+k*k;
	double G=h*sinL-k*cosL;
	double HW=H/W;
	double HWG=HW*G;

	for(i=0;i<18;i++)
			M[i]=0.0;

	M[0*3+1]=2.0*p*HW;
	M[1*3+0]=H*sinL;
	M[1*3+1]=HW*((1.0+W)*cosL+f);
	M[1*3+2]=-HWG*g;
	M[2*3+0]=-H*cosL;
	M[2*3+1]=HW*((1.0+W)*sinL+g);
	M[2*3+2]=HWG*f;
	M[3*3+2]=0.5*HW*S*cosL;
	M[4*3+2]=0.5*HW*S*sinL;
	M[5*3+2]=HWG;
	D6=W*W/(H*H*H*muNU);
	if(step<2) 
		return;

	double cos2L=cos(2.0*L);
	double sin2L=sin(2.0*L);	
	double dWdL=-f*sinL+g*cosL;
	double dGdL=h*cosL+k*sinL;	
	double halfp=0.5/p;
	double cosLW=cosL/W;
	double sinLW=sinL/W;
	double dWdLW=dWdL/W;

	int idM;
		
	for(i=0;i<108;i++)
		dM[i]=0.0;			
	for(i=0;i<6;i++)
		dD6[i]=0.0;
	
	// M中共有18个元素，每个元素分别对六个轨道根数求导，M中的一行3个元素求导得18个元素，iDM=行数*18+列数*6作为M中每个元素求导的初始位置
	// M12
	idM=0*18+1*6;
	dM[idM]=3.0*HW;
	dM[idM+1]=-cosLW*M[0*3+1];
	dM[idM+2]=-sinLW*M[0*3+1];
	dM[idM+5]=-dWdLW*M[0*3+1];
	// M21
	idM=1*18+0*6;
	dM[idM]=halfp*M[1*3+0];
	dM[idM+5]=H*cosL;
	// M22
	idM=1*18+1*6;
	dM[idM]=halfp*M[1*3+1];
	dM[idM+1]=0.5*HW*(3.0+cos2L)-cosLW*M[1*3+1];
	dM[idM+2]=0.5*HW*sin2L-sinLW*M[1*3+1];
	dM[idM+5]=HW*(g*cos2L-f*sin2L-2.0*sinL)-dWdLW*M[1*3+1];
	idM=1*18+2*6;
	dM[idM]=halfp*M[1*3+2];
	dM[idM+1]=-cosLW*M[1*3+2];
	dM[idM+2]=-HWG-sinLW*M[1*3+2];
	dM[idM+3]=-HW*g*sinL;
	dM[idM+4]=HW*g*cosL;
	dM[idM+5]=-HW*g*dGdL-dWdLW*M[1*3+2];
	idM=2*18+0*6;
	dM[idM]=halfp*M[2*3+0];
	dM[idM+5]=H*sinL;
	idM=2*18+1*6;
	dM[idM]=halfp*M[2*3+1];
	dM[idM+1]=0.5*HW*sin2L-cosLW*M[2*3+1];
	dM[idM+2]=0.5*HW*(3.0-cos2L)-sinLW*M[2*3+1];
	dM[idM+5]=HW*(f*cos2L+g*sin2L+2.0*cosL)-dWdLW*M[2*3+1];
	idM=2*18+2*6;
	dM[idM]=halfp*M[2*3+2];
	dM[idM+1]=HWG-cosLW*M[2*3+2];
	dM[idM+2]=-sinLW*M[2*3+2];
	dM[idM+3]=HW*sinL*f;
	dM[idM+4]=-HW*cosL*f;
	dM[idM+5]=HW*f*dGdL-dWdLW*M[2*3+2];
	idM=3*18+2*6;
	dM[idM]=halfp*M[3*3+2];
	dM[idM+1]=-cosLW*M[3*3+2];
	dM[idM+2]=-sinLW*M[3*3+2];
	dM[idM+3]=HW*h*cosL;
	dM[idM+4]=HW*k*cosL;
	dM[idM+5]=-0.5*HW*S*sinL-dWdLW*M[3*3+2];
	idM=4*18+2*6;
	dM[idM]=halfp*M[4*3+2];
	dM[idM+1]=-cosLW*M[4*3+2];
	dM[idM+2]=-sinLW*M[4*3+2];
	dM[idM+3]=HW*h*sinL;
	dM[idM+4]=HW*k*sinL;
	dM[idM+5]=0.5*HW*S*cosL-dWdLW*M[4*3+2];
	idM=5*18+2*6;
	dM[idM+0]=halfp*HWG;
	dM[idM+1]=-cosLW*HWG;
	dM[idM+2]=-sinLW*HWG;
	dM[idM+3]=HW*sinL;
	dM[idM+4]=-HW*cosL;
	dM[idM+5]=HW*dGdL-dWdLW*HWG;

	double C=2.0*D6/W;
	dD6[0]=-1.5*D6/p;
	dD6[1]=C*cosL;
	dD6[2]=C*sinL;
	dD6[5]=C*dWdL;
	if(step<3)
		return;

	int q=0, z=0;
	int iddM;
	double ddWdL=1.0-W;
	double ddWdLW=ddWdL/W;
	double ddGdL=-G;
	for(i=0;i<648;i++)
		ddM[i]=0.0;			
	for(i=0;i<36;i++)
		ddD6[i]=0.0;
	//ddM[i*108+j*36+k*6+q]=ddM[i*108+j*36+q*6+k]
	
	idM=0*18+1*6+0;//第1行2列1纵指标索引
	iddM=0*108+1*36+0*6;
	C=-dM[idM];
	ddM[iddM]=1.5/p*HW;
	ddM[iddM+1]=C*cosLW;
	ddM[iddM+2]=C*sinLW;
	ddM[iddM+5]=C*dWdLW;
	iddM=0*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+2]*cosLW-dM[idM+1]*sinLW;
	ddM[iddM+5]=-dM[idM+5]*cosLW+M[0*3+1]*sinLW-dM[idM+1]*dWdLW;	
	iddM=0*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-M[0*3+1]*cosLW-dM[idM+5]*sinLW-dM[idM+2]*dWdLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-M[0*3+1]*ddWdLW;
	idM=1*18+0*6+0;//第2行1列1纵指标索引	
	iddM=1*108+0*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+5]=halfp*dM[idM+5];
	ddM[iddM+5*6+5]=-H*sinL;
	idM=1*18+1*6+0;//第2行2列1纵指标索引
	iddM=1*108+1*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+5]=halfp*dM[idM+5];	
	iddM=1*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW-HW*sin2L+M[1*3+1]*sinLW-dM[idM+5]*cosLW;
	iddM=1*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+2]*dWdLW+HW*cos2L-M[1*3+1]*cosLW-dM[idM+5]*sinLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-2.0*HW*(g*sin2L+f*cos2L+cosL)-M[1*3+1]*ddWdLW;
	idM=1*18+2*6+0;//第2行3列1纵指标索引
	iddM=1*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=1*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+2]*cosLW-dM[idM+1]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+5]*cosLW+M[1*3+2]*sinLW-dM[idM+1]*dWdLW;
	iddM=1*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-HW*sinL-dM[idM+3]*sinLW;
	ddM[iddM+4]=HW*cosL-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-HW*dGdL-M[1*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=1*108+2*36+3*6+5;
	C=HW*g;
	ddM[iddM]=-dM[idM+3]*dWdLW-C*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW-C*sinL;
	ddM[iddM+2*6]=-2.0*dM[idM+5]*dWdLW-C*ddGdL-M[1*3+2]*ddWdLW;
	idM=idM+0;//第3行1列1纵指标索引
	iddM=2*108+0*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+5]=halfp*dM[idM+5];
	ddM[iddM+5*6+5]=H*cosL;
	idM=2*18+1*6+0;//第3行2列1纵指标索引
	iddM=2*108+1*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=2*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+HW*cos2L+M[2*3+1]*sinLW-dM[idM+5]*cosLW;
	iddM=2*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+2]*dWdLW+HW*sin2L-M[2*3+1]*cosLW-dM[idM+5]*sinLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-HW*2.0*(f*sin2L-g*cos2L+sinL)-M[2*3+1]*ddWdLW;
	idM=2*18+2*6+0;//第3行3列1纵指标索引
	iddM=2*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=2*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=HW*sinL-dM[idM+3]*cosLW;
	ddM[iddM+4]=-HW*cosL-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[2*3+2]*sinLW+HW*dGdL;
	iddM=2*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[2*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=2*108+2*36+3*6+5;
	ddM[iddM]=-dM[idM+3]*dWdLW+HW*f*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW+HW*f*sinL;
	ddM[iddM+2*6]=-2.0*dM[idM+5]*dWdLW+HW*ddGdL*f-M[2*3+2]*ddWdLW;
	idM=3*18+2*6+0;//第4行3列1纵指标索引
	iddM=3*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=3*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW-dM[idM+5]*cosLW+M[3*3+2]*sinLW;
	iddM=3*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[3*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=3*108+2*36+3*6+3;
	ddM[iddM]=HW*cosL;
	ddM[iddM+2]=-dM[idM+3]*dWdLW-HW*h*sinL;
	iddM=3*108+2*36+4*6+4;
	ddM[iddM]=HW*cosL;
	ddM[iddM+1]=-dM[idM+4]*dWdLW-HW*k*sinL;
	ddM[iddM+6+1]=-2.0*dM[idM+5]*dWdLW-0.5*HW*S*cosL-M[3*3+2]*ddWdLW;
	idM=4*18+2*6+0;//第5行3列1纵指标索引
	iddM=4*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=4*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[4*3+2]*sinLW-dM[idM+5]*cosLW;
	iddM=4*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[4*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=4*108+2*36+3*6+3;
	ddM[iddM]=HW*sinL;
	ddM[iddM+2]=-dM[idM+3]*dWdLW+HW*h*cosL;
	iddM=4*108+2*36+4*6+4;
	ddM[iddM]=HW*sinL;
	ddM[iddM+1]=-dM[idM+4]*dWdLW+HW*k*cosL;
	ddM[iddM+6+1]=-2.0*dM[idM+5]*dWdLW-0.5*HW*S*sinL-M[4*3+2]*ddWdLW;
	idM=5*18+2*6+0;//第6行3列1纵指标索引
	iddM=5*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=5*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[5*3+2]*sinLW-dM[idM+5]*cosLW;
	iddM=5*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[5*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=5*108+2*36+3*6+5;
	ddM[iddM]=-dM[idM+3]*dWdLW+HW*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW+HW*sinL;
	ddM[iddM+12]=-2.0*dM[idM+5]*dWdLW+HW*ddGdL-M[5*3+2]*ddWdLW;
	for(i=0;i<6;i++)
	{
		for(j=0;j<3;j++)
			for(q=0;q<6;q++)
				for(z=0;z<q;z++)
					ddM[i*108+j*36+q*6+z]=ddM[i*108+j*36+z*6+q];
	}

	ddD6[0*6+0]=-2.5/p*dD6[0];
	C=-1.5/p;
	ddD6[0*6+1]=C*dD6[1];
	ddD6[0*6+2]=C*dD6[2];
	ddD6[0*6+5]=C*dD6[5];
	ddD6[1*6+0]=2.0*dD6[0]*cosLW;
	ddD6[1*6+1]=dD6[1]*cosLW;
	ddD6[1*6+2]=-dD6[1]*sinLW+2.0*dD6[2]*cosLW;
	ddD6[1*6+5]=-dD6[1]*dWdLW-2.0*D6*sinLW+2.0*dD6[5]*cosLW;
	ddD6[2*6+0]=2.0*dD6[0]*sinLW;
	ddD6[2*6+1]=-dD6[2]*cosLW+2.0*dD6[1]*sinLW;
	ddD6[2*6+2]=-dD6[2]*sinLW+2.0*dD6[2]*sinLW;
	ddD6[2*6+5]=-dD6[2]*dWdLW+2.0*D6*cosLW+2.0*dD6[5]*sinLW;
	ddD6[5*6+0]=2.0*dD6[0]*dWdLW;
	ddD6[5*6+1]=2.0*dD6[1]*dWdLW-dD6[5]*cosLW-2.0*D6*sinLW;
	ddD6[5*6+2]=2.0*dD6[2]*dWdLW-dD6[5]*sinLW+2.0*D6*cosLW;
	ddD6[5*6+5]=dD6[5]*dWdLW+2.0*D6*ddWdLW;	
}

// 不会用到dfpara的值，设置成NULL即可
int dynamics_ee_top(double t, const double* x, double* dx, const double* dfpara)
{
	double M[18] = {0.0};
	double dM[108] = {0.0};
	double D6 = 0.0;
	double dD6[6] = {0.0};
	double ddM[1] = {0.0};
	double ddD6[1] = {0.0};
	MatMD(M, D6, dM, dD6, ddM, ddD6, x, 2);

	double v_lrv[6] = {0.0}; // 六个春分点轨道根数的协态
	V_Copy(v_lrv, &x[7], 6);

	// 根据六个春分点轨道根数的协态，求推力方向alpha并归一化
	double alpha[3] = {0.0}; // 推力方向的单位向量
	double norma;
	int i, j, k;
	for (j=0;j<3;++j)
	{
		for (i=0;i<6;++i)
			alpha[j] -= M[i*3+j]*v_lrv[i];
	}
	norma = V_Norm2(alpha, 3);
	// alpha归一化
	for (j=0;j<3;++j)
		alpha[j] /= norma;

	double tempd, u, m;
	u = 1.0; // 时间最优情况下，推力始终满开
	m = x[6]; // 质量
	for (i=0;i<14;++i)
		dx[i] = 0.0;
	tempd = TmaxNU*u/m;
	// 7个状态变量的导数
	for (i=0;i<6;++i)
	{
		for (j=0;j<3;++j)
			dx[i] += M[i*3+j]*alpha[j];
		dx[i] *= tempd;
	}
	dx[5] += D6;
	dx[6] = -TmaxNU*u/Ispg0NU; // 质量的导数
	// 7个协态变量的导数
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			for (k=0;k<3;++k)
				dx[7+i] -= v_lrv[j]*tempd*dM[j*18+k*6+i]*alpha[k];
		}
		dx[7+i] -= v_lrv[5]*dD6[i];
	}
	dx[13] = -TmaxNU*u/(m*m)*norma;
	return 1;
}

int dynamics_ee_fop(double t, const double* x, double* dx, const double* dfpara)
{
	double M[18] = {0.0};
	double dM[108] = {0.0};
	double D6 = 0.0;
	double dD6[6] = {0.0};
	double ddM[1] = {0.0};
	double ddD6[1] = {0.0};
	MatMD(M, D6, dM, dD6, ddM, ddD6, x, 2);

	double v_lrv[6] = {0.0}; // 六个春分点轨道根数的协态
	V_Copy(v_lrv, &x[7], 6);
	double lm = x[13]; // 质量协态

	double epsi, lam0;
	epsi = dfpara[0];
	lam0 = dfpara[1];

	// 根据六个春分点轨道根数的协态，求推力方向alpha并归一化
	double alpha[3] = {0.0}; // 推力方向的单位向量
	double norma;
	int i, j, k;
	for (j=0;j<3;++j)
	{
		for (i=0;i<6;++i)
			alpha[j] -= M[i*3+j]*v_lrv[i];
	}
	norma = V_Norm2(alpha, 3);
	// alpha归一化
	for (j=0;j<3;++j)
		alpha[j] /= norma;

	double tempd, u, m, rou;
	m = x[6]; // 质量
	rou=1.0-(Ispg0NU*norma/m+lm)/lam0; // 开关函数
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/(2*epsi);
	for (i=0;i<14;++i)
		dx[i] = 0.0;
	tempd = TmaxNU*u/m;
	// 7个状态变量的导数
	for (i=0;i<6;++i)
	{
		for (j=0;j<3;++j)
			dx[i] += M[i*3+j]*alpha[j];
		dx[i] *= tempd;
	}
	dx[5] += D6;
	dx[6] = -TmaxNU*u/Ispg0NU; // 质量的导数
	// 7个协态变量的导数
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			for (k=0;k<3;++k)
				dx[7+i] -= v_lrv[j]*tempd*dM[j*18+k*6+i]*alpha[k];
		}
		dx[7+i] -= v_lrv[5]*dD6[i];
	}
	dx[13] = -TmaxNU*u/(m*m)*norma;
	return 1;
}

