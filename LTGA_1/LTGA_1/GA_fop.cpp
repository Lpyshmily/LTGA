#include "GA_fop.h"
#include "localConst.h"

// n 打靶变量个数 17
/* x[17] 打靶变量
[0] lam0
[1-7] 初始协态变量
[8-11] 等式约束乘子
[12] 不等式约束乘子
[13-15] 引力辅助速度增量dvg
[16] 从出发到引力辅助的时间tm */
/* fvec[17] 打靶偏差 此函数的输出量
[0-5] 末端状态约束
[6] 末端质量协态为零
[7-9] 引力辅助时的位置约束
[10] 引力辅助前后相对速度大小相等
[11] 不等式约束的互补松弛条件
[12-14] tm-时刻速度协态横截条件
[15] 引力辅助静态条件
[16] 协态归一化 */
// iflag为0时表示上方调用函数告知不必算残差
/* sfpara[22]
[0-5] rv0
[6] m0
[7-12] 出发时刻目标天体的位置速度
[13-18] 末端时刻引力辅助天体的位置速度
[19] 转移时间tf
[20] epsi
[21] 是否输出的标志，0表示不输出，1表示输出 */
// 计算成功时，返回0或大于0的值，不成功时返回小于0的值
int fvec_GA_fop_rend(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if (iflag==0)
		return 0;

	int i;

	double x0[14] = {0.0}, rvf[6] = {0.0}, rvm[6] = {0.0};
	V_Copy(x0, sfpara, 7);
	V_Copy(rvf, &sfpara[7], 6);
	V_Copy(rvm, &sfpara[13], 6);
	double tf = sfpara[19];
	double epsi = sfpara[20];
	int outflag = NINT(sfpara[21]);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7);
	double chi[4] = {0.0}, kappa, dvg[3] = {0.0}, tm;
	V_Copy(chi, &x[8], 4);
	kappa = x[12];
	V_Copy(dvg, &x[13], 3);
	tm = x[16];

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;

	double AbsTol[14] = {0.0};
	for(i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;
	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	double work[140] = {0.0};
	int flag, NumPoint;
	flag = ode45(dynamics_rv_fop, x0, dfpara, 0.0, tm, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid1);
	// 计算一部分偏差
	// 引力辅助位置约束,7-9
	double rvm_tm[6] = {0.0};
	rv02rvf(flag, rvm_tm, rvm, tm, muNU);
	V_Minus(&fvec[7], x0, rvm_tm, 3);
	// normvinfin==normvinfout,10
	double vinfin[3], vinfout[3], normvinfin, normvinfout, unit_in[3], unit_out[3];
	V_Minus(vinfin, &x0[3], &rvm_tm[3], 3);
	V_Add(vinfout, vinfin, dvg, 3);
	normvinfin = V_Norm2(vinfin, 3);
	normvinfout = V_Norm2(vinfout, 3);
	V_Divid(unit_in, vinfin, normvinfin, 3);
	V_Divid(unit_out, vinfout, normvinfout, 3);
	fvec[10] = normvinfin - normvinfout;
	// 互补松弛条件,11
	double theta = acos(V_Dot(unit_in, unit_out, 3));
	double rp = muNU_MARS/(normvinfin*normvinfout) * (1/sin(theta/2) - 1)/rminU_MARS; // rp/rmin
	fvec[11] = kappa*(1 - rp);

	double A[3], B[3], C[3], temp, c;
	temp = 1/( 4*sin(theta/2)*sin(theta/2) * (1-sin(theta/2)) );
	for (i=0;i<3;++i)
	{
		A[i] = rp/normvinfin*(temp*(unit_out[i] - cos(theta)*unit_in[i]) - unit_in[i]);
		B[i] = rp/normvinfout*(temp*(unit_in[i] - cos(theta)*unit_out[i]) - unit_out[i]);
		C[i] = rp*( -temp*(1/normvinfout*(unit_in[i] - cos(theta)*unit_out[i]) + 1/normvinfin*(unit_out[i] - cos(theta)*unit_in[i])) + unit_out[i]/normvinfout + unit_in[i]/normvinfin );
	}
	double am_tm[3];
	double rm_tm = V_Norm2(rvm_tm, 3);
	for (i=0;i<3;++i)
		am_tm[i] = -muNU/(rm_tm*rm_tm*rm_tm)*rvm_tm[i];
	c = V_Dot(C, am_tm, 3);

	// tm-时刻速度协态,12-14
	for (i=0;i<3;++i)
		fvec[12+i] = x0[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	// 静态条件偏差,15
	double H1 =  ham_rv_fop(x0, epsi, lam0);
	// 计算新的状态变量和协态变量
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // 速度
		x0[7+i] = x0[7+i] - chi[i]; // 位置协态
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // 速度协态
	}
	double H2 = ham_rv_fop(x0, epsi, lam0);
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rvm_tm[3], 3) + chi[3]*V_Dot(tempu, am_tm, 3) - kappa*c;

	// 第二阶段的积分
	flag = ode45(dynamics_rv_fop, x0, dfpara, tm, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	// 末端状态约束,0-5
	V_Minus(fvec, x0, rvf, 6);
	// 末端质量协态,6
	fvec[6] = x0[13];
	fvec[16] = V_Norm2(x, 13) - 1.0;
	if(outflag>0)
	{
		fvec[0]=x0[6];
	}
	return 0;
}

// 将16个随机猜测值转换为17个打靶变量
// [0]:lam0
// [1-7]:7,初始时刻协态变量
// [8-11]:4,等式约束乘子
// [12]:不等式约束乘子
// [13-15]:3,速度增量
// [16]:引力辅助时间
void guess2lam_GA_fop(const double* guessValue, double* lamValue)
{
	double temp1, temp2;
	int i;
	// 用1个猜测值作为前8个打靶变量的norm，然后根据协态归一化条件计算出5个乘子的norm
	double norm1, norm2;
	norm1 = guessValue[0];
	norm2 = sqrt(1 - norm1*norm1);
	// 用7个猜测值，结合norm1，产生前8个打靶变量
	double alpha1[7];
	V_Copy(alpha1, &guessValue[1], 7);
	for (i=0;i<3;++i)
		alpha1[i]  = alpha1[i]*0.5*M_PI;
	alpha1[3] = (alpha1[3] - 0.5)*M_PI;
	alpha1[4] = alpha1[4]*M_2PI;
	alpha1[5] = (alpha1[5] - 0.5)*M_PI;
	alpha1[6] = alpha1[6]*M_2PI;
	
	lamValue[0] = norm1*sin(alpha1[0]);
	temp1 = norm1*cos(alpha1[0])*cos(alpha1[1]);
	temp2 = temp1*sin(alpha1[2]);
	lamValue[1] = temp2*cos(alpha1[3])*cos(alpha1[4]);
	lamValue[2] = temp2*cos(alpha1[3])*sin(alpha1[4]);
	lamValue[3] = temp2*sin(alpha1[3]);
	temp2 = temp1*cos(alpha1[2]);
	lamValue[4] = temp2*cos(alpha1[5])*cos(alpha1[6]);
	lamValue[5] = temp2*cos(alpha1[5])*sin(alpha1[6]);
	lamValue[6] = temp2*sin(alpha1[5]);
	
	lamValue[7] = norm1*cos(alpha1[0])*sin(alpha1[1]);
	// 用4个猜测值，结合norm2，产生5个约束的乘子
	double alpha2[4];
	V_Copy(alpha2, &guessValue[8], 4);
	alpha2[0] = alpha2[0]*0.5*M_PI;
	alpha2[1] = (alpha2[1] - 0.5)*M_PI;
	alpha2[2] = (alpha2[2] - 0.5)*M_PI;
	alpha2[3] = alpha2[3]*M_2PI;
	temp1 = norm2*cos(alpha2[0]);
	temp2 = temp1*cos(alpha2[1]);
	lamValue[8] = temp2*cos(alpha2[2])*cos(alpha2[3]);
	lamValue[9] = temp2*cos(alpha2[2])*sin(alpha2[3]);
	lamValue[10] = temp2*sin(alpha2[2]);
	lamValue[11] = temp1*sin(alpha2[1]);
	lamValue[12] = norm2*sin(alpha2[0]);
	// 用1个猜测值，产生速度增量幅值
	double vAmplitude = sqrt(muNU_MARS/rminU_MARS)*guessValue[12];
	// 用2个猜测值，结合vAmplitude，产生3个速度分量
	double phi, delta;
	phi = (guessValue[13] - 0.5)*M_PI;
	delta = guessValue[14]*M_2PI;
	lamValue[13] = vAmplitude*cos(phi)*cos(delta);
	lamValue[14] = vAmplitude*cos(phi)*sin(delta);
	lamValue[15] = vAmplitude*sin(phi);
	// 用1个猜测值，产生引力辅助时刻
	lamValue[16] = guessValue[15]*2201*86400/TUnit;
}

// x[17]:[0],lam0;[1-7]初始协态;[8-11]等式约束乘子;[12]不等式约束乘子;[13-15]引力辅助速度增量dvg;[16]引力辅助时刻
// Out[18]:[0],剩余质量;[1-17],17个打靶变量 
int solve_GA_fop(double* Out, const double* rv0, const double* rv1, const double* rvm, double m0, double tf, double epsi, int MaxGuessNum)
{
	double sfpara[22] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rv1, 6);
	V_Copy(&sfpara[13], rvm, 6);
	sfpara[19] = tf;
	sfpara[20] = epsi;
	sfpara[21] = 0.0;

	int n = 17;
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;

	int num = 0;
	int j, info, flag = 0;
	while (num<MaxGuessNum)
	{
		for (j=0; j<16; ++j)
			guessArray[j] = (double)rand()/RAND_MAX;
		guess2lam_GA_fop(guessArray, x);
		x[13] = 7.823746734334966e-002;
		x[14] = 7.901755272921728e-002;
		x[15] = -5.190057707485141e-003;
		x[16] = 1.420237177534943e+001;
		/*
		x[0] = 6.137291124645334e-001;
		x[1] = -2.792803801659735e-001;
		x[2] = -4.608278012924570e-001;
		x[3] = -5.363877225501464e-002;
		x[4] = 3.636848762039098e-001;
		x[5] = -3.348964164277419e-001;
		x[6] = -5.564392847421700e-002;
		x[7] = 1.768777878142729e-001;
		x[8] = -7.358387075241563e-003;
		x[9] = -1.036134836461006e-001;
		x[10] = 6.240424655580455e-002;
		x[11] = -1.905572908013827e-001;
		x[12] = 1.729442111511000e-002;
		x[13] = 7.823746734334966e-002;
		x[14] = 7.901755272921728e-002;
		x[15] = -5.190057707485141e-003;
		x[16] = 1.420237177534943e+001;
		*/
		info = hybrd1(fvec_GA_fop_rend, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
		// std::cout << enorm(n,fvec) << std::endl;
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[21]=1.0;
			int _j = fvec_GA_fop_rend(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
			if(fvec[0]>0.0) // 剩余质量为正，且满足打靶精度要求，停止
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, n);
				break;
			}
			sfpara[21]=0.0;
		}
		num++;
	}
	printf("随机猜测次数：%d\n", num);
	return flag;
}