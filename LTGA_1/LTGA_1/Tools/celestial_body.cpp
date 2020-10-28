#include "celestial_body.h"
#include "base_functions.h"

celestial_body::celestial_body(int t_id, double t_epoch, double t_a, double t_e, double t_i, double t_OMEGA, double t_omega, double t_f0,int center,int flag)
{
	m_id=t_id;m_epoch=t_epoch;m_a=t_a;m_e=t_e;m_inc=t_i;m_OMEGA=t_OMEGA;m_omega=t_omega;
	if(flag==0)
		m_f0=t_f0;
	if(flag==1)
	{
		int iflag=0;
		m_f0=e2f(iflag,m2e(iflag,t_f0,m_e),m_e);
	}
	m_center=center;
	m_n=sqrt(GMMu[m_center]/MuUnit/(m_a*m_a*m_a));
}
// 将mjd时刻天体的经典轨道根数赋值给coe
void celestial_body::GetCoe(double* coe, double mjd, double muNU)
{
	int flag;
	double dt=(mjd-m_epoch)*86400/TUnit;
	double coe0[6] = {m_a,m_e,m_inc,m_OMEGA,m_omega,m_f0};
	coe02coef(flag,coe,coe0,dt,muNU);
}
// 将mjd时刻天体的位置和速度赋值给rv
void celestial_body::GetRV(double* rv, double mjd, double muNU)
{
	int flag;
	double coe[6];
	GetCoe(coe,mjd,muNU);
	coe2rv(flag, rv, coe, muNU);
}