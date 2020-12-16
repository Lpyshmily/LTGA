#ifndef CELESTIAL_BODY_
#define CELESTIAL_BODY_
#include <string>
// using namespace std;

class celestial_body
{
public:
	celestial_body(int t_id, double t_epoch, double t_a, double t_e, double t_i, double t_OMEGA, double t_omega, double t_f0,int center,int flag=0);//flag:0=真近点角,1=平近点角
	void GetCoe(double* coe, double mjd, double muNU);
	void GetRV(double* rv, double mjd, double muNU);

	std::string name;
	int m_id;
	double m_epoch;
	double m_a;
	double m_e;
	double m_inc;
	double m_OMEGA;
	double m_omega;
	double m_f0;
	double m_n; // 平均角速度
	int m_center; // 中心天体编号
};

#endif