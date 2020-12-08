#ifndef _CONSTANTS_
#define _CONSTANTS_
#include <math.h>

//������
const double M_PI=3.1415926535897932384626433832795;
const double M_2PI=6.283185307179586476925286766559;
const double M_PI2=1.570796326794897;
const double D2R=0.017453292519943295769236907684886;
const double R2D=57.295779513082320876798154814105;

//��������ʱ��
const double g0=9.80665;//m/s^2
const double AU = 1.49597870691e11;//m
const double JD_MJD  = 2400000.5;
const double TT_TAI = 32.184;
const double MJDJ2000 = 51544.5;
const double JC2JD= 36525.0;
const double JD2Second=86400.0;

//Ĭ�����������ţ�Ӱ��ƽ˲�������ת��
const int Ncenter = 3;

//�������Ǳ�ţ����峣���Ͱ뾶���ο���STK�ļ�
//0=̫����1=ˮ�ǣ�2=���ǣ�3=����4=���ǣ�5=ľ�ǣ�6=���ǣ�7=�����ǣ�8=�����ǣ�9=ڤ���ǣ�10=����
const double GMMu[11] ={1.32712440018e+20,  2.203209000000e+13,  3.248585920790e+14,  
	                  3.986004415000e+14,  4.282837564100e+13,  1.267127648382e+17,
	                  3.794058536168e+16,  5.794557628118e+15,  6.836534878892e+15,
					  9.769998557980e+11,  4.902800305555e+12   }; //m^3/s^2��Ĭ��Ϊ����ϵͳ����������

const double Ra[11] ={695508000.0,  2439000.0,  6051000.0,  
	                  6378136.3,    3396000.0,  71492000.0,
					  60330000.0,   26200000.0, 25225000.0,
					  1151000.0,    1738000.0                   }; //m

const double J2[11] ={0.0,        6.0e-5,      4.404435e-6,
	                  0.0010826,  0.0019526,   0.014696,       
					  0.016291,   0.0033442,   0.0034105,
					  0.0,        0.20333e-3                    };

const double LUnit=AU; // һ�����ĵ�λ
const double TUnit=sqrt(LUnit*LUnit*LUnit/GMMu[0]); // �볤��Ϊ���ȵ�λʱ����������ڴ�ʱ�䵥λ��Ϊ2PI
const double MuUnit=LUnit/TUnit*LUnit/TUnit*LUnit;
const double VUnit=LUnit/TUnit;
const double AUnit=VUnit/TUnit;
const double Isp=6000;//s
const double Tmax=2.26;//N
const double MUnit=19820;//kg
const double Ispg0=Isp*9.80665;//m/s
const double Ispg0NU=Ispg0/VUnit;
const double muNU=(GMMu[0]/LUnit)*(TUnit/LUnit)*(TUnit/LUnit); // ��ֵΪ1
const double TmaxNU=Tmax/(MUnit*AUnit);

#endif