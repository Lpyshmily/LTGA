#ifndef LOCALCONST_H_
#define LOCALCONST_H_

#include "Tools\tools.h"

extern celestial_body Earth;
extern celestial_body Jupiter;

// ˦�����ǲ�����ֵ
const double RU_Earth = 6378.145e3/LUnit; // ���ǰ뾶
const double rminU_Earth = 6878.145e3/LUnit; // ��С˦�ڰ뾶
const double muNU_Earth = GMMu[3]/GMMu[0]; // ��������ϵ��

const double MaxNum = 1e7;
const int RepeatTime = 10; // ʱ�������ظ�������

#endif