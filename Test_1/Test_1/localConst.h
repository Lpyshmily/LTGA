#ifndef LOCALCONST_H_
#define LOCALCONST_H_

#include "Tools\tools.h"

extern celestial_body Earth;
extern celestial_body Jupiter;

// 甩摆行星参数赋值
const double RU_Earth = 6378.145e3/LUnit; // 火星半径
const double rminU_Earth = 6878.145e3/LUnit; // 最小甩摆半径
const double muNU_Earth = GMMu[3]/GMMu[0]; // 火星引力系数

const double MaxNum = 1e7;
const int RepeatTime = 10; // 时间最优重复求解次数

#endif