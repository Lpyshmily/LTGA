#ifndef LOCALCONST_H_
#define LOCALCONST_H_

#include "Tools\tools.h"

extern celestial_body Mars;

// 甩摆行星参数赋值
const double RU_MARS = 3389.9e3/LUnit; // 火星半径
const double rminU_MARS = 3889.9e3/LUnit; // 最小甩摆半径
const double muNU_MARS = 42828.3e9/GMMu[0]; // 火星引力系数

const double MaxNum = 1e7;
const int RepeatTime = 10; // 时间最优重复求解次数

#endif