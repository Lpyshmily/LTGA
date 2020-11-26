#ifndef LOCALCONST_H_
#define LOCALCONST_H_

#include "Tools\tools.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

// 甩摆行星参数赋值
const double RU_MARS = 3389.9e3/LUnit; // 火星半径
const double rminU_MARS = 3889.9e3/LUnit; // 最小甩摆半径
const double muNU_MARS = 42828.3e9/GMMu[0]; // 火星引力系数

const double MaxNum = 1e7;

#endif