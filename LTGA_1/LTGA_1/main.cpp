#include <iostream>
#include <ctime>
#include "localConst.h"
#include "GA_fop.h"

celestial_body Mars(-1, 60389, 1.52363312, 0.09327933, 1.84785414*D2R, 49.48935357*D2R, 286.67090811*D2R, 328.88755274*D2R, 0, 0); // 真近点角f

void test_GA_fop()
{
	int flag, i;

	// **********
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
		// rv0[i] = rv0[i]/D2PI;
		// rv1[i] = rv1[i]/D2PI;
	}
	double rvm[6] = {0.0};
	Mars.GetRV(rvm, 59534, muNU);
	double m0 = 20000/MUnit;

	// 求解算法的一些参数设置
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[18] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	double tf = 2201*86400/TUnit, epsi = 1.0;
	flag = solve_GA_fop(Out, rv0, rv1, rvm, m0, tf, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("打靶变量值为:\n");
	for (i=1; i<=17; i++)
		printf("%.15e,\n", Out[i]);
	printf("引力辅助时间为%f天\n", Out[17]*TUnit/86400);
}

int main()
{
	printf("Hello!\n");

	clock_t start, stop;
	start = clock();
	
	test_GA_fop();

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}