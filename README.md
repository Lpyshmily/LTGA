# LTGA

Low-thrust Gravity-assist

* [LTGA_1](#LTGA_1)

* [LTGA_2](#LTGA_2)

* [LTGA_3](LTGA_3.md)

* [LTGA_4](LTGA_4.md)

* [LTGA_5](LTGA_5.md)

* [LTGA_6](LTGA_6.md)

* [Test_1](Test_1.md)

* LTGA_1

    用传统间接法求解小推力引力辅助问题。

    地球-火星-木星，燃料最优，交会。

    采用之前得到的打靶变量作为初值，求解同伦参数`epsi=1.0`的情况，可以收敛，但是需要保留比较多的小数位数，才能达到精度要求。

    如果只将收敛解中的引力辅助时间和引力辅助速度增量作为初值，其余13个乘子采用随机猜测，很难收敛。

* LTGA_2

    将引力辅助问题分成两段求解，引力辅助时刻作为一个总体的优化变量，两端轨迹分别为燃料最优飞越问题和燃料最优交会问题。

    采用直角坐标描述。

    计算第二段燃料最优问题之前先计算时间最优问题，末端位置固定。

    在LTGA_1中，真实最优解的比例系数为0.375，剩余质量为 15735.957 kg。

    在本程序中，如果比例系数为0.375，不能完成转移过程；比例系数为0.35时，可以完成两段转移，剩余质量为 15395.208 kg。

    改变比例系数，生成多组结果。`info.txt`

    设置同伦参数为1.0E-5，模拟燃料最优，同样生成多组结果。 `info_2.txt`

    添加第一段求解时间最优飞越问题。但是很多情况下时间最优问题的结果会比较大，应该是收敛到了局部最优解，需要进行多次求解选出全局最优解。

    * 20201103

        1. 对程序进行完善，还没有添加重复求解时间最优问题的功能

        2. 加入重复求解时间最优问题的功能，运行结果为[info_3.txt](LTGA_2/LTGA_2/info_3.txt)。
        
            设置同伦参数为1.0E-5，当比例系数为0.35时，剩余质量最大，为 15822.030198 kg。

            注意：当比例系数为0.43时，单独有一组可行解。搜索是否不够充分？
        
        3. 将时间最优问题的重复求解次数增加为20次
        
            运行结果为[info_4.txt](LTGA_2/LTGA_2/info_4.txt)，和[info_3.txt](LTGA_2/LTGA_2/info_3.txt)相同。

            因此，后续都采用重复10次。
        
        4. 在比例系数为0.34~0.36的范围内，以0.001为间隔，进行更详细的搜索。

            运行结果为[info_5.txt](LTGA_2/LTGA_2/info_5.txt)

            比例系数为0.347时，剩余质量最大，为 15839.872773 kg。
    
    * 20201105
        
        1. 将factor作为目标函数的参数，使用搜索算法计算最优解

            搭建搜索算法的目标函数`test_GA_obj()`，用比例系数0.34、0.37测试成功。

            由于PSO搜索程序是寻找目标函数的最小值，而我们要求的是剩余质量最大，因此目标函数的返回值应该是剩余质量的相反数。

            取消文件输出，进行PSO搜索，采用默认参数设置。

            采用**多核并行加快计算速度**。
    
    * 20201111

        PSO迭代500次运行结果：比例系数为0.346872586479787，剩余质量为15839.8971406158kg。相较于比例系数为0.347的结果，PSO算法的提升较小。