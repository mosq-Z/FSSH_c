// #include <stdio.h>
// #include <stdlib.h>
// #include <complex.h>
// #include "fssh_traj.h"
// #include "fssh_config.h"
#include "fssh_batch.h"


/*
// 打印d的代码  
int main(void)
{
    // 和你Python一致：x从-10到10，500个采样点
    const double x_start = -10.0;
    const double x_end = 10.0;
    const int N = 500;
    double dx = (x_end - x_start) / (N - 1);

    FILE *fp = fopen("output/result_data.txt", "w");
    if (!fp)
    {
        perror("文件打开失败");
        return 1;
    }

    // 表头：x E1 E2 d12_real d12_imag
    fprintf(fp, "x\tE1\tE2\td12_real\n");

    for (int i = 0; i < N; i++)
    {
        double x = x_start + i * dx;
        Mat2c H = H_di(x);
        Mat2c U;
        Vec2c eig;
        get_eig_info(&H, &U, &eig);

        double E1 = creal(eig.E1);
        double E2 = creal(eig.E2);
        double complex d12_val = d12(x);

        // 写入一行数据
        fprintf(fp, "%.6lf\t%.8lf\t%.8lf\t%.8lf\n",
                x, E1, E2, creal(d12_val)/50.0);
    }

    fclose(fp);
    printf("数据已保存到 result_data.txt，可使用Python绘图验证\n");
    return 0;

}
*/


// 一个可以展示中间过程的traj_run_full
/* int main(void)
 {
    
    // 1. 定义轨迹
    TrajContext traj;

    // 初始化参数：初始动量、初始坐标、初始振幅
    double init_k = 13;

    traj_init(&traj, init_k);

    // 2. 打开输出文件，追加写入模式
    FILE* fp = fopen("output/result_data.txt", "w");
    if (fp == NULL)
    {
        perror("无法打开output/result_data.txt");
        return EXIT_FAILURE;
    }
    srand((unsigned int)time(NULL));
    // 写入表头：步数 C1实部 C1虚部 C2实部 C2虚部 x 当前势能面 核受力 动量
    fprintf(fp, "step\tC1_real\tC1_imag\tC2_real\tC2_imag\tx\tsurface\tF\tk\tg12\n");

    // 先输出初始步数据
    fprintf(fp, "%ld\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n",
            traj.step,
            creal(traj.C_t.E1), cimag(traj.C_t.E1),
            creal(traj.C_t.E2), cimag(traj.C_t.E2),
            traj.x, traj.C_evo,
            traj.F, traj.k ,traj.g12);

    // 3. 循环演化，每一步输出数据，直到满足终止条件
    while (!traj_is_stop(&traj))
    {
        traj.step++;
        // 执行三步演化
        traj_step_II(&traj);
        traj_step_III(&traj);
        traj_step_IV(&traj);

        fprintf(fp, "%ld\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n",
                traj.step,
                creal(traj.C_t.E1), cimag(traj.C_t.E1),
                creal(traj.C_t.E2), cimag(traj.C_t.E2),
                traj.x, traj.C_evo,
                traj.F, traj.k ,traj.g12);
    }

    // 4. 关闭文件
    fclose(fp);
    printf("轨迹模拟结束，数据已保存至 output/result_data.txt\n");

    // 可选：打印最终结果
    int res[2];
    traj_get_result(&traj, res);
    printf("最终所处势能面：%d，透射(1)/反射(0)：%d\n", res[0], res[1]);

    return 0;
} */


int main(){
    fssh_cfg_set_ntraj(4000);
    batch_run_from_file("input_data","output/result_data.txt");
}