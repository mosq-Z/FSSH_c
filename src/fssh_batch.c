#include "fssh_batch.h"
#include "fssh_thread_pool.h"

// 单个动量下批量运行所有轨迹，统计四类事件概率：
// TL：低能态透射、RL：低能态反射、TH：高能态透射、RH：高能态反射、other：超出步数未出边界
SingleKResult batch_run_single_k(double k)
{
    SingleKResult res;
    res.k = k;
    // 初始化统计计数器
    long count_TL = 0, count_RL = 0;
    long count_TH = 0, count_RH = 0;
    long count_other = 0;
    int total_traj = g_fssh_cfg.N_Single_K_Traj;

    // 遍历所有模拟轨迹
    for (long i = 0; i < total_traj; i++)
    {
        // 1. 初始化单条轨迹
        TrajContext traj;
        traj_init(&traj, k);

        // 2. 执行整条轨迹演化
        traj_run_full(&traj);

        // 3. 获取单条轨迹结果：res[0] 高低态，res[1] 透射反射
        int traj_res[2];
        traj_get_result(&traj, traj_res);
        int surface_flag = traj_res[0];
        int trans_flag = traj_res[1];

        // 分类统计
        if (traj.x > g_fssh_cfg.X_Max || traj.x < -g_fssh_cfg.X_Max)   //还给我的代码自动升级
        {
            if (surface_flag == 0)
            {
                if (trans_flag == 1)
                    count_TL++;
                else
                    count_RL++;
            }
            else
            {
                if (trans_flag == 1)
                    count_TH++;
                else
                    count_RH++;
            }
        }
        else
        {
            // 达到最大步数仍未跑出边界
            count_other++;
        }
    }

    // 转换为概率（总轨迹数归一化）
    double total = (double)total_traj;
    res.TL = (double)count_TL / total;
    res.RL = (double)count_RL / total;
    res.TH = (double)count_TH / total;
    res.RH = (double)count_RH / total;
    res.other = (double)count_other / total;

    return res;
}

// 遍历动量区间 [k1, k2]，步长gap，批量模拟并写入结果文件
void batch_run_k_range(double k1, double k2, double gap, const char* out_path)   //数组名是数组第一个元素的指针，这里的out_path是一个字符串名字，也是字符串第一个字符的指针
{
    // // 设置随机种子，保证每次批量模拟随机序列不同
    // srand((unsigned int)time(NULL));

    FILE* fp = fopen(out_path, "w");
    if (fp == NULL)
    {
        perror("批量结果文件打开失败");
        return;
    }

    // 写入表头：动量k 低能透射TL 低能反射RL 高能透射TH 高能反射RH 未出边界other
    fprintf(fp, "k\tTL\tRL\tTH\tRH\tother\n");

    // 循环遍历所有动量点
    for (double k = k1; k <= k2 + 1e-12; k += gap)
    {
        SingleKResult k_res = batch_run_single_k_parallel(k);
        // 格式化写入，保留6位小数
        fprintf(fp, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                k_res.k, k_res.TL, k_res.RL, k_res.TH, k_res.RH, k_res.other);

        printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
        k_res.k, k_res.TL, k_res.RL, k_res.TH, k_res.RH, k_res.other);
        
        // 控制台打印进度
        printf("动量 k = %.6f 模拟完成\n", k);
    }

    fclose(fp);
    printf("动量区间批量模拟全部完成，结果已保存至：%s\n", out_path);
}

// 从输入文件逐行读取动量值，批量模拟并输出统计结果
void batch_run_from_file(const char* in_path, const char* out_path)
{
    // srand((unsigned int)time(NULL));

    FILE* fp_in = fopen(in_path, "r");
    if (fp_in == NULL)
    {
        perror("动量输入文件打开失败");
        return;
    }

    FILE* fp_out = fopen(out_path, "w");
    if (fp_out == NULL)
    {
        perror("批量结果输出文件打开失败");
        fclose(fp_in);
        return;
    }

    fprintf(fp_out, "k\tTL\tRL\tTH\tRH\tother\n");

    double k;
    // 逐行读取文件中的动量数值
    while (fscanf(fp_in, "%lf", &k) == 1)
    {
        SingleKResult k_res = batch_run_single_k_parallel(k);
        fprintf(fp_out, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                k_res.k, k_res.TL, k_res.RL, k_res.TH, k_res.RH, k_res.other);

        printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
        k_res.k, k_res.TL, k_res.RL, k_res.TH, k_res.RH, k_res.other);

        printf("已完成动量 k = %.6f 的批量模拟\n", k);
    }

    fclose(fp_in);
    fclose(fp_out);
    printf("文件读取式批量模拟执行完毕，结果已输出\n");
}