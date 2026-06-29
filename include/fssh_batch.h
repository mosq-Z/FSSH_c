// fssh_batch.h
// 本文件用来批量运行fssh_traj中单条轨迹的代码
# ifndef FSSH_BATCH_H
# define FSSH_BATCH_H
# include "fssh_config.h"
# include "fssh_traj.h"

// 单个动量的概率统计结果：k, TL, RL, TH, RH, others
typedef struct SingleKResult{
    double k;
    double TL, RL, TH, RH, other;
} SingleKResult;

// 单个动量下批量运行所有轨迹并统计概率
SingleKResult batch_run_single_k(double k);

// 从k1到k2按间隔遍历批量模拟
void batch_run_k_range(double k1, double k2, double gap, const char* out_path);

// 从test_input文件读取动量批量模拟
void batch_run_from_file(const char* in_path, const char* out_path);

#endif