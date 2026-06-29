// fssh_traj.h
// 本文件定义一条路径的所有实例变量，生成实例，步进实例，终止实例，读取结果的对应分布函数，功能集中在修改实例变量上
#ifndef FSSH_TRAJ_H
#define FSSH_TRAJ_H

#include "H_quantities.h"


// 单条FSSH轨迹所有动态演化变量（完全等价Python实例self.xxx所有属性）
typedef struct TrajContext{
    // 1. 电子态、振幅
    int C_evo;            // 当前所处绝热态 C_evolution_t
    int C_evo_prev;       // 上一步绝热态 C_evolution_t0
    Vec2c C_t;            // 电子振幅向量

    // 2. 核动力学变量
    double x;             // 核坐标 x_Localtion_t
    double k;             // 动量 k_Momentum_t
    double F;             // 当前绝热态受力

    // 3. 势能、变换矩阵、耦合、有效哈密顿
    Mat2c H_di;
    Mat2c U;
    Vec2c H_ad;
    double complex d12;
    Mat2c H_eff;

    // 4. 跃迁相关
    double g12, g21, g_total;
    bool if_shift;        // 是否发生过跃迁
    long step;            // 当前演化步数

    // 5. 体系总能量
    double E_total;
} TrajContext;



void traj_init(TrajContext* traj, double k0);

// 三大核心步进函数（对应stepII、stepIII、stepIV）
void traj_step_II(TrajContext* traj);
void traj_step_III(TrajContext* traj);
void traj_step_IV(TrajContext* traj);

// 单条轨迹完整演化直到停止
void traj_run_full(TrajContext* traj);

// 判断轨迹是否终止
bool traj_is_stop(TrajContext* traj);

// 获取单条轨迹最终统计结果：[低能/高能, 透射/反射]
void traj_get_result(TrajContext* traj, int res[2]);

#endif