// 本文件定义一条路径的所有实例变量，生成实例，步进实例，终止实例，读取结果的对应分布函数，功能集中在修改实例变量上
#include "fssh_traj.h"
#include "fssh_config.h"



// 生成并初始化实例
void traj_init(TrajContext* traj, double k0){
    // 1. 电子态、电子振幅初始化
    double x0 = - g_fssh_cfg.X_Max;
    Vec2c C0 = {1.0 + 0.0*I ,0.0 + 0.0*I};
    traj->C_evo = (int)creal(C0.E1);       // C_evolution_t：当前所处绝热态
    traj->C_evo_prev = traj->C_evo;      // C_evolution_t0：上一步初始态
    traj->C_t = C0;                      // 初始电子复振幅向量

    // 2. 核动力学初始变量
    traj->x = - g_fssh_cfg.X_Max;                        // 初始核坐标
    traj->k = k0;                        // 初始动量
    // 根据初始坐标初始化透热哈密顿
    traj->H_di = H_di(x0);
    // 对角化得到酉变换矩阵U、绝热能级H_ad
    get_eig_info(&traj->H_di, &traj->U, &traj->H_ad);
    // 初始非绝热耦合d12
    traj->d12 = d12(x0);
    // 构造有效哈密顿
    traj->H_eff = H_eff(traj->H_ad, traj->k, g_fssh_cfg.M_Quality, traj->d12);
    // 初始受力：基于初始绝热态+坐标
    traj->F = F(x0, traj->C_evo);

    // 3. 跃迁相关参数初始化
    traj->g12 = 0.0;
    traj->g21 = 0.0;
    traj->g_total = 0.0;
    traj->if_shift = false;             // 初始未发生任何面跃迁
    traj->step = 0L;                    // 初始步数为0

    // 4. 体系初始总能量：动能 + 当前绝热态势能
    traj->E_total = E_total(traj->k, g_fssh_cfg.M_Quality, traj->C_evo, traj->H_ad);
}

// 三大核心步进函数（对应stepII、stepIII、stepIV）
void traj_step_II(TrajContext* traj){
    double temp_F;
    double delta_t = g_fssh_cfg.Delta_T;
    double m = g_fssh_cfg.M_Quality;

    // 1.更新核坐标
    traj->x = x_change_by_t(traj->k, m, traj->x, traj->F, delta_t);
    // 2.更新核动量
    traj->k = k_change_by_t(traj->k, traj->F, delta_t);

    // 3.刷新当前坐标下所有x依赖物理量
    traj->H_di = H_di(traj->x);
    get_eig_info(&traj->H_di, &traj->U, &traj->H_ad);
    traj->d12 = d12(traj->x);
    traj->H_eff = H_eff(traj->H_ad, traj->k, m, traj->d12);

    // 临时获取当前位置受力，用于电子振幅演化
    temp_F = F(traj->x, traj->C_evo);

    // 4.RKG四阶方法演化电子振幅
    traj->C_t = C_change_by_t_RKG(traj->k, m, traj->x, temp_F, traj->H_eff, traj->C_t, delta_t);

    // 5.更新体系总能量
    traj->E_total = E_total(traj->k, m, traj->C_evo, traj->H_ad);
}


void traj_step_III(TrajContext* traj){
    double delta_t = g_fssh_cfg.Delta_T;
    double m = g_fssh_cfg.M_Quality;
    // 两次rand组合，把精度从15位提升到30位，最小间隔 ≈ 9.3e-10
    // double xi = (double)rand() / RAND_MAX
    double xi = ( (double)rand() * (RAND_MAX + 1) + rand() ) / ( (double)(RAND_MAX + 1) * RAND_MAX );

    if (traj->C_evo == 1)
    {
        traj->g12 = g12(traj->C_t, traj->x, traj->k, m, delta_t);
        traj->g_total += traj->g12;
        if (traj->g12 > xi)
        {
            traj->C_evo = 0;
        }
    }
    else
    {
        traj->g21 = g21(traj->C_t, traj->x, traj->k, m, delta_t);
        if (traj->g21 > xi)
        {
            traj->C_evo = 1;
        }
    }
}


void traj_step_IV(TrajContext* traj){
    double m = g_fssh_cfg.M_Quality;
    // 未发生态跃迁，直接跳过
    if (traj->C_evo == traj->C_evo_prev)
    {
        traj->F = F(traj->x, traj->C_evo);
        traj->C_evo_prev = traj->C_evo;
        traj->E_total = E_total(traj->k, m, traj->C_evo, traj->H_ad);
        return;
    }

    // 跃迁发生，修正动量
    double new_k = k_aft(traj->k, m, traj->H_ad, traj->C_evo);
    if (new_k == -10086.0)
    {
        // 动能不足，撤销本次态跃迁
        traj->C_evo = traj->C_evo_prev;
    }
    else
    {
        // 跃迁生效，更新动量，标记跃迁发生
        traj->k = new_k;
        traj->if_shift = true;
    }

    // 更新当前绝热态对应的核受力
    traj->F = F(traj->x, traj->C_evo);
    // 刷新上一步绝热态记录
    traj->C_evo_prev = traj->C_evo;
    // 更新总能量
    traj->E_total = E_total(traj->k, m, traj->C_evo, traj->H_ad);
}

// 单条轨迹完整演化直到停止
void traj_run_full(TrajContext* traj){
    while (true)
    {
        traj->step++;
        traj_step_II(traj);
        traj_step_III(traj);
        traj_step_IV(traj);

        if (traj_is_stop(traj))
        {
            break;
        }
    }
}

// 判断轨迹是否终止
bool traj_is_stop(TrajContext* traj){
    bool stop_flag = false;
    double x_max = g_fssh_cfg.X_Max;
    long step_max = g_fssh_cfg.Step_Max;

    if (traj->x > x_max || traj->x < -x_max)
    {
        stop_flag = true;
    }
    if (traj->step >= step_max)
    {
        stop_flag = true;
    }
    return stop_flag;
}

// 获取单条轨迹最终统计结果：[低能/高能, 透射/反射]   这里小技巧：列表传参，相当于指针传参
void traj_get_result(TrajContext* traj, int res[2]){
    // res[0]:0低能态  1高能态
    if (traj->C_evo == 1)
        res[0] = 0;
    else
        res[0] = 1;

    // res[1]:1透射(x>0)  0反射(x<0)
    if (traj->x > 0.0)
        res[1] = 1;
    else
        res[1] = 0;
}
