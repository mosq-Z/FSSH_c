#include "fssh_config.h"

// 全局配置变量，初始化赋值为你Python里的默认类参数
FSSH_GlobalConfig g_fssh_cfg = {
    .Step_Max = 400000L,
    .X_Max = 5.0,
    .N_Single_K_Traj = 4000,
    .M_Quality = 2000.0,
    .Delta_T = 1e-1
};

// 设置单轨迹最大演化步数
void fssh_cfg_set_stepmax(long val)
{
    if (val > 0)
    {
        g_fssh_cfg.Step_Max = val;
    }
}

// 设置坐标逃逸边界
void fssh_cfg_set_xmax(double val)
{
    if (val > 0.0)
    {
        g_fssh_cfg.X_Max = val;
    }
}

// 设置单个动量下模拟轨迹条数
void fssh_cfg_set_ntraj(int val)
{
    if (val > 0)
    {
        g_fssh_cfg.N_Single_K_Traj = val;
    }
}

// 设置粒子质量
void fssh_cfg_set_mquality(double val)
{
    if (val > 0.0)
    {
        g_fssh_cfg.M_Quality = val;
    }
}

// 设置动力学时间步长
void fssh_cfg_set_deltat(double val)
{
    if (val > 0.0)
    {
        g_fssh_cfg.Delta_T = val;
    }
}

// 重置所有参数为原始默认值
void fssh_cfg_reset_default(void)
{
    g_fssh_cfg.Step_Max = 400000L;
    g_fssh_cfg.X_Max = 5.0;
    g_fssh_cfg.N_Single_K_Traj = 4000;
    g_fssh_cfg.M_Quality = 2000.0;
    g_fssh_cfg.Delta_T = 1e-1;
}