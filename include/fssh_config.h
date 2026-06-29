// fssh_config.h 本文件存放类变量
#ifndef FSSH_CONFIG_H
#define FSSH_CONFIG_H

// 全局模拟超参数
typedef struct {
    long Step_Max;         // 单轨迹最大步数
    float X_Max;          // 坐标逃逸边界±X_Max
    int N_Single_K_Traj;   // 单个动量下模拟轨迹数
    float M_Quality;      // 粒子质量
    float Delta_T;        // 时间步长
} FSSH_GlobalConfig;

// 全局唯一配置实例
extern FSSH_GlobalConfig g_fssh_cfg;

// 配置修改接口（对应Python类Set方法）
void fssh_cfg_set_stepmax(long val);
void fssh_cfg_set_xmax(double val);
void fssh_cfg_set_ntraj(int val);
void fssh_cfg_set_mquality(double val);
void fssh_cfg_set_deltat(double val);
void fssh_cfg_reset_default(void);

#endif