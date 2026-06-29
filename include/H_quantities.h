//此文件对标math1库，用来存放函数

#ifndef H_QUANTITIES
#define H_QUANTITIES
#include <math.h>
#include "matrix_math.h"



// 同时获取矩阵酉变换矩阵和特征值
void get_eig_info(Mat2c* H, Mat2c* U, Vec2c* H_ad);

// 根据坐标x生成透热哈密顿矩阵
Mat2c H_di(double x);

// 获取哈密顿矩阵的绝热特征值
Vec2c H_ad(Mat2c* H);

// 获取特征向量构成的酉变换矩阵
Mat2c U(Mat2c* H);

// 中心差分计算透热哈密顿对x的导数
Mat2c dH_di_dx(double x);

// 计算两态体系非绝热耦合常数d12
double complex d12(double x);

// 计算当前电子态势能面对应的核受力
double F(double x, int C);

// 维尔莱特算法更新核坐标
double x_change_by_t(double k, double m, double x, double F, double del_t);

// 根据受力更新核动量
double k_change_by_t(double k, double F, double del_t);

// 计算粒子总能量：动能加绝热势能
double E_total(double k, double m, int C, Vec2c H_ad);

// 构建含时演化的两态有效哈密顿
Mat2c H_eff(Vec2c H_ad, double k, double m, double complex d12_val);

// RKG方法微分方程右端函数，计算电子振幅导数
Vec2c f_t_c(double k, double m, double x, double F, double f_del_t, Vec2c c);

// 使用RKG四阶方法推进电子振幅向量
Vec2c C_change_by_t_RKG(double k, double m, double x, double F, Mat2c H_eff_now, Vec2c C_t, double del_t);

// 计算电子态1向态2跃迁的跃迁概率
double g12(Vec2c C, double x, double k, double m, double del_t);

// 计算电子态2向态1跃迁的跃迁概率
double g21(Vec2c C, double x, double k, double m, double del_t);

// 电子态跃迁后更新核动量，能量不足返回-10086
double k_aft(double k, double m, Vec2c H_ad, int C);




#endif