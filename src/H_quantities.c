//此文件对标math1库，用来存放函数

#include "H_quantities.h"



//透热哈密顿量
Mat2c H_di(double x)
{
    Mat2c H;


    // ========== 模型一 ==========
    const double A = 0.01;
    const double B = 1.6;
    const double C = 0.005;
    const double D = 1.0;

    double V11_real;
    if (x > 0.0)
        V11_real = A * (1.0 - exp(-B * x));
    else
        V11_real = -A * (1.0 - exp(B * x));

    H.V11 = V11_real + 0.0*I;
    H.V22 = -V11_real + 0.0*I;
    double off = C * exp(-D * x * x);
    H.V12 = off + 0.0*I;
    H.V21 = H.V12;


    // ========== 模型二  ==========
/*     const double A = 0.1;
    const double B = 0.28;
    const double E0 = 0.05;
    const double C = 0.015;
    const double D = 0.06;

    H.V11 = 0.0 + 0.0*I;
    double V22_real = -A * exp(-B * x * x) + E0;
    H.V22 = V22_real + 0.0*I;
    double off = C * exp(-D * x * x);
    H.V12 = off + 0.0*I;
    H.V21 = H.V12;
 */

/*     // ========== 模型三 ==========
    const double A = 6.0e-4;
    const double B = 0.1;
    const double C = 0.9;

    H.V11 = A + 0.0 * I;
    H.V22 = -A + 0.0 * I;

    double V12_real;
    if (x < 0.0)
    {
        V12_real = B * exp(C * x);
    }
    else
    {
        V12_real = B * (2.0 - exp(-C * x));
    }
    // 非对角元为实数，虚部0
    H.V12 = V12_real + 0.0 * I;
    H.V21 = H.V12;
 */
    return H;
}

// 一次性获取酉变换矩阵U 和 绝热哈密顿特征值
void get_eig_info(Mat2c* H, Mat2c* U, Vec2c* H_ad)
{
    matrix_diagonalization(H, U, H_ad);
}

// 获取哈密顿矩阵的绝热特征值（绝热哈密顿对角元）
Vec2c H_ad(Mat2c* H)
{
    Mat2c U_tmp;
    Vec2c eig;
    matrix_diagonalization(H, &U_tmp, &eig);
    return eig;
}

// 获取酉变换矩阵U（特征向量矩阵，按列存放）
Mat2c U(Mat2c* H)
{
    Mat2c U_mat;
    Vec2c eig_tmp;
    matrix_diagonalization(H, &U_mat, &eig_tmp);
    return U_mat;
}

// 1. 中心差分计算透热哈密顿对x的导数 dH/dx
Mat2c dH_di_dx(double x)
{
    float DX_FD = 1e-7;
    Mat2c H_plus  = H_di(x + DX_FD);
    Mat2c H_minus = H_di(x - DX_FD);
    Mat2c dH;
    double fac = 1.0 / (2.0 * DX_FD);

    dH.V11 = (H_plus.V11 - H_minus.V11) * fac;
    dH.V12 = (H_plus.V12 - H_minus.V12) * fac;
    dH.V21 = (H_plus.V21 - H_minus.V21) * fac;
    dH.V22 = (H_plus.V22 - H_minus.V22) * fac;

    return dH;
}

// 2. 计算非绝热耦合 d12 = [U†·dH/dx·U]_{0,1} / (E2-E1)
double complex d12(double x)
{
    Mat2c H = H_di(x);
    Mat2c U;
    Vec2c eig;
    get_eig_info(&H, &U, &eig);

    Mat2c dHdx = dH_di_dx(x);
    Mat2c U_dagger = matrix_dagger(&U);

    // U† * dHdx
    Mat2c tmp1 = matrix_multiplication(&U_dagger, &dHdx);
    // U† * dHdx * U
    Mat2c mat_final = matrix_multiplication(&tmp1, &U);

    double complex E1 = eig.E1;
    double complex E2 = eig.E2;
    double complex denom = E2 - E1;

    return mat_final.V12 / denom;
}

// 3. 计算核受力 F = - dE_C / dx
double F(double x, int C)
{
    double  DX_FD = 1e-8;
    Mat2c H_di_plus = H_di(x + DX_FD);
    Mat2c H_di_minus = H_di(x - DX_FD);
    Vec2c eig_plus  = H_ad(&H_di_plus);
    Vec2c eig_minus = H_ad(&H_di_minus);
    double fac = 1.0 / (2.0 * DX_FD);

    double e_p, e_m;
    if(C == 1)
    {
        e_p = creal(eig_plus.E1);
        e_m = creal(eig_minus.E1);
    }
    else
    {
        e_p = creal(eig_plus.E2);
        e_m = creal(eig_minus.E2);
    }
    double dEdx = (e_p - e_m) * fac;
    return -dEdx;
}

// 4. 维尔莱特位置更新 x_new = x + k*dt/m + 0.5*F/m*dt²
double x_change_by_t(double k, double m, double x, double F, double del_t)
{
    double term1 = k * del_t / m;
    double term2 = 0.5 * (F / m) * del_t * del_t;
    return x + term1 + term2;
}

// 5. 动量更新 k_new = k + F * del_t
double k_change_by_t(double k, double F, double del_t)
{
    return k + F * del_t;
}

// 6. 总能量：动能 + 当前绝热势能面能量
double E_total(double k, double m, int C, Vec2c H_ad)
{
    double kin = k * k / (2.0 * m);
    double pot;
    if(C == 1)
        pot = creal(H_ad.E1);
    else
        pot = creal(H_ad.E2);
    return kin + pot;
}

// 7. 构造有效哈密顿 H_eff
Mat2c H_eff(Vec2c H_ad, double k, double m, double complex d12_val)
{
    Mat2c Heff;
    double v = k / m;
    double complex E1 = H_ad.E1;
    double complex E2 = H_ad.E2;

    double complex imag_unit = 1.0 * I;
    double complex off = -v * d12_val * imag_unit * HBAR;

    Heff.V11 = E1;
    Heff.V22 = E2;
    Heff.V12 = off;
    Heff.V21 = -off;

    return Heff;
}


// 1. RKG方法的ODE右端函数 dc/dt = (-i/ħ)·H_eff · c
Vec2c f_t_c(double k, double m, double x, double F, double f_del_t, Vec2c c)
{
    // 先推进坐标
    double x_f = x_change_by_t(k, m, x, F, f_del_t);
    Mat2c H_di_f = H_di(x_f);
    Mat2c U_f;
    Vec2c H_ad_f;
    get_eig_info(&H_di_f, &U_f, &H_ad_f);
    double complex d12_f = d12(x_f);
    // 推进动量
    double k_f = k_change_by_t(k, F, f_del_t);
    // 构造当前有效哈密顿
    Mat2c H_eff_f = H_eff(H_ad_f, k_f, m, d12_f);
    // H_eff 左乘振幅向量
    Vec2c Hc = matrix_multiplied_by_vector(&H_eff_f, &c);
    double complex coeff = (-1.0 * I) / HBAR;
    Vec2c ans = {coeff * Hc.E1, coeff * Hc.E2};
    return ans;
}

// 2. RKG积分推进电子振幅向量
Vec2c C_change_by_t_RKG(double k, double m, double x, double F, Mat2c H_eff_now, Vec2c C_t, double del_t)
{
    // k1 = dt * (-i/ħ) * H_eff_now @ C_t
    Vec2c H_ct = matrix_multiplied_by_vector(&H_eff_now, &C_t);
    double complex coeff = (-1.0 * I) / HBAR;
    Vec2c k1 = {del_t * coeff * H_ct.E1, del_t * coeff * H_ct.E2};

    // C_t + k1/2
    Vec2c C_k1_half = {C_t.E1 + k1.E1 / 2.0, C_t.E2 + k1.E2 / 2.0};
    Vec2c k2 = f_t_c(k, m, x, F, del_t / 2.0, C_k1_half);
    k2.E1 *= del_t;
    k2.E2 *= del_t;
    
    double SQRT2 = 1.4142135623730950488016887242097;

    // C_t + k1/2 + (sqrt(2)-1)*k2
    double alpha = SQRT2 - 1.0;
    Vec2c C_k3 = {
        C_t.E1 + k1.E1/2.0 + alpha * k2.E1,
        C_t.E2 + k1.E2/2.0 + alpha * k2.E2
    };
    Vec2c k3 = f_t_c(k, m, x, F, del_t / 2.0, C_k3);
    k3.E1 *= del_t;
    k3.E2 *= del_t;

    // C_t + (1-sqrt2/2)*k2 + (sqrt2/2)*k3
    double beta1 = 1.0 - SQRT2 / 2.0;
    double beta2 = SQRT2 / 2.0;
    Vec2c C_k4 = {
        C_t.E1 + beta1 * k2.E1 + beta2 * k3.E1,
        C_t.E2 + beta1 * k2.E2 + beta2 * k3.E2
    };
    Vec2c k4 = f_t_c(k, m, x, F, del_t, C_k4);
    k4.E1 *= del_t;
    k4.E2 *= del_t;

    // 组合系数
    double c2 = 2.0 - SQRT2;
    double c3 = 2.0 + SQRT2;
    double fac = 1.0 / 6.0;

    Vec2c ans;
    ans.E1 = C_t.E1 + fac * (k1.E1 + c2 * k2.E1 + c3 * k3.E1 + k4.E1);
    ans.E2 = C_t.E2 + fac * (k1.E2 + c2 * k2.E2 + c3 * k3.E2 + k4.E2);
    return ans;
}

// 3. 跃迁概率 g12：态1 → 态2
double g12(Vec2c C, double x, double k, double m, double del_t)
{
    double complex d12_val = d12(x);
    double complex d21_conj = conj(-d12_val);
    double complex c1 = C.E1;
    double complex c2 = C.E2;
    double v = k / m;

    double complex term = c1 * conj(c2) * v * d21_conj;
    double real_part = creal(term);
    double abs_c1_sq = cabs(c1) * cabs(c1);

    return del_t * (-2.0 * real_part) / abs_c1_sq;
}

// 4. 跃迁概率 g21：态2 → 态1
double g21(Vec2c C, double x, double k, double m, double del_t)
{
    double complex d12_val = d12(x);
    double complex c1 = C.E1;
    double complex c2 = C.E2;
    double v = k / m;

    double complex term = c2 * conj(c1) * v * d12_val;
    double real_part = creal(term);
    double abs_c2_sq = cabs(c2) * cabs(c2);

    return del_t * (-2.0 * real_part) / abs_c2_sq;
}

// 5. 跃迁后更新动量，能量不足返回 -10086
double k_aft(double k, double m, Vec2c H_ad, int C)
{
    double E1 = creal(H_ad.E1);
    double E2 = creal(H_ad.E2);
    double kin_old = k * k / (2.0 * m);
    double E_new;

    if (C == 0)
    {
        // 从态0跃迁到态1，需要吸收能量 E2-E1
        E_new = kin_old - (E2 - E1);
    }
    else
    {
        // 从态1跃迁到态0，释放能量 E2-E1
        E_new = kin_old + (E2 - E1);
    }

    if (E_new >= 0.0)
    {
        double k_abs = sqrt(2.0 * E_new * m);
        // 保持原来动量符号
        if (k >= 0)
            return k_abs;
        else
            return -k_abs;
    }
    else
    {
        // 能量不足，禁止跃迁
        return -10086.0;
    }
}