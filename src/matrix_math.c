//本文件存放矩阵预算与矩阵定义
#include "matrix_math.h"
#include <complex.h>
#include <math.h>
#include <stdbool.h>

//矩阵乘矩阵
Mat2c matrix_multiplication(Mat2c* pM1,Mat2c* pM2){
    double complex V11 = pM1->V11;
    double complex V12 = pM1->V12;
    double complex V21 = pM1->V21;
    double complex V22 = pM1->V22;
    double complex V11_ = pM2->V11;
    double complex V12_ = pM2->V12;
    double complex V21_ = pM2->V21;
    double complex V22_ = pM2->V22;
    Mat2c ans = {V11*V11_+V12*V21_ , V11*V12_+V12*V22_ , V21*V11_+V22*V21_ , V21*V12_+V22*V22_};
    return ans;
}

//矩阵乘向量
Vec2c matrix_multiplied_by_vector(Mat2c* pM,Vec2c* pVc){
    double complex V11 = pM->V11;
    double complex V12 = pM->V12;
    double complex V21 = pM->V21;
    double complex V22 = pM->V22;
    double complex E1 = pVc->E1;
    double complex E2 = pVc->E2;
    Vec2c ans = {V11*E1+V12*E2 , V21*E1+V22*E2};
    return ans;
}

//矩阵转置共轭
Mat2c matrix_dagger(Mat2c* A){
    Mat2c res;
    res.V11 = conj(A->V11);
    res.V12 = conj(A->V21);
    res.V21 = conj(A->V12);
    res.V22 = conj(A->V22);
    return res;
}

//矩阵对角化
void matrix_diagonalization(Mat2c* pM,Mat2c* U,Vec2c* H_ad){
    double complex V11 = pM->V11;
    double complex V12 = pM->V12;
    double complex V21 = pM->V21;
    double complex V22 = pM->V22;
    double del = (V11+V22)*(V11+V22) - 4*(V11*V22-V21*V12);
    double sqrt_delta = sqrt(del);
    double tr = V11 + V22;
    double E1 = (tr - sqrt_delta) / 2.0;
    double E2 = (tr + sqrt_delta) / 2.0;

    H_ad->E1 = E1 + 0.0*I;
    H_ad->E2 = E2 + 0.0*I;

    double u00, u01, u10, u11;
    // 方程 (V11-E)x + V12 y = 0
    if (fabs(creal(V12)) > 1e-12)
    {
        u00 = V12;
        u10 = E1 - V11;
        u01 = V12;
        u11 = E2 - V11;
    }
    else
    {
        u00 = 1.0; u10 = 0.0;
        u01 = 0.0; u11 = 1.0;
    }

    // 归一化向量
    double norm0 = sqrt(u00*u00 + u10*u10);
    double norm1 = sqrt(u01*u01 + u11*u11);
    u00 /= norm0;
    u10 /= norm0;
    u01 /= norm1;
    u11 /= norm1;

    U->V11 = u00;
    U->V12 = u01;
    U->V21 = u10;
    U->V22 = u11;
}





