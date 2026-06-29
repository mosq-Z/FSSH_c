//本文件存放矩阵预算与矩阵定义
#ifndef MATRIX_MATH 
#define MATRIX_MATH 
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define HBAR 1.0

//二维矩阵结构
typedef struct Mat2c
{
    double complex V11,V12,V21,V22;
}Mat2c;

//二维向量结构
typedef struct Vec2c
{
    double complex E1,E2;
}Vec2c;

//矩阵乘矩阵
Mat2c matrix_multiplication(Mat2c* pM1,Mat2c* pM2);
//矩阵乘向量
Vec2c matrix_multiplied_by_vector(Mat2c* pM,Vec2c* pVc);
//矩阵转置共轭
Mat2c matrix_dagger(Mat2c* A);
//矩阵对角化
void matrix_diagonalization(Mat2c* pM,Mat2c* U,Vec2c* H_ad);



#endif
