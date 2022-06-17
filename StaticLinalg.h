
#include<stdint.h>
#ifndef STATICLINALG_H_STATICLINALG_H
#define STATICLINALG_H_STATICLINALG_H


typedef struct Static{
    uint16_t row;
    uint16_t col;
    double mat[10][10];
}sSMAT_t;

sSMAT_t StaticLinalg_zeros(uint16_t fp_u16ROW, uint16_t fp_u16COL);
void StaticLinalg_empty(sSMAT_t *fp_sMat);
void StaticLinalg_print2D(sSMAT_t *fp_sMat);

double StaticLinalg_det3x3(double fp_dMat[3][3]);
double StaticLinalg_trace3x3(double fp_dMat[3][3]);
sSMAT_t StaticLinalg_adj3x3(double fp_dMat[3][3]);

// Scalar Matrix Operations
void StaticLinalg_adds(sSMAT_t *fp_sMat, double fp_dNum);
void StaticLinalg_subs(sSMAT_t *fp_sMat, double fp_dNum);
void StaticLinalg_muls(sSMAT_t *fp_sMat, double fp_dNum);
void StaticLinalg_divs(sSMAT_t *fp_sMat, double fp_dNum);

// Matrix Operations
sSMAT_t StaticLinalg_addMat(sSMAT_t *fp_sMat1, sSMAT_t *fp_sMat2);
sSMAT_t StaticLinalg_subMat(sSMAT_t *fp_sMat1, sSMAT_t *fp_sMat2);
sSMAT_t StaticLinalg_mulMat(sSMAT_t *fp_sMat1, sSMAT_t *fp_sMat2);
sSMAT_t StaticLinalg_dot(sSMAT_t *fp_sMat1, sSMAT_t *fp_sMat2);


// Linear Algebra Operations
sSMAT_t StaticLinalg_eye(uint16_t fp_u16NUM);
void StaticLinalg_transpose(sSMAT_t *fp_dMat);
double StaticLinalg_trace(sSMAT_t *fp_dMat);
double StaticLinalg_trace3x3mat(sSMAT_t *fp_dMat);
double StaticLinalg_det3x3mat(sSMAT_t *fp_dMat);
sSMAT_t StaticLinalg_adj3x3mat(sSMAT_t *fp_dMat);

#endif
