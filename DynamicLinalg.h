#include<stdint.h>
#ifndef DYNAMICLINALG_H_DYNAMICLINALG_H
#define DYNAMICLINALG_H_DYNAMICLINALG_H

typedef struct Dynamic{
    double **mat;
    uint16_t row;
    uint16_t col;
}sDMAT_t;

sDMAT_t DynamicLinalg_create2D(uint16_t row, uint16_t col);
void DynamicLinalg_free2D(sDMAT_t *fp_sMat);
void DynamicLinalg_print2D(sDMAT_t *fp_sMat);
void DynamicLinalg_adds(sDMAT_t *fp_sMat, double fp_dNUM);
void DynamicLinalg_subs(sDMAT_t *fp_sMat, double fp_dNUM);
void DynamicLinalg_muls(sDMAT_t *fp_sMat, double fp_dNUM);
void DynamicLinalg_dvs(sDMAT_t *fp_sMat, double fp_dNUM);

sDMAT_t DynamicLinalg_addMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2);
sDMAT_t DynamicLinalg_subMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2);
sDMAT_t DynamicLinalg_mulMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2);
sDMAT_t DynamicLinalg_dot(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2);
sDMAT_t DynamicLinalg_eye(uint16_t fp_u16NUM);
void DynamicLinalg_transpose(sDMAT_t *fp_sMat);
double DynamicLinalg_trace(sDMAT_t *fp_sMat);

sDMAT_t DynamicLinalg_adj3x3(double fp_dMat[3][3]);
sDMAT_t DynamicLinalg_adj3x3mat(sDMAT_t *fp_sMat);
double DynamicLinalg_det3x3mat(sDMAT_t *fp_sMat);
#endif