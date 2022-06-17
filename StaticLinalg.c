
#include "StaticLinalg.h"
#include<stdlib.h>
#include<stdio.h>



sSMAT_t StaticLinalg_zeros(uint16_t fp_u16ROW, uint16_t fp_u16COL){

    sSMAT_t fp_dMat;
    fp_dMat.row = fp_u16ROW;
    fp_dMat.col = fp_u16COL;
    for (uint16_t i = 0; i<fp_u16ROW; i++)
        for (uint16_t j = 0; j<fp_u16COL; j++)
            fp_dMat.mat[i][j] = 0;

    return fp_dMat;
}

void StaticLinalg_empty(sSMAT_t *fp_sMat){
    for(uint16_t i = 0; i<fp_sMat->row; i++){
        for(uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] = 0;
        }
    }
}

void StaticLinalg_print2D(sSMAT_t *fp_sMat){
    for(uint16_t i = 0; i<fp_sMat->row; i++){
        for(uint16_t j = 0; j<fp_sMat->col; j++){
            printf(" %f   ",fp_sMat->mat[i][j]);
        }
        printf("\n");
    }
}





/***************** Simple Array Operations *******************/
double StaticLinalg_det3x3(double fp_dMat[3][3]){
    double det = 0;
    det += fp_dMat[0][0]*((fp_dMat[1][1]*fp_dMat[2][2])-(fp_dMat[2][1]*fp_dMat[1][2]));
    det -= fp_dMat[0][1]*((fp_dMat[1][0]*fp_dMat[2][2])-(fp_dMat[2][0]*fp_dMat[1][2]));
    det += fp_dMat[0][2]*((fp_dMat[1][0]*fp_dMat[2][1])-(fp_dMat[2][0]*fp_dMat[1][1]));
    return det;
}
double StaticLinalg_trace3x3(double fp_dMat[3][3]){
    return fp_dMat[0][0] + fp_dMat[1][1] + fp_dMat[2][2];
}
sSMAT_t StaticLinalg_adj3x3(double fp_dMat[3][3]){
    double temp[3][3];
    sSMAT_t adj = StaticLinalg_zeros(3,3);

    temp[0][0] = ((fp_dMat[1][1]*fp_dMat[2][2])-(fp_dMat[2][1]*fp_dMat[1][2]));
    temp[0][1] = -((fp_dMat[1][0]*fp_dMat[2][2])-(fp_dMat[2][0]*fp_dMat[1][2]));
    temp[0][2] =  ((fp_dMat[1][0]*fp_dMat[2][1])-(fp_dMat[2][0]*fp_dMat[1][1]));

    temp[1][0] = -((fp_dMat[0][1]*fp_dMat[2][2])-(fp_dMat[2][1]*fp_dMat[0][2]));
    temp[1][1] = ((fp_dMat[0][0]*fp_dMat[2][2])-(fp_dMat[2][0]*fp_dMat[0][2]));
    temp[1][2] = -((fp_dMat[0][0]*fp_dMat[2][1])-(fp_dMat[2][0]*fp_dMat[0][1]));

    temp[2][0] = ((fp_dMat[0][1]*fp_dMat[1][2])-(fp_dMat[1][1]*fp_dMat[0][2]));
    temp[2][1] = -((fp_dMat[0][0]*fp_dMat[1][2])-(fp_dMat[1][0]*fp_dMat[0][2]));
    temp[2][2] = ((fp_dMat[0][0]*fp_dMat[1][1])-(fp_dMat[1][0]*fp_dMat[0][1]));

    for (uint16_t i = 0; i < 3; i++){
        for (uint16_t j = 0; j < 3; j++){
            adj.mat[i][j] = temp[i][j];
        }
    }
    return adj;
}



/*********** Scalar matrix operation function ***********/

void StaticLinalg_adds(sSMAT_t *fp_sMat, double fp_dNum){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] += fp_dNum;
        }
    }
}
void StaticLinalg_subs(sSMAT_t *fp_sMat, double fp_dNum){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] -= fp_dNum;
        }
    }
}
void StaticLinalg_muls(sSMAT_t *fp_sMat, double fp_dNum){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] *= fp_dNum;
        }
    }
}
void StaticLinalg_divs(sSMAT_t *fp_sMat, double fp_dNum){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] /= fp_dNum;
        }
    }
}
/***************** Matrix operation function ******************/
sSMAT_t StaticLinalg_addMat(sSMAT_t *mat1, sSMAT_t *mat2){
    sSMAT_t result;
    if ((mat1->row == mat2->row) && (mat1->col == mat2->col)){
        result = StaticLinalg_zeros(mat1->row, mat1->col);
        result.row = mat1->row;
        result.col = mat2->col;
        for (uint16_t i = 0; i<mat1->row; i++){
            for (uint16_t j = 0; j<mat1->col;j++){
                result.mat[i][j] = mat1->mat[i][j] + mat2->mat[i][j];
            }
        }
    }
    return result;
}

sSMAT_t StaticLinalg_subMat(sSMAT_t *mat1, sSMAT_t *mat2){
    sSMAT_t result;
    if ((mat1->row == mat2->row) && (mat1->col == mat2->col)){
        result = StaticLinalg_zeros(mat1->row, mat1->col);
        result.row = mat1->row;
        result.col = mat2->col;
        for (uint16_t i = 0; i<mat1->row; i++){
            for (uint16_t j = 0; j<mat1->col;j++){
                result.mat[i][j] = mat1->mat[i][j] - mat2->mat[i][j];
            }
        }
    }
    return result;
}
sSMAT_t StaticLinalg_mulMat(sSMAT_t *mat1, sSMAT_t *mat2){
    sSMAT_t result;
    if ((mat1->row == mat2->row) && (mat1->col == mat2->col)){
        result = StaticLinalg_zeros(mat1->row, mat1->col);
        result.row = mat1->row;
        result.col = mat2->col;
        for (uint16_t i = 0; i<mat1->row; i++){
            for (uint16_t j = 0; j<mat1->col;j++){
                result.mat[i][j] = mat1->mat[i][j] * mat2->mat[i][j];
            }
        }
    }
    return result;
}

sSMAT_t StaticLinalg_dot(sSMAT_t *mat1, sSMAT_t *mat2){
    sSMAT_t result;
    if ((mat1->col == mat2->row)){
        result = StaticLinalg_zeros(mat1->row, mat2->col);
        result.row = mat1->row;
        result.col = mat2->col;
        for (uint16_t i = 0; i < result.row; i++) {
            for (uint16_t j = 0; j < result.col; j++) {
                for (uint16_t k = 0; k < mat1->col; k++)
                    result.mat[i][j] += mat1->mat[i][k] * mat2->mat[k][j];
            }
        }
    }
    return result;

}
/******************Linear Algebra functions*******************************/

sSMAT_t StaticLinalg_eye(uint16_t fp_u16NUM){
    sSMAT_t identity = StaticLinalg_zeros(fp_u16NUM, fp_u16NUM);
    for (uint16_t i = 0; i < fp_u16NUM; i++)
        identity.mat[i][i] = 1;
    return identity;
}

void StaticLinalg_transpose(sSMAT_t *fp_dMat){
    sSMAT_t temp = StaticLinalg_zeros(fp_dMat->col, fp_dMat->row);

    for (uint16_t i = 0; i < fp_dMat->row; i++) {
        for (uint16_t j = 0; j < fp_dMat->col; j++) {
            temp.mat[j][i] = fp_dMat->mat[i][j];
        }
    }

    fp_dMat->row = temp.col;
    fp_dMat->col = temp.row;
    for (uint16_t i = 0; i < fp_dMat->row; i++) {
        for (uint16_t j = 0; j < fp_dMat->col; j++) {
             fp_dMat->mat[i][j] = temp.mat[i][j];
        }
    }

}

double StaticLinalg_trace(sSMAT_t *fp_dMat){
    double sum = 0;
    for (uint16_t i = 0; i < fp_dMat->row; ++i){
        sum += fp_dMat->mat[i][i];
    }
    return sum;
}

double StaticLinalg_trace3x3mat(sSMAT_t *fp_dMat){
    return fp_dMat->mat[0][0] + fp_dMat->mat[1][1] + fp_dMat->mat[2][2];
}


double StaticLinalg_det3x3mat(sSMAT_t *fp_dMat){
    double det = 0;
    det += fp_dMat->mat[0][0]*((fp_dMat->mat[1][1]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][1]*fp_dMat->mat[1][2]));
    det -= fp_dMat->mat[0][1]*((fp_dMat->mat[1][0]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][0]*fp_dMat->mat[1][2]));
    det += fp_dMat->mat[0][2]*((fp_dMat->mat[1][0]*fp_dMat->mat[2][1])-(fp_dMat->mat[2][0]*fp_dMat->mat[1][1]));
    return det;
}
sSMAT_t StaticLinalg_adj3x3mat(sSMAT_t *fp_dMat){
    double temp[3][3];
    sSMAT_t adj = StaticLinalg_zeros(3,3);

    temp[0][0] = ((fp_dMat->mat[1][1]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][1]*fp_dMat->mat[1][2]));
    temp[0][1] = -((fp_dMat->mat[1][0]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][0]*fp_dMat->mat[1][2]));
    temp[0][2] =  ((fp_dMat->mat[1][0]*fp_dMat->mat[2][1])-(fp_dMat->mat[2][0]*fp_dMat->mat[1][1]));

    temp[1][0] = -((fp_dMat->mat[0][1]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][1]*fp_dMat->mat[0][2]));
    temp[1][1] = ((fp_dMat->mat[0][0]*fp_dMat->mat[2][2])-(fp_dMat->mat[2][0]*fp_dMat->mat[0][2]));
    temp[1][2] = -((fp_dMat->mat[0][0]*fp_dMat->mat[2][1])-(fp_dMat->mat[2][0]*fp_dMat->mat[0][1]));

    temp[2][0] = ((fp_dMat->mat[0][1]*fp_dMat->mat[1][2])-(fp_dMat->mat[1][1]*fp_dMat->mat[0][2]));
    temp[2][1] = -((fp_dMat->mat[0][0]*fp_dMat->mat[1][2])-(fp_dMat->mat[1][0]*fp_dMat->mat[0][2]));
    temp[2][2] = ((fp_dMat->mat[0][0]*fp_dMat->mat[1][1])-(fp_dMat->mat[1][0]*fp_dMat->mat[0][1]));

    for (uint16_t i = 0; i < 3; i++){
        for (uint16_t j = 0; j < 3; j++){
            adj.mat[i][j] = temp[i][j];
        }
    }
    return adj;
}