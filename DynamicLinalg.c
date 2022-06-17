#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<DynamicLinalg.h>

sDMAT_t DynamicLinalg_create2D(uint16_t row, uint16_t col){
    sDMAT_t fp_sMat;
    fp_sMat.mat =NULL;
    if (row>=1 && col >=1) {
        fp_sMat.row = row;
        fp_sMat.col = col;
        fp_sMat.mat = (double **) calloc(fp_sMat.row, sizeof(double *));
        for (uint16_t i = 0; i < fp_sMat.row; i++) {
            fp_sMat.mat[i] = (double *) calloc(fp_sMat.col, sizeof(double));
        }
    }
    return fp_sMat;
}

void DynamicLinalg_free2D(sDMAT_t *fp_sMat){
    for (uint16_t i = 0; i<fp_sMat->row;i++){
        free(fp_sMat->mat[i]);
        fp_sMat->mat[i] = NULL;
    }
    free(fp_sMat->mat);
    fp_sMat->mat = NULL;
}

void DynamicLinalg_print2D(sDMAT_t *fp_sMat){
    for(uint16_t i = 0; i<fp_sMat->row; i++){
        for(uint16_t j = 0; j<fp_sMat->col; j++){
            printf(" %f   ",fp_sMat->mat[i][j]);
        }
        printf("\n");
    }
}


/****************************** Scalar matrix operation function **********************************/

void DynamicLinalg_adds(sDMAT_t *fp_sMat, double fp_dNUM){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] += fp_dNUM;
        }
    }
}
void DynamicLinalg_subs(sDMAT_t *fp_sMat, double fp_dNUM){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] -= fp_dNUM;
        }
    }
}
void DynamicLinalg_muls(sDMAT_t *fp_sMat, double fp_dNUM){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] *= fp_dNUM;
        }
    }
}
void DynamicLinalg_dvs(sDMAT_t *fp_sMat, double fp_dNUM){
    for (uint16_t i = 0; i<fp_sMat->row; i++){
        for (uint16_t j = 0; j<fp_sMat->col; j++){
            fp_sMat->mat[i][j] /= fp_dNUM;
        }
    }
}

/****************************** Matrix operation function **********************************/
sDMAT_t DynamicLinalg_addMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2){
    sDMAT_t result;
    result.mat=NULL;
    if ((fp_sMat1->row == fp_sMat2->row) && (fp_sMat1->col == fp_sMat2->col)){
        result = DynamicLinalg_create2D(fp_sMat1->row, fp_sMat1->col);
        result.row = fp_sMat1->row;
        result.col = fp_sMat2->col;
        for (uint16_t i = 0; i<fp_sMat1->row; i++){
            for (uint16_t j = 0; j<fp_sMat1->col;j++){
                result.mat[i][j] = fp_sMat1->mat[i][j] + fp_sMat2->mat[i][j];
            }
        }
    }
    return result;
}
sDMAT_t DynamicLinalg_subMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2){
    sDMAT_t result;
    result.mat=NULL;
    if ((fp_sMat1->row == fp_sMat2->row) && (fp_sMat1->col == fp_sMat2->col)){
        result = DynamicLinalg_create2D(fp_sMat1->row, fp_sMat1->col);
        result.row = fp_sMat1->row;
        result.col = fp_sMat2->col;
        for (uint16_t i = 0; i<fp_sMat1->row; i++){
            for (uint16_t j = 0; j<fp_sMat1->col;j++){
                result.mat[i][j] = fp_sMat1->mat[i][j] - fp_sMat2->mat[i][j];
            }
        }
    }
    return result;
}
sDMAT_t DynamicLinalg_mulMat(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2){
    sDMAT_t result;
    result.mat=NULL;
    if ((fp_sMat1->row == fp_sMat2->row) && (fp_sMat1->col == fp_sMat2->col)){
        result = DynamicLinalg_create2D(fp_sMat1->row, fp_sMat1->col);
        result.row = fp_sMat1->row;
        result.col = fp_sMat2->col;
        for (uint16_t i = 0; i<fp_sMat1->row; i++){
            for (uint16_t j = 0; j<fp_sMat1->col;j++){
                result.mat[i][j] = fp_sMat1->mat[i][j] * fp_sMat2->mat[i][j];
            }
        }
    }
    return result;

}
sDMAT_t DynamicLinalg_dot(sDMAT_t *fp_sMat1, sDMAT_t *fp_sMat2){
    sDMAT_t result;
    result.mat=NULL;
    if ((fp_sMat1->col == fp_sMat2->row)){
        result = DynamicLinalg_create2D(fp_sMat1->row, fp_sMat2->col);
        result.row = fp_sMat1->row;
        result.col = fp_sMat2->col;
        for (uint16_t i = 0; i < result.row; i++) {
            for (uint16_t j = 0; j < result.col; j++) {
                for (uint16_t k = 0; k < fp_sMat1->col; k++)
                    result.mat[i][j] += fp_sMat1->mat[i][k] * fp_sMat2->mat[k][j];
            }
        }
    }
    return result;

}
/******************************Linear Algebra functions******************************************/

sDMAT_t DynamicLinalg_eye(uint16_t fp_u16NUM){
    sDMAT_t identity;
    identity.mat =NULL;
    if (fp_u16NUM>0) {
        identity = DynamicLinalg_create2D(fp_u16NUM, fp_u16NUM);
        for (uint16_t i = 0; i < fp_u16NUM; i++) {
            identity.mat[i][i] = 1;
        }
    }
    return identity;
}

void DynamicLinalg_transpose(sDMAT_t *fp_sMat){
    double **temp =  (double**) calloc(fp_sMat->col, sizeof(double*));
    for(uint16_t i=0; i<fp_sMat->col; i++){
        temp[i] = (double*) calloc(fp_sMat->row, sizeof(double));
    }
    for (uint16_t i = 0; i < fp_sMat->row; ++i) {
        for (uint16_t j = 0; j < fp_sMat->col; ++j) {
            temp[j][i] = fp_sMat->mat[i][j];
        }
    }
    DynamicLinalg_free2D(fp_sMat);
    fp_sMat->mat = temp;
    uint16_t swap = fp_sMat->row;
    fp_sMat->row = fp_sMat->col;
    fp_sMat->col = swap;
}



double DynamicLinalg_trace(sDMAT_t *fp_sMat){
    double sum = 0;
    for (uint16_t i = 0; i < fp_sMat->row; ++i){
        sum += fp_sMat->mat[i][i];
    }
    return sum;
}

sDMAT_t DynamicLinalg_adj3x3(double fp_dMat[3][3]){
    double temp[3][3];
    sDMAT_t adj = DynamicLinalg_create2D(3,3);
    adj.row = 3;
    adj.col = 3;
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

sDMAT_t DynamicLinalg_adj3x3mat(sDMAT_t *fp_sMat){
    double temp[3][3];
    sDMAT_t adj = DynamicLinalg_create2D(3,3);
    adj.row = 3;
    adj.col = 3;
    temp[0][0] = ((fp_sMat->mat[1][1]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][1]*fp_sMat->mat[1][2]));
    temp[0][1] = -((fp_sMat->mat[1][0]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][0]*fp_sMat->mat[1][2]));
    temp[0][2] =  ((fp_sMat->mat[1][0]*fp_sMat->mat[2][1])-(fp_sMat->mat[2][0]*fp_sMat->mat[1][1]));

    temp[1][0] = -((fp_sMat->mat[0][1]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][1]*fp_sMat->mat[0][2]));
    temp[1][1] = ((fp_sMat->mat[0][0]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][0]*fp_sMat->mat[0][2]));
    temp[1][2] = -((fp_sMat->mat[0][0]*fp_sMat->mat[2][1])-(fp_sMat->mat[2][0]*fp_sMat->mat[0][1]));

    temp[2][0] = ((fp_sMat->mat[0][1]*fp_sMat->mat[1][2])-(fp_sMat->mat[1][1]*fp_sMat->mat[0][2]));
    temp[2][1] = -((fp_sMat->mat[0][0]*fp_sMat->mat[1][2])-(fp_sMat->mat[1][0]*fp_sMat->mat[0][2]));
    temp[2][2] = ((fp_sMat->mat[0][0]*fp_sMat->mat[1][1])-(fp_sMat->mat[1][0]*fp_sMat->mat[0][1]));

    for (uint16_t i = 0; i < 3; i++){
        for (uint16_t j = 0; j < 3; j++){
            adj.mat[i][j] = temp[i][j];
        }
    }
    return adj;
}
double DynamicLinalg_det3x3mat(sDMAT_t *fp_sMat){
    double det = 0;
    det += fp_sMat->mat[0][0]*((fp_sMat->mat[1][1]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][1]*fp_sMat->mat[1][2]));
    det -= fp_sMat->mat[0][1]*((fp_sMat->mat[1][0]*fp_sMat->mat[2][2])-(fp_sMat->mat[2][0]*fp_sMat->mat[1][2]));
    det += fp_sMat->mat[0][2]*((fp_sMat->mat[1][0]*fp_sMat->mat[2][1])-(fp_sMat->mat[2][0]*fp_sMat->mat[1][1]));
    return det;
}
