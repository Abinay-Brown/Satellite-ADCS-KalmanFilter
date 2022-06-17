
#include "Frames.h"
#include<math.h>
#include<stdio.h>

// Global Variables

double delta_t;     // delta time
double sig_omega;   // Covariance of Rate Gyro
double sig_beta;    // Covariance of Bias vector

// The Variables below will change their values every iteration
double P_k[6][6];   // Error Covariance Matrix
double K[6][6];     // Kalman Gain
sVEC_t beta_k;      // Angular Velocity Bias vector
sQRN_t QRN_k;       // Attitude Quaternion (Corrected through filtering)
sVEC_t omega_k;     // Angular Velocity (Corrected through filtering)

void Kalman_Init(){

    // Initialize delta_t (Change Later)
    delta_t = 0.5;

    // Initialize the covariances. (Change Later)
    sig_omega = 0;
    sig_beta = 0;

    // Initialize Error covariance matrix. (Change Later)
    P_k[0][0] = 0;

    // Initialize Bias vector. (Change Later)
    beta_k.i = 0;

    // Initialize Quaternion. (Change Later)
    QRN_k.q1 = 0;

}
/**
 * 6-state Multiplicative Extended Kalman Filtering (MEKF)
 * Routine
 * @param   sVEC_t mi IGRF-13 Magnetic Field strength vector.
 * @param   sVEC_t si Sun vector.
 * @param   sVEC_t mb Magnetometer Body frame measurement.
 * @param   sVEC_t sb Sun Sensor body frame measurement.
 * @param   sVEC_t omega_m Angular Velocity measurement.
 * @param   sDCM_t A Attitude matrix to transform Inertial to Body frame.
 * @param   sig_mag standard deviation of magnetometer measurements
 * @param   sig_sun standard deviation og sun sensor measurements
 */
void Kalman_Filter(sVEC_t mi, sVEC_t si, sVEC_t mb, sVEC_t sb, sVEC_t omega_m, sDCM_t A, double sig_mag, double sig_sun){

    // 1. Prediction

    sVEC_t beta_k_1; // Angular Velocity bias vector
    sQRN_t QRN_k_1; // Quaternion

    // Calculating the Corrected Angular Vector
    omega_k.i = omega_m.i - beta_k.i;
    omega_k.j = omega_m.j - beta_k.j;
    omega_k.k = omega_m.k - beta_k.k;

    double omega_norm = sqrt(pow(omega_k.i,2) + pow(omega_k.j,2) + pow(omega_k.k,2));
    sVEC_t psi_k;
    psi_k.i = sin(0.5*omega_norm*delta_t)*omega_k.i/omega_norm;
    psi_k.j = sin(0.5*omega_norm*delta_t)*omega_k.j/omega_norm;
    psi_k.k = sin(0.5*omega_norm*delta_t)*omega_k.k/omega_norm;
    
    double cos_term = cos(0.5*omega_norm*delta_t);
    double TH[4][4];

    TH[0][0] = cos_term; TH[0][1] = psi_k.k; TH[0][2] = -psi_k.j; TH[0][3] = psi_k.i;
    TH[1][0] = -psi_k.k; TH[1][1] = cos_term; TH[1][2] = psi_k.i; TH[1][3] = psi_k.j;
    TH[2][0] = psi_k.j; TH[2][1] = -psi_k.i; TH[2][2] = cos_term; TH[2][3] = psi_k.k;
    TH[3][0] = -psi_k.i; TH[3][1] = -psi_k.j; TH[3][2] = -psi_k.k; TH[3][3] = cos_term;
    
    // Propagating Quaternion state and Bias Vector
    QRN_k_1.q1 = (TH[0][0]*QRN_k.q1)+(TH[0][1]*QRN_k.q2)+(TH[0][2]*QRN_k.q3)+(TH[0][3]*QRN_k.q0);
    QRN_k_1.q2 = (TH[1][0]*QRN_k.q1)+(TH[1][1]*QRN_k.q2)+(TH[1][2]*QRN_k.q3)+(TH[1][3]*QRN_k.q0);
    QRN_k_1.q3 = (TH[2][0]*QRN_k.q1)+(TH[2][1]*QRN_k.q2)+(TH[2][2]*QRN_k.q3)+(TH[2][3]*QRN_k.q0);
    QRN_k_1.q0 = (TH[3][0]*QRN_k.q1)+(TH[3][1]*QRN_k.q2)+(TH[3][2]*QRN_k.q3)+(TH[3][3]*QRN_k.q0);

    beta_k_1.i = beta_k.i;
    beta_k_1.j = beta_k.j;
    beta_k_1.k = beta_k.k;

    // Propagating the error covariance Matrix.
    double P_k_1[6][6]; // Error Covariance Matrix
    double PHI[6][6], PHI_T[6][6]; // State Transition Matrix
    double G[6][6], G_T[6][6]; // process Noise covariance Matrix
    double Q[6][6]; // Process Noise covariance Matrix

    double factor1 = sin(omega_norm*delta_t)/omega_norm;
    double factor2 = (1-cos(omega_norm*delta_t))/pow(omega_norm,2);
    double factor3 = factor2;
    double factor4 = ((omega_norm*delta_t)-sin(omega_norm*delta_t))/(pow(omega_norm, 3));
    cos_term = cos(omega_norm*delta_t);

    double S_omg[3][3];
    double S_omg_sq[3][3];

    // Constructing the S(omega_k) matrix.
    S_omg[0][0] = 0; S_omg[0][1] = -omega_k.k; S_omg[0][2] = omega_k.j;
    S_omg[1][0] = omega_k.k; S_omg[1][1] = 0; S_omg[1][2] = -omega_k.i;
    S_omg[2][0] = -omega_k.j; S_omg[2][1] = omega_k.i; S_omg[2][2] = 0;

    // Constructing the square of the S(omega_k) matrix.
    S_omg_sq[0][0] = -pow(omega_k.k,2)-pow(omega_k.j,2);
    S_omg_sq[0][1] = omega_k.j*omega_k.i;
    S_omg_sq[0][2] = omega_k.k*omega_k.i;

    S_omg_sq[1][0] = omega_k.j*omega_k.i;
    S_omg_sq[1][1] = -pow(omega_k.k,2)-pow(omega_k.i,2);
    S_omg_sq[1][2] = omega_k.k*omega_k.j;

    S_omg_sq[2][0] = omega_k.k*omega_k.i;
    S_omg_sq[2][1] = omega_k.k*omega_k.j;
    S_omg_sq[2][2] = -pow(omega_k.j,2)-pow(omega_k.i,2);

    // plugging in the phi11, phi12, phi21, phi22 into the PHI matrix directly.
    for (int i = 0; i<3; i++){
        for (int j = 0; j<3; j++){
            // Inserting phi11
            if(i == j){
                PHI[i][j] = 1 - (S_omg[i][j]*factor1) + (S_omg_sq[i][j]*factor2);
            }else{
                PHI[i][j] = 0 - (S_omg[i][j]*factor1) + (S_omg_sq[i][j]*factor2);
            }
            // Inserting phi12
            if(i == j){
                PHI[i][j+3] = (S_omg[i][j]*factor3) - delta_t - (S_omg_sq[i][j]*factor4);
            }else{
                PHI[i][j+3] = (S_omg[i][j]*factor3) - (S_omg_sq[i][j]*factor4);
            }
            // Inserting phi21
            PHI[i+3][j] = 0;
            // Inserting phi22
            if(i == j){
                PHI[i+3][j+3] = 1;
            }else{
                PHI[i+3][j+3] = 0;
            }
        }
    }
    // Finding PHI Transpose
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++) {
            PHI_T[j][i] = PHI[i][j];
        }
    }
    factor1 = (pow(sig_omega,2)*delta_t)+(pow(sig_beta,2)*pow(delta_t,3)/3);
    factor2 = -(0.5*pow(sig_beta,2)*pow(delta_t,2));
    factor4 = pow(sig_beta,2)*delta_t;

    // Constructing Q, G, and G transpose matrix.
    for (int i = 0; i<3; i++){
        for (int j = 0; j<3; j++){
            if (i==j){
                Q[i][j] = factor1;
                Q[i+3][j+3] = factor4;
                Q[i][j+3] = factor2;
                Q[i+3][j] = factor2;

                G[i][j] = -1;
                G[i+3][j+3] = 1;

                G_T[i][j] = -1;
                G_T[i+3][j+3] = 1;
            }
            else{
                Q[i][j] = 0;
                G[i][j] = 0;
                G_T[i][j] = 0;
            }
        }
    }
    // matrix multiplication for finding the propagated error Covariance matrix.

    double P_PHI_T[6][6];
    double Q_G_T[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            P_PHI_T[i][j] = 0;
            Q_G_T[i][j] = 0;
            for (int k = 0; k < 6; k++){
                P_PHI_T[i][j] += P_k[i][k] * PHI_T[k][j];
                Q_G_T[i][j] += Q[i][k] * G_T[k][j];
            }
        }
    }
    double PHI_P_PHI_T[6][6];
    double G_Q_G_T[6][6];

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            PHI_P_PHI_T[i][j] = 0;
            G_Q_G_T[i][j] = 0;
            for (int k = 0; k < 6; k++){
                PHI_P_PHI_T[i][j] += PHI[i][k] * P_PHI_T[k][j];
                G_Q_G_T[i][j] += G[i][k] * Q_G_T[k][j];
            }
            P_k_1[i][j] = PHI_P_PHI_T[i][j] + G_Q_G_T[i][j];
        }
    }
    ///////////////////////////////////////////////////////////////////
    //2. Measurement Update.
    double H[6][6], H_T[6][6];
    double R[6][6];
    double arr[6][6];
    // Constructing the H matrix and its transpose.
    H[0][0] = 0; H[0][1] = -2*mi.k; H[0][2] = 2*mi.j;
    H[1][0] = 2*mi.k; H[1][1] = 0; H[1][2] = -2*mi.i;
    H[2][0] = -2*mi.j; H[2][1] = 2*mi.i; H[2][2] = 0;

    H[3][0] = 0; H[3][1] = -2*si.k; H[3][2] = 2*si.j;
    H[4][0] = 2*si.k; H[4][1] = 0; H[4][2] = -2*si.i;
    H[5][0] = -2*mi.j; H[5][1] = 2*si.i; H[5][2] = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 3; j < 6; j++) {
            H[i][j] = 0;
        }
    }
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++) {
            H_T[j][i] = H[i][j];
        }
    }
    // Constructing the Measurement Covariance matrix.
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++) {
            if ((i == j) && (i<3)){
                R[i][j] = sig_mag*sig_mag;
            }
            else if ((i == j) && (i>=3)){
                R[i][j] = sig_sun*sig_sun;
            }
            else{
                R[i][j] = 0;
            }
        }
    }

    // Calculating the Kalman Gain Matrix.
    double P_H_T[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            P_H_T[i][j] = 0;
            for (int k = 0; k < 6; k++){
                P_H_T[i][j] += P_k_1[i][k] * H_T[k][j];
            }
        }
    }
    double H_P_H_T[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            H_P_H_T[i][j] = 0;
            for (int k = 0; k < 6; k++){
                H_P_H_T[i][j] += H[i][k] * P_H_T[k][j];
            }
            arr[i][j] = H_P_H_T[i][j] + R[i][j];
        }
    }
    // Finding inverse part of the K gain matrix
    double I[6][6];
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++){
            if (i==j){
                I[i][j] = 1;
            }
        }
    }
    // Gauss elimination
    double c;
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            if (i > j) {
                c = arr[i][j] / arr[j][j];
                for (int k = 0; k < 6; k++) {
                    arr[i][k] -= (c * arr[j][k]);
                    I[i][k] -= (c * I[j][k]);
                }
            }
        }
    }


    for (int j = 0; j < 6; j++) {
        c = arr[j][j];
        for (int k = 0; k < 6; k++) {
            arr[j][k] /= c;
            I[j][k] /= c;
        }
    }

    for (int j = 5; j >= 0; j--) {
        for (int i = 5; i >= 0; i--) {
            if (i<j){
                c = arr[i][j] / arr[j][j];
                for (int k = 5; k >= 0; k--) {
                    arr[i][k] -=(c * arr[j][k]);
                    I[i][k] -= (c * I[j][k]);
                }
            }
        }
    }
    // Now the Identity I matrix contains the inverse
    double H_T_I[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            H_T_I[i][j] = 0;
            for (int k = 0; k < 6; k++){
                H_T_I[i][j] += H_T[i][k] * I[k][j];
            }
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            K[i][j] = 0;
            for (int k = 0; k < 6; k++){
                K[i][j] += P_k_1[i][k] * H_T_I[k][j];
            }
        }
    }

    // Update the Error Covariance Matrix
    double I_K_H[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            I_K_H[i][j] = 0;
            for (int k = 0; k < 6; k++){
                I_K_H[i][j] += K[i][k] * H[k][j];
            }
            if (i==j){
                I_K_H[i][j] = 1-I_K_H[i][j]; // The Identity matrix part.
            }
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            P_k[i][j] = 0;
            for (int k = 0; k < 6; k++){
                P_k[i][j] += I_K_H[i][k] * P_k_1[k][j];
            }
        }
    }

    // Finding Corrected State Vector.

    // Calculate Estimated Vectors (matrix multiply Attitude matrix and inertial frame vectors).
    sVEC_t m_est;
    sVEC_t s_est;
    sVEC_t m_err, s_err;
    m_est.i = (A.mat[0][0]*mi.i) + (A.mat[0][1]*mi.j) + (A.mat[0][2]*mi.k);
    m_est.j = (A.mat[1][0]*mi.i) + (A.mat[1][1]*mi.j) + (A.mat[1][2]*mi.k);
    m_est.k = (A.mat[2][0]*mi.i) + (A.mat[2][1]*mi.j) + (A.mat[2][2]*mi.k);

    s_est.i = (A.mat[0][0]*si.i) + (A.mat[0][1]*si.j) + (A.mat[0][2]*si.k);
    s_est.j = (A.mat[1][0]*si.i) + (A.mat[1][1]*si.j) + (A.mat[1][2]*si.k);
    s_est.k = (A.mat[2][0]*si.i) + (A.mat[2][1]*si.j) + (A.mat[2][2]*si.k);

    m_err.i = mb.i - m_est.i;
    m_err.j = mb.j - m_est.j;
    m_err.k = mb.k - m_est.k;

    s_err.i = sb.i - s_est.i;
    s_err.j = sb.j - s_est.j;
    s_err.k = sb.k - s_est.k;

    sQRN_t qerr;
    sVEC_t berr;
    qerr.q1 =(K[0][0]*m_err.i)+(K[0][1]*m_err.j)+(K[0][2]*m_err.k)+(K[0][3]*s_err.i)+(K[0][4]*s_err.j)+(K[0][5]*s_err.k);
    qerr.q2 =(K[1][0]*m_err.i)+(K[1][1]*m_err.j)+(K[1][2]*m_err.k)+(K[1][3]*s_err.i)+(K[1][4]*s_err.j)+(K[1][5]*s_err.k);
    qerr.q3 =(K[2][0]*m_err.i)+(K[2][1]*m_err.j)+(K[2][2]*m_err.k)+(K[2][3]*s_err.i)+(K[2][4]*s_err.j)+(K[2][5]*s_err.k);

    berr.i =(K[3][0]*m_err.i)+(K[3][1]*m_err.j)+(K[3][2]*m_err.k)+(K[3][3]*s_err.i)+(K[3][4]*s_err.j)+(K[3][5]*s_err.k);
    berr.j =(K[4][0]*m_err.i)+(K[4][1]*m_err.j)+(K[4][2]*m_err.k)+(K[4][3]*s_err.i)+(K[4][4]*s_err.j)+(K[4][5]*s_err.k);
    berr.k =(K[5][0]*m_err.i)+(K[5][1]*m_err.j)+(K[5][2]*m_err.k)+(K[5][3]*s_err.i)+(K[5][4]*s_err.j)+(K[5][5]*s_err.k);

    qerr.q0 = sqrt(1 - ((qerr.q1*qerr.q1) + (qerr.q2*qerr.q2) + (qerr.q3*qerr.q3)));

    // Performing Quaternion cross product.
    double q_cross[4][4];
    q_cross[0][0] = qerr.q0; q_cross[0][1] = qerr.q3; q_cross[0][2] = -qerr.q2; q_cross[0][3] = qerr.q1;
    q_cross[1][0] = -qerr.q3; q_cross[1][1] = qerr.q0; q_cross[1][2] = qerr.q1; q_cross[1][3] = qerr.q2;
    q_cross[2][0] = qerr.q2; q_cross[2][1] = -qerr.q1; q_cross[2][2] = qerr.q0; q_cross[2][3] = qerr.q3;
    q_cross[3][0] = -qerr.q1; q_cross[3][1] = -qerr.q2; q_cross[3][2] = -qerr.q3; q_cross[3][3] = qerr.q0;

    // Final Corrected Quaternion
    QRN_k.q1 = (q_cross[0][0]*QRN_k_1.q1) + (q_cross[0][1]*QRN_k_1.q2) + (q_cross[0][2]*QRN_k_1.q3) + (q_cross[0][3]*QRN_k_1.q0);
    QRN_k.q2 = (q_cross[1][0]*QRN_k_1.q1) + (q_cross[1][1]*QRN_k_1.q2) + (q_cross[1][2]*QRN_k_1.q3) + (q_cross[1][3]*QRN_k_1.q0);
    QRN_k.q3 = (q_cross[2][0]*QRN_k_1.q1) + (q_cross[2][1]*QRN_k_1.q2) + (q_cross[2][2]*QRN_k_1.q3) + (q_cross[2][3]*QRN_k_1.q0);
    QRN_k.q0 = (q_cross[3][0]*QRN_k_1.q1) + (q_cross[3][1]*QRN_k_1.q2) + (q_cross[3][2]*QRN_k_1.q3) + (q_cross[3][3]*QRN_k_1.q0);

    // Final Corrected Rate-Gyro Bias Vector
    beta_k.i = beta_k_1.i + berr.i;
    beta_k.j = beta_k_1.j + berr.j;
    beta_k.k = beta_k_1.k + berr.k;
}

void inv6x6(){
    double arr[6][6];
    arr[0][0]=56;arr[0][1]=76;arr[0][2]=19;arr[0][3]=93;arr[0][4]=30;arr[0][5]=87;
    arr[1][0]=45;arr[1][1]=34;arr[1][2]=82;arr[1][3]=34;arr[1][4]=69;arr[1][5]=57;
    arr[2][0]=32;arr[2][1]=92;arr[2][2]=23;arr[2][3]=56;arr[2][4]=49;arr[2][5]=64;
    arr[3][0]=34;arr[3][1]=23;arr[3][2]=12;arr[3][3]=73;arr[3][4]=56;arr[3][5]=17;
    arr[4][0]=12;arr[4][1]=37;arr[4][2]=65;arr[4][3]=45;arr[4][4]=67;arr[4][5]=47;
    arr[5][0]=91;arr[5][1]=61;arr[5][2]=28;arr[5][3]=76;arr[5][4]=12;arr[5][5]=23;


    double I[6][6];
    for (int i = 0; i<6; i++){
        for (int j = 0; j<6; j++){
            if (i==j){
                I[i][j] = 1;
            }
        }
    }
    // Gauss elimination
    double c;
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            if (i > j) {
                c = arr[i][j] / arr[j][j];
                for (int k = 0; k < 6; k++) {
                    arr[i][k] -= (c * arr[j][k]);
                    I[i][k] -= (c * I[j][k]);
                }
            }
        }
    }


    for (int j = 0; j < 6; j++) {
        c = arr[j][j];
        for (int k = 0; k < 6; k++) {
            arr[j][k] /= c;
            I[j][k] /= c;
        }
    }

    for (int j = 5; j >= 0; j--) {
        for (int i = 5; i >= 0; i--) {
            if (i<j){
                c = arr[i][j] / arr[j][j];
                for (int k = 5; k >= 0; k--) {
                    arr[i][k] -=(c * arr[j][k]);
                    I[i][k] -= (c * I[j][k]);
                }
            }
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            printf(" %0.4f ",I[i][j]);
        }
        printf("\n");
    }

}

void Magnetometer_measure(){

}
void SunSensor_measure(){

}
void rateGyro_measure(){

}