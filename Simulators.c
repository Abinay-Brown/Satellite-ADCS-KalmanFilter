#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "Simulators.h"
#include "Frames.h"
#include "Physical_Constants.h"
#include "Controllers.h"
#include "IGRF.h"
#include "SunVector.h"

sSTATE1_t Simulators_DiffEqn1(sSTATE1_t state, sVEC_t TRQ, sVEC_t MOI){
    sSTATE1_t state_deriv;

    // Velocity
    state_deriv.POS.i = state.VEL.i;
    state_deriv.POS.j = state.VEL.j;
    state_deriv.POS.k = state.VEL.k;

    // Acceleration
    double pnorm = sqrt(pow(state.POS.i, 2) + pow(state.POS.j, 2) + pow(state.POS.k, 2));
    state_deriv.VEL.i =-(mu/pow(pnorm,3)) * state.POS.i;
    state_deriv.VEL.j =-(mu/pow(pnorm,3)) * state.POS.j;
    state_deriv.VEL.k =-(mu/pow(pnorm,3)) * state.POS.k;

    // MRP derivative
    double s = sqrt(pow(state.ATT.s1, 2) + pow(state.ATT.s2, 2) + pow(state.ATT.s3, 2));
    double s1 = state.ATT.s1;
    double s2 = state.ATT.s2;
    double s3 = state.ATT.s3;

    state_deriv.ATT.s1 = (1-pow(s,2)+(2*pow(s1,2))) * state.OMG.i;
    state_deriv.ATT.s1 += 2*((s1*s2)-s3) * state.OMG.j;
    state_deriv.ATT.s1 += 2*((s1*s3)+s2) * state.OMG.k;
    state_deriv.ATT.s1 /= 4;

    state_deriv.ATT.s2 = 2*((s2*s1)+s3) * state.OMG.i;
    state_deriv.ATT.s2 += (1-pow(s,2)+(2*pow(s2,2))) * state.OMG.j;
    state_deriv.ATT.s2 += 2*((s2*s3)-s1) * state.OMG.k;
    state_deriv.ATT.s2 /= 4;

    state_deriv.ATT.s3 = 2*((s3*s1)-s2) * state.OMG.i;
    state_deriv.ATT.s3 += 2*((s3*s2)+s1) * state.OMG.j;
    state_deriv.ATT.s3 += (1-pow(s,2)+(2*pow(s3,2))) * state.OMG.k;
    state_deriv.ATT.s3 /= 4;

    // Angular Acceleration Alpha
    // Input torque here, otherwise it is non accelerating.
    state_deriv.OMG.i = TRQ.i / MOI.i;
    state_deriv.OMG.j = TRQ.j / MOI.j;
    state_deriv.OMG.k = TRQ.k / MOI.k;

    return state_deriv;
}

void Simulators_Run1(){
    // Simulating B-dot Control
    FILE *fin;
    fin = fopen("C:\\Users\\abina\\Desktop\\Endurosat\\ADCS_MAIN\\IdealDetumbler_log.txt", "w");

    // Set initial conditions.
    sSTATE1_t state;
    state.POS.i = 4370.07; state.POS.j = -3776.848; state.POS.k = 3583.967;
    state.VEL.i = 5.74329; state.VEL.j = 2.42846; state.VEL.k = -4.44542;

    state.ATT.s1 = 0.4; state.ATT.s2 = 0.2; state.ATT.s3 = -0.1;
    state.OMG.i = -0.16; state.OMG.j = 0.07; state.OMG.k = -0.11;

    fprintf(fin,"%0.3f,%0.3f,%0.3f", state.POS.i, state.POS.j, state.POS.k);
    fprintf(fin,",%0.3f,%0.3f,%0.3f", state.VEL.i, state.VEL.j, state.VEL.k);
    fprintf(fin,",%0.3f,%0.3f,%0.3f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
    fprintf(fin,",%0.3f,%0.3f,%0.3f, %0.1f\n", state.OMG.i, state.OMG.j, state.OMG.k,0);

    double t0  = 0;         // Initial time.
    double h = 0.5;         // 0.5 seconds step size.
    double t1 = 300 * 60;   // Around 200 minutes.

    sVEC_t TRQ = Controllers_IdealDetumbler(state.ATT, state.OMG, 3000, 0.000005);
    TRQ.i *= -1;TRQ.j *= -1; TRQ.k *= -1;
    sVEC_t MOI; MOI.i = 0.0275; MOI.j = 0.0275; MOI.k = 0.0055;
    sSTATE1_t k1, k2, k3, k4, inp;
    // Runge Kutta-4 Method

    // Propagate.
    for (float t = h; t <= t1; t+=h){
        // integrating
        k1 = Simulators_DiffEqn1(state, TRQ, MOI);
        inp.POS.i = (state.POS.i + (h*k1.POS.i/2));
        inp.POS.j = (state.POS.j + (h*k1.POS.j/2));
        inp.POS.k = (state.POS.k + (h*k1.POS.k/2));

        inp.VEL.i = (state.VEL.i + (h*k1.VEL.i/2));
        inp.VEL.j = (state.VEL.j + (h*k1.VEL.j/2));
        inp.VEL.k = (state.VEL.k + (h*k1.VEL.k/2));

        inp.ATT.s1 = (state.ATT.s1 + (h*k1.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k1.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k1.ATT.s3/2));

        inp.OMG.i = (state.OMG.i + (h*k1.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k1.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k1.OMG.k/2));

        k2 = Simulators_DiffEqn1(inp, TRQ, MOI);

        inp.POS.i = (state.POS.i + (h*k2.POS.i/2));
        inp.POS.j = (state.POS.j + (h*k2.POS.j/2));
        inp.POS.k = (state.POS.k + (h*k2.POS.k/2));

        inp.VEL.i = (state.VEL.i + (h*k2.VEL.i/2));
        inp.VEL.j = (state.VEL.j + (h*k2.VEL.j/2));
        inp.VEL.k = (state.VEL.k + (h*k2.VEL.k/2));

        inp.ATT.s1 = (state.ATT.s1 + (h*k2.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k2.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k2.ATT.s3/2));

        inp.OMG.i = (state.OMG.i + (h*k2.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k2.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k2.OMG.k/2));

        k3 = Simulators_DiffEqn1(inp, TRQ, MOI);

        inp.POS.i = (state.POS.i + (h*k3.POS.i));
        inp.POS.j = (state.POS.j + (h*k3.POS.j));
        inp.POS.k = (state.POS.k + (h*k3.POS.k));

        inp.VEL.i = (state.VEL.i + (h*k3.VEL.i));
        inp.VEL.j = (state.VEL.j + (h*k3.VEL.j));
        inp.VEL.k = (state.VEL.k + (h*k3.VEL.k));

        inp.ATT.s1 = (state.ATT.s1 + (h*k3.ATT.s1));
        inp.ATT.s2 = (state.ATT.s2 + (h*k3.ATT.s2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k3.ATT.s3));

        inp.OMG.i = (state.OMG.i + (h*k3.OMG.i));
        inp.OMG.j = (state.OMG.j + (h*k3.OMG.j));
        inp.OMG.k = (state.OMG.k + (h*k3.OMG.k));

        k4 = Simulators_DiffEqn1(inp, TRQ, MOI);

        state.POS.i = state.POS.i + (h*(k1.POS.i+(2*k2.POS.i)+(2*k3.POS.i)+k4.POS.i)/6);
        state.POS.j = state.POS.j + (h*(k1.POS.j+(2*k2.POS.j)+(2*k3.POS.j)+k4.POS.j)/6);
        state.POS.k = state.POS.k + (h*(k1.POS.k+(2*k2.POS.k)+(2*k3.POS.k)+k4.POS.k)/6);

        state.VEL.i = state.VEL.i + (h*(k1.VEL.i+(2*k2.VEL.i)+(2*k3.VEL.i)+k4.VEL.i)/6);
        state.VEL.j = state.VEL.j + (h*(k1.VEL.j+(2*k2.VEL.j)+(2*k3.VEL.j)+k4.VEL.j)/6);
        state.VEL.k = state.VEL.k + (h*(k1.VEL.k+(2*k2.VEL.k)+(2*k3.VEL.k)+k4.VEL.k)/6);

        state.ATT.s1 = state.ATT.s1 + (h*(k1.ATT.s1+(2*k2.ATT.s1)+(2*k3.ATT.s1)+k4.ATT.s1)/6);
        state.ATT.s2 = state.ATT.s2 + (h*(k1.ATT.s2+(2*k2.ATT.s2)+(2*k3.ATT.s2)+k4.ATT.s2)/6);
        state.ATT.s3 = state.ATT.s3 + (h*(k1.ATT.s3+(2*k2.ATT.s3)+(2*k3.ATT.s3)+k4.ATT.s3)/6);

        state.OMG.i = state.OMG.i + (h*(k1.OMG.i+(2*k2.OMG.i)+(2*k3.OMG.i)+k4.OMG.i)/6);
        state.OMG.j = state.OMG.j + (h*(k1.OMG.j+(2*k2.OMG.j)+(2*k3.OMG.j)+k4.OMG.j)/6);
        state.OMG.k = state.OMG.k + (h*(k1.OMG.k+(2*k2.OMG.k)+(2*k3.OMG.k)+k4.OMG.k)/6);
        // State is Updated Now for the next time step.
        // Write to file. Or provide any inputs (torque).

        fprintf(fin,"%0.3f,%0.3f,%0.3f", state.POS.i, state.POS.j, state.POS.k);
        fprintf(fin,",%0.3f,%0.3f,%0.3f", state.VEL.i, state.VEL.j, state.VEL.k);
        fprintf(fin,",%0.3f,%0.3f,%0.3f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        fprintf(fin,",%0.3f,%0.3f,%0.3f, %0.1f\n", state.OMG.i, state.OMG.j, state.OMG.k,t);
        //printf(",%0.3f,%0.3f,%0.3f\n", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        printf(",%0.6f,%0.6f,%0.6f\n", state.OMG.i, state.OMG.j, state.OMG.k);
        TRQ = Controllers_IdealDetumbler(state.ATT, state.OMG, 3000, 0.0000005);
        TRQ.i *= -1;TRQ.j *= -1; TRQ.k *= -1;
    }
    fclose(fin);

}

void Simulators_Bdot(){
    // Simulating B-dot Control
    FILE *fin;
    FILE *lla;  // Read LLA Geocentric (Usually given by GPS)
    fin = fopen("C:\\Users\\abina\\Desktop\\Endurosat\\ADCS_MAIN\\Bdot_log.txt", "w");
    lla = fopen("C:\\Users\\abina\\Desktop\\Endurosat\\ADCS_MAIN\\ISS_ORBIT_LLA.txt", "r");

    // Set initial conditions.
    sSTATE1_t state;

    state.ATT.s1 = 0.4; state.ATT.s2 = 0.2; state.ATT.s3 = -0.1;
    state.OMG.i = 0.2; state.OMG.j = -0.18; state.OMG.k = -0.2;

    fprintf(fin,"%0.6f,%0.6f,%0.6f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
    fprintf(fin,",%0.6f,%0.6f,%0.6f\n", state.OMG.i, state.OMG.j, state.OMG.k);

    double h = 0.1;         // 0.1 seconds step size.
    double t1 = 900 * 60;
    double JDN, yy, mm, dd, hh, mn, sc;
    double lat, lon, alt;

    IGRF_init();
    sVEC_t TRQ; sVEC_t Sat; sBField_t B; sVEC_t Bold, Bnew; sVEC_t Kp; double Max; sDCM_t DCM;
    Kp. i = pow(10,9); Kp.j = pow(10,9); Kp.k = pow(10,9); // Setting Proportional gain to 30000
    sVEC_t MOI; MOI.i = 0.0275; MOI.j = 0.0275; MOI.k = 0.0275; // Moment of Inertia
    TRQ.i = 0;TRQ.j = 0; TRQ.k = 0;
    Max = 0.3;   // Am^2 Maximum Magnetic Moment

    sSTATE1_t k1, k2, k3, k4, inp;
    // Runge Kutta-4 Method

    // Propagate Attitude and Omega only.
    for (float t = 0; t <= t1; t+=h){
        // Get Initial Control Torque
        // Bfield new
        fscanf(lla,"%lf,%lf,%lf,%lf", &JDN, &lat, &lon, &alt);
        B = IGRF_computeMagneticField(alt+R_e, lat, lon, 13, JDN);
        DCM = Frames_MRPtoDCM(state.ATT);
        Bold = Frames_DCMxVEC(DCM, B.Bi);
        Bold.i += (rand()%200); Bold.j += (rand()%200); Bold.k += (rand()%200); // Add 200nT noise
        // integrating
        k1 = Simulators_DiffEqn1(state, TRQ, MOI);
        inp.ATT.s1 = (state.ATT.s1 + (h*k1.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k1.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k1.ATT.s3/2));
        inp.OMG.i = (state.OMG.i + (h*k1.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k1.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k1.OMG.k/2));
        k2 = Simulators_DiffEqn1(inp, TRQ, MOI);
        inp.ATT.s1 = (state.ATT.s1 + (h*k2.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k2.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k2.ATT.s3/2));
        inp.OMG.i = (state.OMG.i + (h*k2.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k2.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k2.OMG.k/2));
        k3 = Simulators_DiffEqn1(inp, TRQ, MOI);
        inp.ATT.s1 = (state.ATT.s1 + (h*k3.ATT.s1));
        inp.ATT.s2 = (state.ATT.s2 + (h*k3.ATT.s2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k3.ATT.s3));
        inp.OMG.i = (state.OMG.i + (h*k3.OMG.i));
        inp.OMG.j = (state.OMG.j + (h*k3.OMG.j));
        inp.OMG.k = (state.OMG.k + (h*k3.OMG.k));
        k4 = Simulators_DiffEqn1(inp, TRQ, MOI);
        state.ATT.s1 = state.ATT.s1 + (h*(k1.ATT.s1+(2*k2.ATT.s1)+(2*k3.ATT.s1)+k4.ATT.s1)/6);
        state.ATT.s2 = state.ATT.s2 + (h*(k1.ATT.s2+(2*k2.ATT.s2)+(2*k3.ATT.s2)+k4.ATT.s2)/6);
        state.ATT.s3 = state.ATT.s3 + (h*(k1.ATT.s3+(2*k2.ATT.s3)+(2*k3.ATT.s3)+k4.ATT.s3)/6);
        state.OMG.i = state.OMG.i + (h*(k1.OMG.i+(2*k2.OMG.i)+(2*k3.OMG.i)+k4.OMG.i)/6);
        state.OMG.j = state.OMG.j + (h*(k1.OMG.j+(2*k2.OMG.j)+(2*k3.OMG.j)+k4.OMG.j)/6);
        state.OMG.k = state.OMG.k + (h*(k1.OMG.k+(2*k2.OMG.k)+(2*k3.OMG.k)+k4.OMG.k)/6);
        // State is Updated Now for the next time step.
        // Write to file. Or provide any inputs (torque).
        // Bfield new
        fscanf(lla,"%lf,%lf,%lf,%lf", &JDN, &lat, &lon, &alt);
        B = IGRF_computeMagneticField(alt+R_e, lat, lon, 13, JDN);
        DCM = Frames_MRPtoDCM(state.ATT);
        Bnew = Frames_DCMxVEC(DCM, B.Bi);
        Bnew.i += (rand()%200); Bnew.j += (rand()%200); Bnew.k += (rand()%200); // Add 200nT noise
        Sat.i = lat; Sat.j = lon; Sat.k = alt + R_e;
        //printf(",%0.3f,%0.3f,%0.3f\n", Bnew.i, Bnew.j, Bnew.k);
        TRQ = Controllers_BdotControl(JDN,Sat,state.ATT,0.1, Bold, Bnew, Kp, Max);
        TRQ.i *= -1;TRQ.j *= -1; TRQ.k *= -1;

        fprintf(fin,"%0.6f,%0.6f,%0.6f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        fprintf(fin,",%0.6f,%0.6f,%0.6f\n", state.OMG.i, state.OMG.j, state.OMG.k);
        //printf(",%0.3f,%0.3f,%0.3f\n", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        printf(",%0.6f,%0.6f,%0.6f\n", state.OMG.i, state.OMG.j, state.OMG.k);
        //printf(",%0.6f,%0.6f,%0.6f\n", TRQ.i, TRQ.j, TRQ.k);

    }
    fclose(lla);

    fclose(fin);

}

void Simulators_Pointing(){
    FILE *fin, *lla;
    fin = fopen("C:\\Users\\abina\\Desktop\\Endurosat\\ADCS_MAIN\\PointingControl_log.txt", "w");
    lla = fopen("C:\\Users\\abina\\Desktop\\Endurosat\\ADCS_MAIN\\ISS_ORBIT_LLA.txt", "r");
    // Set initial conditions.
    sSTATE1_t state;
    state.POS.i = 4370.07; state.POS.j = -3776.848; state.POS.k = 3583.967;
    state.VEL.i = 5.74329; state.VEL.j = 2.42846; state.VEL.k = -4.44542;

    state.ATT.s1 = 0.4; state.ATT.s2 = 0.2; state.ATT.s3 = -0.1;
    state.OMG.i = 0; state.OMG.j = 0.0; state.OMG.k = 0;

    fprintf(fin,"%0.3f,%0.3f,%0.3f", state.POS.i, state.POS.j, state.POS.k);
    fprintf(fin,",%0.3f,%0.3f,%0.3f", state.VEL.i, state.VEL.j, state.VEL.k);
    fprintf(fin,",%0.3f,%0.3f,%0.3f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
    fprintf(fin,",%0.3f,%0.3f,%0.3f, %0.1f\n", state.OMG.i, state.OMG.j, state.OMG.k,0);

    double t0  = 0;         // Initial time.
    double h = 0.1;         // 0.5 seconds step size.
    double t1 = 900 * 60;   // Around 200 minutes.
    sVEC_t MOI; MOI.i = 0.0275; MOI.j = 0.0275; MOI.k = 0.0275;

    // Input Torque
    double JDN, lat, lon, alt;
    sVEC_t OPT1, OPT2; sPOINT_t OPT3;
    sGAIN_t gain;sVEC_t TRQ;
    double Max = 0.000005;
    sSTATE1_t k1, k2, k3, k4, inp;
    // Runge Kutta-4 Method

    // Propagate.
    for (float t = 0; t <= t1; t+=h){
        gain.K = 0.000001;
        gain.P.i = 0 ;
        gain.P.j = 0 ;
        gain.P.k = 0 ;

        Max = 0.00005;
        fscanf(lla,"%lf,%lf,%lf,%lf", &JDN, &lat, &lon, &alt);
        TRQ = Controllers_IdealPD(JDN, 2, state.ATT, state.OMG, gain, Max, OPT1, OPT2, OPT3);
        TRQ.i *= -1;TRQ.j *= -1; TRQ.k *= -1;

        // integrating
        k1 = Simulators_DiffEqn1(state, TRQ, MOI);
        inp.POS.i = (state.POS.i + (h*k1.POS.i/2));
        inp.POS.j = (state.POS.j + (h*k1.POS.j/2));
        inp.POS.k = (state.POS.k + (h*k1.POS.k/2));
        inp.VEL.i = (state.VEL.i + (h*k1.VEL.i/2));
        inp.VEL.j = (state.VEL.j + (h*k1.VEL.j/2));
        inp.VEL.k = (state.VEL.k + (h*k1.VEL.k/2));
        inp.ATT.s1 = (state.ATT.s1 + (h*k1.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k1.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k1.ATT.s3/2));
        inp.OMG.i = (state.OMG.i + (h*k1.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k1.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k1.OMG.k/2));
        k2 = Simulators_DiffEqn1(inp, TRQ, MOI);
        inp.POS.i = (state.POS.i + (h*k2.POS.i/2));
        inp.POS.j = (state.POS.j + (h*k2.POS.j/2));
        inp.POS.k = (state.POS.k + (h*k2.POS.k/2));
        inp.VEL.i = (state.VEL.i + (h*k2.VEL.i/2));
        inp.VEL.j = (state.VEL.j + (h*k2.VEL.j/2));
        inp.VEL.k = (state.VEL.k + (h*k2.VEL.k/2));
        inp.ATT.s1 = (state.ATT.s1 + (h*k2.ATT.s1/2));
        inp.ATT.s2 = (state.ATT.s2 + (h*k2.ATT.s2/2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k2.ATT.s3/2));
        inp.OMG.i = (state.OMG.i + (h*k2.OMG.i/2));
        inp.OMG.j = (state.OMG.j + (h*k2.OMG.j/2));
        inp.OMG.k = (state.OMG.k + (h*k2.OMG.k/2));
        k3 = Simulators_DiffEqn1(inp, TRQ, MOI);
        inp.POS.i = (state.POS.i + (h*k3.POS.i));
        inp.POS.j = (state.POS.j + (h*k3.POS.j));
        inp.POS.k = (state.POS.k + (h*k3.POS.k));
        inp.VEL.i = (state.VEL.i + (h*k3.VEL.i));
        inp.VEL.j = (state.VEL.j + (h*k3.VEL.j));
        inp.VEL.k = (state.VEL.k + (h*k3.VEL.k));
        inp.ATT.s1 = (state.ATT.s1 + (h*k3.ATT.s1));
        inp.ATT.s2 = (state.ATT.s2 + (h*k3.ATT.s2));
        inp.ATT.s3 = (state.ATT.s3 + (h*k3.ATT.s3));
        inp.OMG.i = (state.OMG.i + (h*k3.OMG.i));
        inp.OMG.j = (state.OMG.j + (h*k3.OMG.j));
        inp.OMG.k = (state.OMG.k + (h*k3.OMG.k));
        k4 = Simulators_DiffEqn1(inp, TRQ, MOI);
        state.POS.i = state.POS.i + (h*(k1.POS.i+(2*k2.POS.i)+(2*k3.POS.i)+k4.POS.i)/6);
        state.POS.j = state.POS.j + (h*(k1.POS.j+(2*k2.POS.j)+(2*k3.POS.j)+k4.POS.j)/6);
        state.POS.k = state.POS.k + (h*(k1.POS.k+(2*k2.POS.k)+(2*k3.POS.k)+k4.POS.k)/6);
        state.VEL.i = state.VEL.i + (h*(k1.VEL.i+(2*k2.VEL.i)+(2*k3.VEL.i)+k4.VEL.i)/6);
        state.VEL.j = state.VEL.j + (h*(k1.VEL.j+(2*k2.VEL.j)+(2*k3.VEL.j)+k4.VEL.j)/6);
        state.VEL.k = state.VEL.k + (h*(k1.VEL.k+(2*k2.VEL.k)+(2*k3.VEL.k)+k4.VEL.k)/6);
        state.ATT.s1 = state.ATT.s1 + (h*(k1.ATT.s1+(2*k2.ATT.s1)+(2*k3.ATT.s1)+k4.ATT.s1)/6);
        state.ATT.s2 = state.ATT.s2 + (h*(k1.ATT.s2+(2*k2.ATT.s2)+(2*k3.ATT.s2)+k4.ATT.s2)/6);
        state.ATT.s3 = state.ATT.s3 + (h*(k1.ATT.s3+(2*k2.ATT.s3)+(2*k3.ATT.s3)+k4.ATT.s3)/6);
        state.OMG.i = state.OMG.i + (h*(k1.OMG.i+(2*k2.OMG.i)+(2*k3.OMG.i)+k4.OMG.i)/6);
        state.OMG.j = state.OMG.j + (h*(k1.OMG.j+(2*k2.OMG.j)+(2*k3.OMG.j)+k4.OMG.j)/6);
        state.OMG.k = state.OMG.k + (h*(k1.OMG.k+(2*k2.OMG.k)+(2*k3.OMG.k)+k4.OMG.k)/6);
        // State is Updated Now for the next time step.
        // Write to file. Or provide any inputs (torque).

        //fprintf(fin,"%0.3f,%0.3f,%0.3f", state.POS.i, state.POS.j, state.POS.k);
        //fprintf(fin,",%0.3f,%0.3f,%0.3f", state.VEL.i, state.VEL.j, state.VEL.k);
        fprintf(fin,"%0.3f,%0.3f,%0.3f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        fprintf(fin,",%0.3f,%0.3f,%0.3f\n", state.OMG.i, state.OMG.j, state.OMG.k);
        //printf("%0.3f,%0.3f,%0.3f", state.ATT.s1, state.ATT.s2, state.ATT.s3);
        //printf(",%0.3f,%0.3f,%0.3f\n", state.OMG.i, state.OMG.j, state.OMG.k);
        sPOINT_t PNT = Frames_SunPointingDCM(JDN, TRQ, 0);
        sDCM_t DCM = Frames_MRPtoDCM(state.ATT);

        //printf("%0.9f,%0.9f,%0.9f,   ", TRQ.i, TRQ.j, TRQ.k);
        printf("%0.9f,%0.9f,%0.9f,   ", state.OMG.i, state.OMG.j, state.OMG.k);
        printf("%0.3f,%0.3f,%0.3f,   ", DCM.mat[0][0], DCM.mat[0][1], DCM.mat[0][2]);
        printf("%0.3f,%0.3f,%0.3f \n", PNT.DCM.mat[0][0], PNT.DCM.mat[0][1], PNT.DCM.mat[0][2]);
    }
    fclose(fin);
    fclose(lla);
}