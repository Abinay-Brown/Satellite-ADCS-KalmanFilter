#include "Frames.h"
#ifndef SIMULATORS_H_SIMULATORS_H
#define SIMULATORS_H_SIMULATORS_H

typedef struct state1{
    // Inertial Frame components
    sVEC_t POS; // in km
    sVEC_t VEL; // in km/s
    // Attitude Representation
    sMRP_t ATT; // in MRPs
    sVEC_t OMG; // in rad/s
}sSTATE1_t;
sSTATE1_t Simulators_DiffEqn1(sSTATE1_t state, sVEC_t TRQ, sVEC_t MOI);
void Simulators_Run1();
void Simulators_Bdot();
void Simulators_Pointing();
#endif
