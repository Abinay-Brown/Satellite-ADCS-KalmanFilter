

#ifndef KALMAN_H_KALMAN_H
#define KALMAN_H_KALMAN_H

void Kalman_Init();
void Kalman_Filter();


void inv6x6(); // This is just to test if the inverse is working.
// Not used in the Kalman routine directly.


void magnetometer_measure();
void sunSensor_measure();
void rateGyro_measure();
#endif