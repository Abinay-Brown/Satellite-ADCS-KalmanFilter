cmake_minimum_required(VERSION 3.19)
project(ADCS_MAIN C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(ADCS_MAIN
        main.c

        Physical_Constants.h

        StaticLinalg.c
        StaticLinalg.h

        DynamicLinalg.c
        DynamicLinalg.h


        DateTime.c
        DateTime.h

        IGRF.c
        IGRF.h
        IGRFcoeff.c
        IGRFcoeff.h
        Legendre.c
        Legendre.h

        SunVector.c
        SunVector.h

        Frames.c
        Frames.h

        Estimators.c
        Estimators.h

        Controllers.c
        Controllers.h

        Kalman.c
        Kalman.h

        Simulators.c
        Simulators.h

        Statistics.c
        Statistics.h
        )
