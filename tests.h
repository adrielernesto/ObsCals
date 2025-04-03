#ifndef TESTS_H
#define TESTS_H

#include "tov.h"
#include "gamma.h"
#include "cilindric.h"
#include "gammar.h"

void test_sqm()
{
    char  namesin[2] [1024] = {
        "Inputs/sqm_B0_Bag65.dat",
        "Inputs/sqm_B5e17_Bag75.dat"
    };

    char namesout[4] [1024] = {
            "Outputs/gobs_sqm_B0_Bag65.dat",
            "Outputs/gobs_sqm_B5e17_Bag75.dat"
    };

    double rhomin = 0, rmax = 30, hmax = 1, hmin = 5e-3, units = 1;

    try {
        for(int i = 0; i < 2; i++)
            Gamma::main_process(namesin[i],namesout[i],rhomin,hmin,hmax,rmax,units);
    } catch (char const * e) {
        printf("%s\n",e);
    }
    system("cd Outputs && python3 comparaciones_sqm.py");
}




void test_gamma()
{
    char namein1[80] = "Inputs/bosons_T0_B0_a1.dat";
    char nameout1[80] = "Outputs/gobs_bosons_T0_B0_a1.dat";

    char namein5[80] = "Inputs/bosons_T0_B517_a1.dat";
    char nameout5[80] = "Outputs/gobs_bosons_T0_B517_a1.dat";

    double rhomin = 5e+13, hmax = 1, hmin = 1e-3, rmax = 16, units = 1;

    try {
        Gamma::main_process(namein1,nameout1,rhomin,hmin,hmax,rmax,units);
        Gamma::main_process(namein5,nameout5,rhomin,hmin,hmax,rmax,units);
    //    Gamma::main_process(namein2,nameout2,0,5,10,1000,units);
    } catch (char const * e) {
        printf("%s\n",e);
    }
    system("cd Outputs && python3 comparaciones_gamma.py");
}


void test_gammar()
{
    char namein1[80] = "Inputs/bosons_T0_B0_a1.dat";
    char nameout1[80] = "Outputs/gobs_bosons_T0_B0_a1.dat";

    char namein5[80] = "Inputs/bosons_T0_B517_a1.dat";
    char nameout5[80] = "Outputs/gobs_bosons_T0_B517_a1.dat";

    double rhomin = 5e+13, h = 0.001, rmax = 20, units = 1;

    try {
        Gammar::main_process(namein1,nameout1,rhomin,h,rmax,units);
        Gammar::main_process(namein5,nameout5,rhomin,h,rmax,units);
    //    Gamma::main_process(namein2,nameout2,0,5,10,1000,units);
    } catch (char const * e) {
        printf("%s\n",e);
    }
    system("cd Outputs && python3 comparaciones_gamma.py");
}

void test_tov()
{
    char namein1[80] = "Inputs/bosons_T0_B0_a1.dat";
    char nameout1[80] = "Outputs/tobs_bosons_T0_B0_a1.dat";

    char namein5[80] = "Inputs/bosons_T0_B0_a5.dat";
    char nameout5[80] = "Outputs/tobs_bosons_T0_B0_a5.dat";

  //  char namein2[80] = "Inputs/eos_tesis_B0.dat";
  //  char nameout2[80] = "Outputs/tobs_eos_tesis_B0.dat";

    double rhomin = 5e+13, hmax = 1, hmin = 1e-3, rmax = 50, units = 1;

    try {
        tov::main_process(namein1,nameout1,rhomin,hmin,hmax,rmax,units);
        tov::main_process(namein5,nameout5,rhomin,hmin,hmax,rmax,units);
    //    tov::main_process(namein2,nameout2,0,5,10,1000,units);
    } catch (char const * e) {
        printf("%s\n",e);
    }

    system("cd Outputs && python3 comparaciones_tov.py");
}


void test_cilindric()
{
    char namein[80] = "Inputs/eos_tesis_B5d12.dat";
    char nameout[80] = "Outputs/cobs_eos_tesis_B5d12.dat";

    char namein2[80] = "Inputs/eos_tesis_B0.dat";
    char nameout2[80] = "Outputs/cobs_eos_tesis_B0.dat";
    double rhomin = 0, hmin = 0.11, hmax = 5, rmax = 10000, units = 1;

    try {
        cilindric::main_process(namein,nameout,rhomin,hmin,hmax,rmax,units);
        cilindric::main_process(namein2,nameout2,rhomin,hmin,hmax,rmax,units);
    } catch (char const * e) {
        printf("%s\n",e);
    }

    system("cd Outputs && python3 comparaciones_cil.py");
}

#endif // TESTS_H
