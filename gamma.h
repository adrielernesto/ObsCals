#ifndef GAMMA_H
#define GAMMA_H

#include "constfunc.h"

namespace Gamma {

    struct GIntegrator{
        mutable double g = 1;
        const Spline3 E ;
        const Spline3 El ;
        const Spline3 Rho ;
        const Spline3 Rhol ;
        GIntegrator(const EdE &);
        void operator() (const double, const double [] , double []) const ;
     };

    void main_process(char namein[],char nameout[],double rhomin,
                     double hmin,double hmax,double rmax, int units);

};





#endif // GAMMA_H
