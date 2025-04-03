#ifndef TOV_H
#define TOV_H

#include "constfunc.h"

namespace tov {

    struct TIntegrator
    {
        const Spline3 E, Rho;
        TIntegrator(const EdE &);
        void operator() (double, const double[], double[]) const;
    };

    void main_process(char [], char [], double, double , double , double, int);

};


#endif // TOV_H
