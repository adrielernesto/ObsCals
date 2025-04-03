#ifndef GENERAL_CILINDRIC_H
#define GENERAL_CILINDRIC_H

#include "spline.h"
#include "math.h"
#include "time.h"
#include "file_processor.h"
#include "constfunc.h"

namespace general_cilindric {

    const double MS  = 1.477;  //masa del sol en km

    class GCIntegrator
    {
        public :
            const Spline3 E;
            const Spline3 Pl;
            const Spline3 Rho;
            const Spline3 Nb;
            GCIntegrator(const EdE &);
            void operator() (double, const double *, double *) const;
    };

    void main_process(char [], char [], double, double , double , double, int);
};

#endif