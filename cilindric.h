#ifndef CILINDRIC_H
#define CILINDRIC_H

#include "spline.h"
#include "math.h"
#include "time.h"
#include "file_processor.h"
#include "constfunc.h"

namespace cilindric {

    const double MS  = 1.477;  //masa del sol en km

    class CIntegrator
    {
        public :
            const Spline3 E;
            const Spline3 Pl;
            const Spline3 Rho;
            const Spline3 Nb;
            CIntegrator(const EdE &);
            void operator() (double, const double *, double *) const;
    };

    void main_process(char [], char [], double, double , double , double, int);
};
















#endif // CILINDRIC_H
