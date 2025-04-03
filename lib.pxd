cdef extern from "tov.h" namespace "tov":
    void tov_main "tov::main_process" (char [], char [], double, double , double , double, int);

cdef extern from "gamma.h" namespace "Gamma":
    void gamma_main "Gamma::main_process" (char namein[],char nameout[],double rhomin,
            double hmin,double hmax,double rmax, int units);

cdef extern from "cilindric.h" namespace "cilindric":
    void cilindric_main "cilindric::main_process" (char namein[],char nameout[],double rhomin,
            double hmin,double hmax,double rmax, int units);

cdef extern from "general_cilindric.h" namespace "general_cilindric":
    void general_cilindric_main "general_cilindric::main_process" (char namein[],char nameout[],double rhomin,
            double hmin,double hmax,double rmax, int units);
