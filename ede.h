#ifndef EDE_H
#define EDE_H


class EdE
{
public:
    EdE();
    double * E;
    double * P;
    double * Pl;
    double * Rho;
    double * Nb;
    static double from_MeV4_to_km2;
    static double from_km2_to_MeV4;
    static double from_MeV4_to_MeVfm3;
    static double from_MeVfm3_to_MeV4;
    static double from_MeVfm3_to_km2;
    int N;
    EdE(const double *,const  double *,const  double *,const  double *,const  double *, int);
    void convert_from_MeV4_to_km2();
    void convert_from_MeVfm3_to_MeV4();
    void convert_from_MeVfm3_to_km2();
    EdE(const EdE &);
    ~EdE();
};

#endif // EDE_H
