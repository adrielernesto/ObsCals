#include "ede.h"

double EdE::from_MeV4_to_km2 = 1.717e-13 ;
double EdE::from_km2_to_MeV4 = 5.824e+12;
double EdE::from_MeVfm3_to_MeV4 = 7.68035e+6;
double EdE::from_MeV4_to_MeVfm3 = 1.303e-7;
double EdE::from_MeVfm3_to_km2 = 1.31872e-19;

EdE::EdE()
{
    E = nullptr;
    Pl = nullptr;
    P = nullptr;
    Rho = nullptr;
    Nb = nullptr;
    N = 0;
}

EdE::EdE(const double * E0,const  double * Pl0,const  double * P0,
         const  double * Rho0,const  double * Nb0, int n)
{
    int cres = 1;
   while (P0[cres] >= P0[cres + 1])
      cres++;
    N = n - cres + 1;
    P = new double[N+1];
    E = new double [N+1];
    Pl = new double [N+1];
    Rho = new double [N+1];
    Nb = new double [N+1];


    for (int i = 1; i <= N; i++)
    {
        P[i] = P0[i + cres - 1];
        Pl[i] =  Pl0[i + cres - 1];
        Rho[i] = Rho0[i + cres - 1];
        Nb[i] =  Nb0[i + cres - 1];
        E[i] =  E0[i + cres - 1];
    }

}


EdE::EdE(const EdE & ede)
{
    N = ede.N;
    E = new double [N + 1];
    P = new double [N+1];
    Pl = new double [N+1];
    Rho = new double [N+1];
    Nb = new double [N+1];

    for (int i = 1; i <= N; i++)
    {
        E[i] = ede.E[i];
        Pl[i] = ede.Pl[i];
        P[i] = ede.P[i];
        Rho[i] = ede.Rho[i];
        Nb[i] = ede.Nb[i];
    }
}

void EdE::convert_from_MeV4_to_km2()
{
    for (int i = 1; i <= N; i++)
    {
        P[i] *= from_MeV4_to_km2;
        E[i] *= from_MeV4_to_km2;
        Pl[i] *= from_MeV4_to_km2;
        Nb[i] *= from_MeV4_to_MeVfm3;
    }
}



void EdE::convert_from_MeVfm3_to_km2()
{
    for (int i = 1; i <= N; i++)
    {
        P[i] *= from_MeVfm3_to_km2;
        E[i] *= from_MeVfm3_to_km2;
        Pl[i] *= from_MeVfm3_to_km2;
    }
}

void EdE::convert_from_MeVfm3_to_MeV4()
{
    for (int i = 1; i <= N; i++)
    {
        P[i] *= from_MeVfm3_to_MeV4;
        E[i] *= from_MeVfm3_to_MeV4;
        Pl[i] *= from_MeVfm3_to_MeV4;
        Nb[i] *= from_MeVfm3_to_MeV4;
    }
}

EdE::~EdE()
{
    delete [] E;
    delete [] Pl;
    delete [] P;
    delete [] Rho;
    delete [] Nb;
}


