#include "tov.h"

namespace tov {

    TIntegrator::TIntegrator(const EdE & EoS) :
        E(EoS.P, EoS.E, EoS.N) , Rho(EoS.P, EoS.Rho, EoS.N)
    {
        //
    }

    void TIntegrator::operator()(double r, const double * y, double * yf) const
    {
        //y = [p,m,mb]
        double p = y[0] , m = y[1];
        double e = E(p);
        if (r < 1e-5 )
            yf[0]  =  yf[2]=  yf[1] = 0;
        else
        {
            yf[0] = -R0*(p+e)*(m + 4*pi*p*pow(r*5.07e+15,3)/MS )/(r*r-2*m*R0*r);
            yf[2] = 4*pi*r*r*Rho(p)/sqrt(1-2*R0*m/r) *5.61e+41/MS;
        }
         yf[1] = 4*pi*r*r*e*pow(5.07e+15,3)/MS;
    }


    void main_process(char namein[],char nameout[],double rhomin,
                     double hmin,double hmax,double rmax, int units)
    {
        double start_time = clock();
        EdE EoS = get_ede(namein);

        switch (units)
        {
            case 1 :
            {
                //Ya están en MeV⁴
                break;
            }
            case 2 :
            {
                //Se convierten de MeV/fm³ a MeV⁴
                EoS.convert_from_MeVfm3_to_MeV4();
                break;
            }
            default: throw "Error in units";
        }

        int init;
        for (init = 1; init <= EoS.N; init++)
            if ( EoS.P[init] > 0 )
                break;
        /* Se encuentra el primer valor de presión positiva */

        double pmin = max(0.0, EoS.P[1]);

        TIntegrator ft = TIntegrator(EoS);

        FILE * fout = fopen(nameout,"w");
        fprintf(fout, "#%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
                "E[MeV/fm^3]","Rho[g/cm^3]","M/MS","R[km]","Mb/MS","I[MSkm^2]","z"
                );

        double r, h, M, Mb, p0, e0, rho0, z, I, eps = 1e-3, tol[3], dy[3];

        for (int i = init + 1; i <= EoS.N; i++)
        {
            if ( (rho0 = EoS.Rho[i]) < rhomin)
                continue;
            printf("Star # %10d\n",i);
            r = hmin;
            e0 = EoS.E[i];
            Mb = 4.0/3.0*pi*r*r*r*rho0*1e+15/MSg ;
            M = 4.0/3.0*pi*r*r*r*e0*EdE::from_MeV4_to_MeVfm3*1e+54 /MS;
            h = hmin;
            p0 = EoS.P[i];


            double y[3] = {p0, M, Mb};

            while (true)
            {

                M = y[1]; Mb = y[2];
                for(int it = 0; it < 3; it++)
                    tol[it] = eps*y[it] + 1e-4;
                tol[0] += pmin;

              //  ft(r,y,dy);
              //  if (y[0] + hmin * dy[0] <= pmin)
              //      break;
                rkqs<TIntegrator>(y,r,ft,h,3,tol,hmin,hmax);

                if (y[0] <= pmin) break;

                if ( r > rmax)
                {
                    printf("Failed to converge with Ec = %10.3e [MeV/fm^3]\n",
                           e0*EdE::from_MeV4_to_MeVfm3);
                    goto next_start;
                }
            }

            if (false) next_start : continue;

            I = 1.0 / 5.0 * M * r * r;
            z = 1.0 / sqrt(1 - 2 * R0 * M / r) - 1;

            fprintf(fout, "%10.5e\t%10.5e\t%10.5lf\t%10.5lf\t%10.5lf\t%10.5e\t%10.5e\n",
                        e0*EdE::from_MeV4_to_MeVfm3,
                        rho0,
                        M,
                        r,
                        Mb,
                        I,
                        z
                    //"E[MeV/fm^3]","Rho[g/cm^3]","M/MS","R[km]","Mb/MS","I[MSkm^2]","z"
                                 );
        }
        fclose(fout);

        double time = (clock() - start_time) / CLOCKS_PER_SEC;

        printf("-----------------------------------------------------------------------\n"); //71
        printf("%-20s:%50s\n","Input file",namein);
        printf("%-20s:%50s\n","Output file",nameout);
        printf("%-20s:%50d\n","Used lines",EoS.N);
        printf("%-20s:%50d\n","Used time (sec)",(int) time);
        printf("-----------------------------------------------------------------------\n"); //71




    }


};
