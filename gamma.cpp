#include "gamma.h"

namespace Gamma {

    GIntegrator::GIntegrator(const EdE & EoS) :
        E(EoS.P, EoS.E, EoS.N), El(EoS.Pl, EoS.E, EoS.N),
        Rho(EoS.P,EoS.Rho,EoS.N), Rhol(EoS.Pl,EoS.Rho,EoS.N)
    {
        g = 1;
    }

    void GIntegrator::operator() (double r, const double y[], double dydr[]) const
    {
        // y = {p, pl, m, mb}
        double p = y[0], pl = y[1], m = y[2];
        double e = E(p), el = El(pl);
        if (r == 0)
        {
            dydr[0] = dydr[1] = dydr[2] = dydr[3] = 0;
            return;
        }
        double power = pow(1 - 2*R0*m/r,g);

        dydr[0] = - (p + e)*(r/2 + 4*pi*r*r*r*p*G*5.07*5.07e+30 - r/2*power)/(r*r*power);
        dydr[1] = - (pl + el)*(r/2 + 4*pi*r*r*r*pl*G*5.07*5.07e+30 - r/2*power)/(r*r*power);
        dydr[2] = 2*pi/MS*r*r*(e+el)*pow(5.07e+15,3)*g;
        dydr[3] = 2*pi/MS*r*r*(Rho(p)+Rhol(pl))*5.609e+41/sqrt(power);
      //  printf("p,e,pl,el = %lf %lf %lf %lf\n",p,e,pl,el);

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
        int init = 1;
        for (init = 1; init <= EoS.N; init++)
            if (EoS.P[init] > 0 && EoS.Pl[init] > 0)
                break;
        double pmin = max(EoS.P[1],0.0);
        GIntegrator fg(EoS);
        FILE * fout = fopen(nameout,"w");
        fprintf(fout, "#%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
               "E[MeV/fm^3]","Rho[g/cm^3]","R[km]","Z[km]","M/MS","Mb/MS"
               );
        double r, M, Mb, h, eps = 1e-3, tol[4], p0, pl0, e0, rho0;
        for(int i = init + 1; i <= EoS.N; i ++)
        {
            if ((rho0 = EoS.Rho[i]) < rhomin )
                continue;
            printf("Star # %10d\n",i);
            r = hmin;
            p0 = EoS.P[i];
            pl0 = EoS.Pl[i];
            e0 = EoS.E[i]* EdE::from_MeV4_to_MeVfm3;
            fg.g = (p0 > 0 ? pl0 / p0 : 1);
            M = 4.0 / 3.0 * pi * r * r * r * e0 / MS * 1e+54 ;
            Mb = 4.0/3.0*pi*r*r*r*rho0*1e+15/MSg ;
            double y[4] = {p0, pl0, M, Mb};
            h = hmin;
            while (true)
            { // fg(r,y,tol);
              //  printf("%lf\n",tol[0]);
                M = y[2];
                Mb = y[3];
                tol[0] = pmin + eps*y[0];
                tol[1] = pmin + eps*y[1];
                tol[2] = eps + eps*y[2];
                tol[3] = eps + eps*y[3];
                rkqs<GIntegrator>(y,r,fg,h,4,tol,hmin,hmax);
                if (r > rmax)
                {
                    printf("Failed to converge with Ec = %10.2e MeV/fm^3\n",e0);
                    goto next_star;
                }
                if (y[0] <= pmin || y[1] <= 0)
                    break;
            }
            if (false) next_star : continue;
            fprintf(fout,"%10.5e\t%10.5e\t%10.5lf\t%10.5lf\t%10.5lf\t%10.5lf\n",
                  //  "E[MeV/fm^3]","Rho[g/cm^3]","R[km]","Z[km]","M/MS","Mb/MS"
                    e0,rho0,r,r*fg.g,M,Mb
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
