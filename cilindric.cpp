#include "cilindric.h"

namespace cilindric {

    CIntegrator::CIntegrator(const EdE & ede) : E(ede.P,ede.E,ede.N),
        Pl(ede.P,ede.Pl,ede.N),  Rho(ede.P,ede.Rho,ede.N), Nb(ede.P,ede.Nb,ede.N)
    {
        //
    }

    void CIntegrator::operator()(double r, const double * y, double * yf) const
    {
            /* y = {P,f,a,t,df,dt,l,lb}  *
            P : presión perpendicular,
            f: coeficiente metrico relativo al tiempo,
            a: coeficiente metrico relativo al radio,
            t: coeficiente metrico relativo a z,
            df,dt : sus derivadas,
            l,lb: M/Z , Mb/Z
            r: coordenada radial.
            */
            if(r == 0)
            {
                yf[0] = yf[1] = yf[2] = yf[3] = yf[4] = yf[5] = yf[6] = yf[7] = 0;
                return ;
            }
            double p = y[0], f = y[1], a = y[2], t = y[3], df = y[4], dt = y[5];
            double pl = Pl(p) , e = E(p) , rho = Rho(p) , exp2a = exp(2*a);
            double dp = -p*(df+dt) - df*e + dt*pl;
            double da =  dt + df - 4*pi*r*exp2a*pl +  4*pi*r*exp2a*e  ;
            double ddf =  4*pi*exp2a*(e + pl + 2*p + r*df*e) - df*4*pi*r*pl*exp2a - df/r ;

            double ddt = -4*pi*exp2a*(e-r*dt*e+pl+r*dt*pl-2*p) - dt/r;
            double dl = 4*pi*r*exp(a+t+f)*(e-2*p-pl)/MS    ;
            double dlb = 4*pi*rho*r*exp(a+f+t)/MSg  * 1e+15;
            yf[0] = dp; yf[1] = df ; yf[2] = da; yf[3] = dt; yf[4] = ddf;
            yf[5] = ddt ; yf[6] = dl ; yf[7] = dlb;
            if (fabs(dp) > 1e200 || fabs(da) > 1e200 || fabs(ddf) > 1e200 || fabs(ddt) > 1e200 || fabs(dl) > 1e200 || fabs(dlb) > 1e200)
                printf("Overflow en derivadas en r = %lf.\n",r);
        }



    void main_process(char namein[],char nameout[],double rhomin,
                     double hmin,double hmax,double rmax, int units)
    {
        double start_time = clock();
        EdE EoS = get_ede(namein);

        switch (units) {
            case 1 :
            {
                EoS.convert_from_MeV4_to_km2();
                break;
            }
            case 2 :
            {
                EoS.convert_from_MeVfm3_to_km2();
                break;
            }
            default: throw "Error in units";
        }

        int init;
        for (init = 1; init <= EoS.N; init++)
            if ( EoS.P[init] > 0 && EoS.Pl[init] > 0)
                break;
        /* Se encuentra el primer par de presiones positivas */

        FILE * fout = fopen(nameout,"w");

        double pmin = max(0.0, EoS.P[1]);

        CIntegrator fc = CIntegrator(EoS);

        fprintf(fout,"#%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
                "E[MeV/fm^3]","Nb[fm^-3]","Rho[g/cm^3]","R[km]","M/Z[MS/km]","Mb/Z[MS/km]");

        double r, h, eps = 1e-6, lb, l, p0, e0, pl0, rho0, nb0, dy[8];
        /*
        *  r : distancia al centro.
        *  h : paso de la integracion, o sea, el dr.
        *  hmin <= h <= hmax
        *  ytol: es la tolerancia absoluta deseada; o sea, en cada paso h cambia de valor para
        *        mantener el error menor que ytol.
        *  eps: especie de error relativo.
        *  l , lb : M/Z , Mb/Z
        */

        for (int i = init + 1; i <= EoS.N; i++)
        {
            if (EoS.Rho[i] < rhomin)
                continue;
            printf("Star # %10d\n",i);
            h = hmax;
            r = 0.0;
            l = lb = 0.0;
            e0 = EoS.E[i];
            p0 = EoS.P[i];
            pl0 = EoS.Pl[i];
            nb0 = EoS.Nb[i];
            rho0 = EoS.Rho[i];
            //y = {P,f,a,t,df,dt,l,lb}
            double y[8] = {p0, 0, 0, 0, 0, 0, l, lb};
            double ytol[8] = {eps*p0 + pmin};
            for (int it = 1; it < 8; it++)
                ytol[it] = fabs(eps*y[it]) + eps;

            while (true)
            {
                l = y[6]; lb = y[7];
                rkqs<CIntegrator>(y,r,fc,h,8,ytol,hmin,hmax);
                //  fc(r,y,dy);
                if (y[0] <= pmin)//(y[0] + dy[0] * hmin <= pmin )
                    break;
                if ( r > rmax)
                {
                    printf("Failed to converge with E_c = %.2e MeV/fm^3\n",
                                  e0*EdE::from_MeV4_to_MeVfm3 *EdE::from_km2_to_MeV4);
                    goto next_star;
                }

            }

            if (false) next_star : continue;

            if (l <= 0) break;

            fprintf(fout, "%10.5e\t%10.5e\t%10.5e\t%10.5lf\t%10.5e\t%10.5e\n",
                        e0*EdE::from_MeV4_to_MeVfm3 *EdE::from_km2_to_MeV4,
                        nb0,
                        rho0,
                        r,
                        l,
                        lb
                    );
            //"E[MeV/fm^3]","Nb[fm^-3]","Rho[g/cm^3]","R[km]","M/Z[MS/km]","Mb/Z[MS/km]"
            printf("\r\r\r\r\r\r\r\r\r\r");
        }


       fclose(fout);

       double total_time = ( clock() - start_time) / CLOCKS_PER_SEC;

       printf("----------------------------------------------------------------------------\n"); //71
       printf("%-20s:%55s\n","Input file",namein);
       printf("%-20s:%55d\n","Used lines", EoS.N);
       printf("%-20s:%55d\n","Used time (sec)",(int) total_time);
       printf("----------------------------------------------------------------------------\n"); //71
    }
};
