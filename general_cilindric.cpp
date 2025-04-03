#include "general_cilindric.h"

namespace general_cilindric {

    GCIntegrator::GCIntegrator(const EdE & ede) : E(ede.P,ede.E,ede.N),
        Pl(ede.P,ede.Pl,ede.N),  Rho(ede.P,ede.Rho,ede.N), Nb(ede.P,ede.Nb,ede.N)
    {
        //
    }

    void GCIntegrator::operator()(double r, const double * y, double * yf) const
    {
            /* y = {P,lambda,nu,dlambda,dnu,l,lb}  
            */
            if(r == 0)
            {
                yf[0] = yf[1] = yf[2] = yf[3] = yf[4] = yf[5] = yf[6] = 0;
                return ;
            }
            double p = y[0], lmb = y[1], nu = y[2], dlmb = y[3], dnu = y[4];
            double pl = Pl(p) , e = E(p) , rho = Rho(p) , exp2lmbmnu = exp(2*nu-2*lmb);

            double dp = - (pl+ e)*dlmb  - (-pl+ p)*dnu ;
            double ddnu = 8*exp(-2*lmb +2*nu) *pi * p - dlmb*dlmb;
            double ddlmb = 4*exp(-2*lmb +2*nu)  *pi * e - dlmb / r + 0.5*(dlmb*dlmb + ddnu*ddnu);

            double dl = 4*pi*r*exp(2*nu - 2*lmb)*(e-2*p-pl)/MS    ;
            double dlb = 4*pi*rho*r*exp(2*nu - 2*lmb)/MSg  * 1e+15;

            yf[0] = dp; yf[1] = dlmb ; yf[2] = dnu ; yf[3] = ddlmb; yf[4] = ddnu;
            yf[5] = dl ; yf[6] = dlb;

            for(int i = 0; i < 7; i++)
                if (abs(yf[i]) > 1e200)
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

        GCIntegrator fc = GCIntegrator(EoS);

        //fprintf(fout,"#%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
        //        "E[MeV/fm^3]","Nb[fm^-3]","Rho[g/cm^3]","R[km]","M/Z[MS/km]","Mb/Z[MS/km]","M[Ms]","Mb[MS]");

        double r, h, eps = 1e-6, lb, l, p0, e0, pl0, rho0, nb0, dy[7];
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
            double y[7] = {p0, 0, 0, 0, 0,  l, lb};
            double ytol[7] = {eps*p0 + pmin};
            for (int it = 1; it < 7; it++)
                ytol[it] = fabs(eps*y[it]) + eps;

            while (true)
            {
                l = y[5]; lb = y[6];
                rkqs<GCIntegrator>(y,r,fc,h,7,ytol,hmin,hmax);
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

            fprintf(fout, "%10.5e\t%10.5e\t%10.5e\t%10.5lf\t%10.5e\t%10.5e\t%lf\t%lf\n",
                        e0*EdE::from_km2_to_MeV4,
                        nb0,
                        rho0,
                        r,
                        l,
                        lb,
                        l*r,
                        lb*r
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

}