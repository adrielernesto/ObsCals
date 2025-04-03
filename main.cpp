#include "cilindric.h"
#include "tov.h"
#include "gamma.h"
#include "general_cilindric.h"

void get_dat(char namein1[],char nameout1[],double &rhomin,
             double & hmin, double & hmax, double & rmax, int & units, int & metrics, double & g);


int main()
{
    char namein[1024], nameout[1024];
    int metrics, units;
    double rhomin, rmax, hmin, hmax,g = 1;

    get_dat(namein,nameout,rhomin,hmin,hmax,rmax,units,metrics,g);

    try
    {
        switch (metrics)
        {
            case 1 :
            {
                Gamma::main_process(namein,nameout,rhomin,hmin,hmax,rmax,units);
                break;
            }
            case 2 :
            {
                cilindric::main_process(namein,nameout,rhomin,hmin,hmax,rmax,units);
                break;
            }
            case 3 :
            {
                tov::main_process(namein,nameout,rhomin,hmin,hmax,rmax,units);
                break;
            }
            case 4:
            {
                general_cilindric::main_process(namein,nameout,rhomin,hmin,hmax,rmax,units);
            }
            
            default: throw "Invalid metric selection";
        }
    } catch (char const * e) {
        printf("%s\n",e);
    }


    return 0;
}

void get_dat(char namein[],char nameout[],double &rhomin,
             double & hmin, double & hmax, double & rmax, int & units, int & metrics,double & g)
{
    printf("Type input file:");
    scanf("%s",namein);
    
    printf("**The input file must be organized in 5 columns, for energy density, parallel pressure, perpendicular"
           "pressure, rest energy density (in g/cm^3) and baryon concentration.**\n");

    printf("Type output file:");
    scanf("%s",nameout);

    printf("Type minimun density:");
    scanf("%lf",&rhomin);

    printf("Type maximun radius expected:");
    scanf("%lf",&rmax);

    printf("Type minimun stepsize in integration:");
    scanf("%lf",&hmin);

    printf("Type maximun stepsize in integration:");
    scanf("%lf",&hmax);

    printf("Type 1 if you use MeV^4 y MeV^3 for energy density and "
          "concentration of baryons and 2 if you use MeV/fm^3 and 1/fm^3 instead:");
    scanf("%d",&units);

    printf("Type:\n1 to use gamma structure equations\n"
            "2for cilindric (Daryel)\n"
            "3 for TOV\n"
            "4 for general cilindric:");
    scanf("%d",&metrics);


}

















