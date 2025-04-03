#include "spline.h"


Spline3::Spline3()
{
    N = xmin = xmax = ymin = ymax = 0;
    x = new double[N+1];
    y = new double[N+1];
    y2 = new double[N+1];
}


void Spline3::spline(const double x[],const double y[] ,const int n , const double yp1,const double ypn, double y2[]) const
{
    int i,k;
    double p,qn,sig,un, u[n+1];
    if(yp1 > 1e+30) y2[1] = u[1] = 0.0;
    else {
        y2[1] = -0.5;
        u[1] = (3.0/(x[2] - x[1]))*((y[2] - y[1])/(x[2] - x[1]) - yp1);
    }
    for(i = 2;i <= n; i++){
        sig=(x[i] - x[i-1])/(x[i+1]-x[i-1]);
        p = sig*y2[i-1] + 2.0;
        y2[i] = (sig-1.0)/p;
        u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i+1])/(x[i]-x[i+1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
    }
    if(ypn > 1e+30) qn=un=0.0;
    else {
        qn = 0.5;
        un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n] - x[n-1]));
    }
    y2[n] = (un - qn*u[n-1])/(qn*y2[n-1]+1.0);
    for(k=n-1;k>=1;k--)
        y2[k] = y2[k]*y2[k+1]+u[k];

}


double Spline3::splint(const double xa[],const double ya[], const double y2a[],const int n,const double x) const
{
   /* if(x < xmin || x > xmax){
        printf("Valor fuera de rango en Spline3::splint()");
        if(x < xmin) printf(" (too small)\n");
        else printf(" (too big).\n");
    }
    */

    if (x < xmin) return -0;
    if (x > xmax) return ya[N];

    int klo,khi,k;
    double h,b,a;
    klo = 1;
    khi = n;
    while(khi -klo > 1) {
        k = (khi+klo) /2;
        if(xa[k] > x) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    double y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(h*h)/6.0;
    return y;
}


Spline3::Spline3(const double * xa, const double * ya, const int n){
    N = n;
    x = new double[N+1];
    y = new double[N+1];
    y2 = new double[N+1];
    for(int i = 0; i <= N ; i++)
    {
        x[i] = xa[i];
        y[i] = ya[i];
    }
    spline(xa,ya,N,0,0,y2);
    xmin = x[1];
    xmax = x[N];
    ymin = y[1];
    ymax = y[N];
}

double Spline3::operator() (const double xx) const
{
    return splint(x,y,y2,N,xx);
}

Spline3::Spline3(const Spline3 &ori)
{
    N = ori.N;
    delete [] x;
    delete [] y;
    delete [] y2;
    x = new double[N+1];
    y = new double[N+1];
    y2 = new double [N+1];
    for (int i = 1; i <= N ; i++)
    {
        x[i] = ori.x[i];
        y[i] = ori.y[i];
        y2[i] = ori.y2[i];
    }
    xmin = ori.xmin; xmax = ori.xmax;
    ymin = ori.ymin; ymax = ori.ymax;
}

Spline3::~Spline3()
{
    delete [] x;
    delete [] y;
    delete [] y2;
}
