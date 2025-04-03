#ifndef ODEINT_H
#define ODEINT_H

#include "math.h"


inline double max(double a, double b) { return a > b ? a : b;}
inline double min(double a, double b) { return a < b ? a : b;}

template<class T>//Método de Runge-Kutta de orden 5 con los coeficientes de Kash Carp que que se va a usar dentro de rkqs.
void rkck_embedded(const double y[], double yf[],double x, const T &f,double h,double yerr[],int N){
    /*Recibe un array y con los valores iniciales de las variables dependientes, el valor inicial x de la variable
    independiente, el objeto f que evalúa el miembro derecho del sistema de edos, N es el número de variables
    dependientes, al final yf es reemplazado por los valores de las variables dependientes en el punto x+h, y yerr
    por las estimaciones del error de truncamiento en el proceso.
    */
    static float a2=1.0/5.0 , a3=0.3,a4=0.6,a5=1,a6=7.0/8.0 ,c1=37.0/378.0,c2=0,c3=250.0/621,
        c4=125.0/594.0,c5=0,c6=512.0/1771.0,dc1=c1-2825.0/27648,dc2=0,dc3=c3-18575.0/48384,
        dc4=c4-13525.0/55296,dc5=c5-277.0/14336,dc6=c6-0.25,b21=0.2,b31=3.0/40,b32=9.0/40,
        b41=0.3,b42=-0.9,b43=1.2,b51=-11.0/54,b52=2.5,b53=-70.0/27 ,b54=35.0/27,b61=1631.0/55296,b62=175.0/512,
        b63=575.0/13824,b64=44275.0/110592,b65=253.0/4096;

        double ytemp[N],k1[N],k2[N],k3[N],k4[N],k5[N],k6[N];

        f(x,y,k1);
        for(int i= 0;i<N;i++)
            ytemp[i] = y[i] + h*k1[i]*b21;
        f(x+a2*h,ytemp,k2);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b31*h*k1[i]+b32*h*k2[i];
        f(x+a3*h,ytemp,k3);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b41*h*k1[i]+b42*h*k2[i]+b43*h*k3[i];
        f(x+a4*h,ytemp,k4);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b51*h*k1[i]+b52*h*k2[i]+b53*h*k3[i]+b54*h*k4[i];
        f(x+a5*h,ytemp,k5);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b61*h*k1[i]+b62*h*k2[i]+b63*h*k3[i]+b64*h*k4[i]+b65*h*k5[i];
        f(x+a6*h,ytemp,k6);

        for(int i=0;i<N;i++)
            yf[i]=y[i]+h*(c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]);

        for(int i=0;i<N;i++)
            yerr[i]=dc1*k1[i]+dc2*k2[i]+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i];
}




template<class T>
void rkqs(double *y, double &x,const T &f, double &htry,int N,double *tol,double hmin,double hmax){
    /*Este método utiliza el Runge-Kutta de orden 5 para calcular integrar y estimar el error de truncamiento, va cambiando
    h para que el error se mantenga menor que tol.
    Argumentos :
    y : array con los valores iniciales de las variables independientes.
    x : referencia al valor inicial de la variable independiente.
    f : referencia al objeto que evalúa el miembro derecho del sistema de edos, debe tener un operator().
    N : cantidad de variables dependientes.
    tol : array con la cota superior deseada para el error de tuncamiento en cada variable dependiente.
    htry : referencia al valor que se propone utilizar par el stepsize y que será modificado en consecuencia.
    hmin - hmax: valores min y max que puede tener el stepsize, respectivamente.
    */
    double errmax , htemp , ytemp[N] , yerr[N] , h = htry;

    while(true){
        rkck_embedded<T>(y,ytemp,x,f,h,yerr,N);//se utiliza el rkck para obtener y en x+h y el error se guarda en yerr
        errmax = 0.0;
        for(int i=0;i<N;i++) errmax=max(errmax,abs(yerr[i]/tol[i]));//se comprueba se algún error excede la tolerancia deseada.
        errmax += 1e-12;//se le suma esto para garantizar que no sea nulo.
        if(errmax <= 1.0) break;//si todo está en regla, se para el ciclo.
        htemp = 0.95*h*pow(errmax,-0.25);//Si no se disminuye h sin que se haga más pequeño que el valor mínimo permotido
        h= max(htemp,hmin);//
        if (h == hmin) break;//En caso de que el valor mínimo no garantice que errmax <= 1, igual se detiene el ciclo y se continúa

    }
    htemp = min(0.95*h*pow(errmax,-0.2),hmax);//htemp  no puede superar a hmax.
    //se estima cuánto puede aumentar h, en caso de que el ciclo se halla detenido en la línea 71, htemp>=h; si se detuvo
    //en la línea 74,htemp <= hmin.
    x += h;//se aumenta x según el valor de h que se utilizó.
    htry = min(htemp,hmin);//se modifica htry para la próxima evaluación de rkqs.
    for(int i = 0; i < N;i++)
        y[i] = ytemp[i];
}






















template<class T>//Método de Runge-Kutta de orden 5 con los coeficientes de Kash Carp.
void rkck(double y[],double  & x, const T &f,double h,int N){
    /*Recibe un array y con los valores iniciales de las variables dependientes, el valor inicial x de la variable
    independiente, el objeto f que evalúa el miembro derecho del sistema de edos, N es el número de variables
    dependientes, al final y es reemplazado por los valores de las variables dependientes en el punto x+h,
    y x por x + h.
    */
    static float a2=1.0/5.0 , a3=0.3,a4=0.6,a5=1,a6=7.0/8.0 ,c1=37.0/378.0,c2=0,c3=250.0/621,
        c4=125.0/594.0,c5=0,c6=512.0/1771.0,dc1=c1-2825.0/27648,dc2=0,dc3=c3-18575.0/48384,
        dc4=c4-13525.0/55296,dc5=c5-277.0/14336,dc6=c6-0.25,b21=0.2,b31=3.0/40,b32=9.0/40,
        b41=0.3,b42=-0.9,b43=1.2,b51=-11.0/54,b52=2.5,b53=-70.0/27 ,b54=35.0/27,b61=1631.0/55296,b62=175.0/512,
        b63=575.0/13824,b64=44275.0/110592,b65=253.0/4096;

        double ytemp[N],k1[N],k2[N],k3[N],k4[N],k5[N],k6[N];

        f(x,y,k1);
        for(int i= 0;i<N;i++)
            ytemp[i] = y[i] + h*k1[i]*b21;
        f(x+a2*h,ytemp,k2);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b31*h*k1[i]+b32*h*k2[i];
        f(x+a3*h,ytemp,k3);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b41*h*k1[i]+b42*h*k2[i]+b43*h*k3[i];
        f(x+a4*h,ytemp,k4);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b51*h*k1[i]+b52*h*k2[i]+b53*h*k3[i]+b54*h*k4[i];
        f(x+a5*h,ytemp,k5);

        for(int i=0;i<N;i++)
            ytemp[i]=y[i]+b61*h*k1[i]+b62*h*k2[i]+b63*h*k3[i]+b64*h*k4[i]+b65*h*k5[i];
        f(x+a6*h,ytemp,k6);

        for(int i=0;i<N;i++)
            y[i]=y[i]+h*(c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]);
        x = x + h;
}













#endif // ODEINT_H
