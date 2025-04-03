#ifndef SPLIN3_H
#define SPLIN3_H


#include "stdio.h"


/*Los objetos Spline3 utilizan el método de spline cúbico para interpolar.
Reciben en el contructor dos arrays x1,x2,..,xN ; y1,y2,..,yN, con los puntos exactos;
y cuenta con el operrador () sobrecargado que recibe un double xx y devuelve el resultado de la
interpolación en ese punto.
Las funciones spline y splint fueron tomadas del libro Numerical Recipies in C y yo lo convertí
en un objeto para una mejor manipulación :) */

class  Spline3{
    private:
        int N;
        double *x ,*y ,*y2;
        void spline(const double *,const double* , const int , const double ,const double , double *) const;
        double splint(const double *,const double *, const double *, const int ,const double ) const;
    public:
        double xmin , xmax ,ymin , ymax;
        Spline3();
        Spline3(const double *, const double *, const int);
        Spline3(const Spline3 &);
        double operator() (const double xx) const;
        ~Spline3();
};





#endif // SPLIN3_H
