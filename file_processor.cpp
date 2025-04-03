#include "file_processor.h"


int number_of_lines(FILE * file)
{
    rewind(file);
    int n = 0;
    while ( ! feof(file) )
        if ( 10 == fgetc(file) )
            n++;
    rewind(file);
    return n;
}


void get_columns(FILE* file,double E[] , double Pl[], double P[], double Rho[], double Nb[],int N)
{
    rewind(file);
    for (int i = 1; i <= N; i++)
        fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\n",
               E+i, Pl+i, P+i, Rho+i,Nb+i);
    rewind(file);
}



EdE get_ede(const char * name)
{
    FILE * f = fopen(name,"r");
    if ( ! f)
        throw "Input file not opend";
    int n = number_of_lines(f);
    double E[n+1], Pl[n+1], P[n+1], Rho[n+1], Nb[n+1];
    get_columns(f,E,Pl,P,Rho,Nb,n);
    fclose(f);
    return EdE(E,Pl,P,Rho,Nb,n);
}



























