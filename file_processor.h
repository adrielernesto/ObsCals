#ifndef FILE_PROCESSOR_H
#define FILE_PROCESSOR_H


#include <cstdio>
#include "ede.h"

int number_of_lines(FILE * file);

void get_columns(FILE*,double[], double [], double[], double[], double[],int);

EdE get_ede(const char *);

#endif // FILE_PROCESSOR_H
