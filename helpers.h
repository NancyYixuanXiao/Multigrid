#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
  int N;
  int Lmax;
  int size[20];
  double a[20];
  double alpha;;
} param_t;

int getIndex(int x, int y, int width);

#endif
