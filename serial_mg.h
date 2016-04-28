#ifndef SERIAL_MG_H
#define SERIAL_MG_H

#include "helpers.h"

int multigrid(int Lmax, double (*boundary_func)(int));

#endif
