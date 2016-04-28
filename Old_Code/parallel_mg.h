#ifndef PARALLEL_MG_H
#define PARALLEL_MG_H

#include "helpers.h"

int parallel_multigrid(int Lmax, double (*boundary_func)(int));

#endif
