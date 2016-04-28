/*======================================================================================
Rich Brower Sat May 28 00:23:17 EDT 2011

C program based on simple FMV code of  S. F. McCormick, Multigrid Methods, Frontiers in Applied! Mathematics,
vol. 3, SIAM Books, Philadelphia, 1987.The code is intentionally put into a single file with no input parameters
to keep is concise. Of course a full code may have more traditional C structures!

We solve the 2-d Laplace problem:    A phi = b

 (A phi)((x,y) = [4 phi(x,y) - phi(x+1,y) - phi(x-1,y) - phi(x,y+1)  -phi(x,y+1)]/a^2 + m^2 phi(x,y) = b(x,y)

 Multiply by scale = 1/(4 + a^2 m^2)

     phi  =  (1 - scale*A) phi +  scale*b = phi + res
          =     scale * (phi(x+1,y) + phi(x-1,y) + phi(x,y+1) + phi(x,y+1))  +  a^2 scale*b(x,y)

where we use rescaled: res =  a^2 scale b - scale* A phi

with scale(lev) = 1/(4 + pow(2,lev)^2 m^2)      or  the lattice spacing is   a(lev) = pow(2,lev)

So the code is really sovling the new matrix iteration:

      phi(x,y) = [(1- scale(lev) * A)phi](x,y) + res(x,y)

               =  scale(lev) * (phi(x+1,y) + phi(x-1,y) + phi(x,y+1) + phi(x,y+1)) + res

At present relex iteration uses Jacobi method.
======================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "helpers.h"
#include "serial_mg.h"
#include "parallel_mg.h"

double const_func(int x);
double sin_func(int x);
struct timespec diff(struct timespec start, struct timespec end);
int main(int argc, char **argv)
{
  struct timespec time1, time2, total;

  // PARALLEL
  printf("Running Parallelized Multigrid\n");
  clock_gettime(CLOCK_MONOTONIC, &time1);
  parallel_multigrid(7, sin_func);
  clock_gettime(CLOCK_MONOTONIC, &time2);
  total = diff(time1, time2);
  printf("  Total time: %ld.%ld seconds\n", total.tv_sec, total.tv_nsec);

  // NON-PARALLEL
  // printf("Running Serial Multigrid\n");
  // clock_gettime(CLOCK_MONOTONIC, &time1);
  // serial_multigrid(7);
  // clock_gettime(CLOCK_MONOTONIC, &time2);
  // total = diff(time1, time2);
  // printf("  Serial Total time: %ld.%ld seconds\n", total.tv_sec, total.tv_nsec);

  return 0;
}

double const_func(int x) {
  return 4.0;
}

double sin_func(int x) {
  return sin( ( (double)(2.0 * M_PI) /256 ) * x );
}

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}
