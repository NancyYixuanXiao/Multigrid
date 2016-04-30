#include <omp.h>
#include "parallel_mg.h"
#include "helpers.h"

void relax_parallel(double *phi, double *res, int lev, int niter, param_t p);
void proj_res_parallel(double *res_c, double *rec_f, double *phi_f, int lev,param_t p);
void inter_add_parallel(double *phi_f, double *phi_c, int lev, param_t p);
double GetResRoot_parallel(double *phi, double *res, int lev, param_t p);

int parallel_multigrid(int Lmax, double (*boundary_func)(int)) {
  double *phi[20], *res[20];
  param_t p;
  int nlev;
  int i,j,lev;

  //set parameters________________________________________
  p.Lmax = Lmax; // max number of levels
  p.N = 2*(int)pow(2,p.Lmax);  // MUST BE POWER OF 2
//  p.m = 0.01;
  nlev = 2; // NUMBER OF LEVELS:  nlev = 0 give top level alone
  if(nlev  > p.Lmax){
    printf("  ERROR: More levels than available in lattice! \n");
    return 1;
  }

  double resmag = 1.0; // not rescaled.
  int ncycle = 0;
  int n_per_lev = 10;

  printf("  V cycle for %d by %d lattice with nlev = %d out of max %d \n", p.N, p.N, nlev, p.Lmax);

  // initialize arrays__________________________________
  p.size[0] = p.N;
  p.a[0] = 1.0;
  p.alpha = 1.0;
  for(lev = 1;lev< p.Lmax+1; lev++) {
    p.size[lev] = (p.size[lev-1])/2;
    p.a[lev] = 2.0 * p.a[lev-1];
    // p.scale[lev] = 1.0/(4.0 + p.m*p.m*p.a[lev]*p.a[lev]);
  }

  for(lev = 0;lev< p.Lmax+1; lev++) {
    phi[lev] = (double *) malloc(p.size[lev]*p.size[lev] * sizeof(double));
    res[lev] = (double *) malloc(p.size[lev]*p.size[lev] * sizeof(double));
    for(i = 0; i < p.size[lev]*p.size[lev]; i++) {
      phi[lev][i] = 0.0;
      res[lev][i] = 0.0;
    }
  }

  // set up the boundary conditions
  for (i=0; i<p.N; i++) {
    double tmp = (*boundary_func)(i);
    phi[0][i] = tmp;  // top edge, left to right
    // res[p.N * p.N - (p.N * (i+1))] = boundary_func(i);
    phi[0][p.N * (p.N - i - 1)] = tmp;  // left edge, bottom to top
    phi[0][(p.N * (i+1)) - 1] = tmp;    // right edge, top to bottom
    phi[0][(p.N * p.N) - i - 1] = tmp;   // bottom edge, right to left
  }

  FILE *nfile = fopen("res_data.dat", "w+");
  for (i=0; i<p.N; i++) {
    for (j=0; j<p.N; j++) {
      fprintf(nfile, "%i %i %f\n", i, j, res[0][i + j*p.N]);
    }
  }

  // set up the heat source
  res[0][(1*(p.N/4)) + ((1*(p.N/4)))*p.N] = p.alpha;
  res[0][(1*(p.N/4)) + ((3*(p.N/4)))*p.N] = p.alpha;
  res[0][(3*(p.N/4)) + ((1*(p.N/4)))*p.N] = p.alpha;
  res[0][(3*(p.N/4)) + ((3*(p.N/4)))*p.N] = p.alpha;

  // iterate to solve_____________________________________
  resmag = 1.0; // not rescaled.
  ncycle = 0;
  n_per_lev = 10;
  resmag = GetResRoot_parallel(phi[0],res[0],0,p);
  printf("    At the %d cycle the mag residue is %g \n", ncycle, resmag);

  // while(resmag > 0.00001 && ncycle < 10000)
  while(resmag > 0.00001) {
    ncycle +=1;
    for(lev = 0;lev<nlev; lev++) {  //go down
      relax_parallel(phi[lev],res[lev],lev, n_per_lev,p); // lev = 1, ..., nlev-1
      proj_res_parallel(res[lev + 1], res[lev], phi[lev], lev,p);    // res[lev+1] += P^dag res[lev]
    }

    for(lev = nlev;lev >= 0; lev--) { //come up
      relax_parallel(phi[lev],res[lev],lev, n_per_lev,p);   // lev = nlev -1, ... 0;
      if(lev > 0) inter_add_parallel(phi[lev-1], phi[lev], lev, p);   // phi[lev-1] += error = P phi[lev] and set phi[lev] = 0;
    }
    resmag = GetResRoot_parallel(phi[0],res[0],0,p);
    if (ncycle % 100 == 0) printf("    At the %d cycle the mag residue is %g \n", ncycle, resmag);
  }

  FILE *file = fopen("parallel_data.dat", "w+");
  for (i=0; i<p.N; i++) {
    for (j=0; j<p.N; j++) {
      fprintf(file, "%i %i %f\n", i, j, phi[0][i + j*p.N]);
    }
  }

  return 0;
}

void relax_parallel(double *T, double *b, int lev, int niter, param_t p) {
  // printf("relax2: level %i\n", lev);
  int i, j;
  int L = p.size[lev];
  int n=L; int m=L;
  double alpha=p.alpha;

  for(int iter=0; iter < niter; iter++) {
      #pragma omp parallel for shared(T, b) private(j)
      for (i=1; i<n-1; i++) {
          for (j=1; j<m-1; j+=2) {
              int shift;
              if (i%2 == 0) { shift = 1; }
              else { shift = 0; }
              int index       = getIndex(i, j+shift, m);
              int pos_l  = getIndex(i-1, j+shift, m);
              int pos_r = getIndex(i+1, j+shift, m);
              int pos_u    = getIndex(i, j+1+shift, m);
              int pos_d  = getIndex(i, j-1+shift, m);

              T[index] = (1.0 - alpha) * T[index] +
                         ((alpha/4.0) *
                         (T[pos_l] + T[pos_r] +
                          T[pos_u]   + T[pos_d] )) + b[index];
          }
      }

      // update second half of the board
      #pragma omp parallel for shared(T, b) private(j)
      for (i=1; i<n-1; i++) {
          for (j=2; j<m-1; j+=2) {
              int shift;
              if (i%2 == 0) { shift = -1; }
              else { shift = 0; }
              int index       = getIndex(i, j+shift, m);
              int pos_l  = getIndex(i-1, j+shift, m);
              int pos_r = getIndex(i+1, j+shift, m);
              int pos_u    = getIndex(i, j+1+shift, m);
              int pos_d  = getIndex(i, j-1+shift, m);

              T[index] = (1.0 - alpha) * T[index] +
                         ((alpha/4.0) *
                         (T[pos_l] + T[pos_r] +
                          T[pos_u]   + T[pos_d] )) + b[index];
          }
      }
  }
  return;
}

void proj_res_parallel(double *res_c, double *res_f, double *phi_f, int lev, param_t p)
{
  int L, Lc, i, j;
  L = p.size[lev];
  double r[L*L]; // temp residue
  Lc = p.size[lev+1];  // course level

  //get residue
  #pragma omp parallel for private(j) shared(r, phi_f, res_f, p, L)
  for(i = 1; i< L-1; i++) {
    for(j = 1; j< L-1; j++) {
        int index       = getIndex(i, j, L);
        int pos_l  = getIndex(i-1, j, L);
        int pos_r = getIndex(i+1, j, L);
        int pos_u    = getIndex(i, j+1, L);
        int pos_d  = getIndex(i, j-1, L);

        r[index] = (-1.0*p.alpha*phi_f[index]) +
            ((p.alpha/4.0) *
            (phi_f[pos_l] + phi_f[pos_r] +
             phi_f[pos_u]   + phi_f[pos_d] )) + res_f[index];
      /*
      left  = (x == 0)   ? res_f[    y*L] : phi_f[(x-1) +  y   *L];
      right = (x == L-1) ? res_f[x + y*L] : phi_f[(x+1) +  y   *L];
      up    = (y == 0)   ? res_f[x      ] : phi_f[ x    + (y-1)*L];
      down  = (y == L-1) ? res_f[x + y*L] : phi_f[ x    + (y+1)*L];
      r[x + y*L] = res_f[x + y*L] -  phi_f[x + y*L] + p.scale[lev]*(left + right + up + down);
      */

    }
  }

  //project residue
  int x, y;
  #pragma omp parallel for private(x, y) shared(p, L, Lc)
  for(x = 1; x < Lc-1; x++) {
    for(y = 1; y < Lc-1; y++) {
      // printf("  project res %i,%i\n", x, y);
      res_c[x + y*Lc] = 0.25*(r[ 2*x      +  2*y   *L] +
                              r[(2*x + 1) +  2*y   *L] +
                              r[ 2*x      + (2*y+1)*L] +
                              r[(2*x+1)%L + (2*y+1)*L]);
    }
  }
  return;
}

void inter_add_parallel(double *phi_f,double *phi_c,int lev,param_t p)
{
  int L, Lc, x, y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1];

  #pragma omp parallel for private(x, y) shared(p, L, Lc)
  for(x = 0; x < Lc; x++) {
    for(y = 0; y < Lc; y++) {
    	phi_f[ 2*x      +  2*y   *L] += phi_c[x + y*Lc];
    	phi_f[(2*x + 1) +  2*y   *L] += phi_c[x + y*Lc];
    	phi_f[ 2*x      + (2*y+1)*L] += phi_c[x + y*Lc];
    	phi_f[(2*x + 1) + (2*y+1)*L] += phi_c[x + y*Lc];
    }
  }

  //set to zero so phi = error
  #pragma omp parallel for private(x, y) shared(p, L, Lc)
  for(x = 0; x < Lc; x++) {
    for(y = 0; y < Lc; y++) {
      phi_c[x + y*Lc] = 0.0;
    }
  }
  return;
}

double GetResRoot_parallel(double *phi, double *res, int lev, param_t p)
{
  //true residue
  int i, j, L = p.size[lev];
  double residual = 0.0;

  #pragma omp parallel for private(j) shared(phi, res, p, L) \
  reduction(+:residual)
  for (i=1; i<L-1; i++) {
      for (j=1; j<L-1; j++) {
          double r;
          int index       = getIndex(i, j, L);
          int pos_l  = getIndex(i-1, j, L);
          int pos_r = getIndex(i+1, j, L);
          int pos_u    = getIndex(i, j+1, L);
          int pos_d  = getIndex(i, j-1, L);

          r = (-1.0*p.alpha*phi[index]) +
              ((p.alpha/4.0) *
              (phi[pos_l] + phi[pos_r] +
               phi[pos_u]   + phi[pos_d] )) + res[index];
          residual += r*r;
      }
  }
/*
  for(x = 0; x < L; x++) {
    for(y = 0; y<L; y++) {
      left  = (x == 0)   ? res[    y*L] : phi[(x-1) +  y   *L];
      right = (x == L-1) ? res[x + y*L] : phi[(x+1) +  y   *L];
      up    = (y == 0)   ? res[x      ] : phi[ x    + (y-1)*L];
      down  = (y == L-1) ? res[x + y*L] : phi[ x    + (y+1)*L];
      residue = res[x + y*L]/p.scale[lev] - phi[x + y*L]/p.scale[lev] + (left + right + up + down);
      ResRoot += residue*residue; // true residue
    }
  }
*/
  return sqrt(residual);
}
