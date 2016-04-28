#include <iostream>
#include <cmath>
#include "redblack.h"
#include <omp.h>
using namespace std;
#define TOL 1.0e-6

/* Red-Black 2D algorithm without parallelization

    In: T, b, n, m, itermax (optional)
        n,m are # rows / colummns
*/
void redBlackSerial(double * T, double * b,
                    int n, int m, double alpha, int itermax) {
    int i, j;
    int index, pos_l, pos_r, pos_u, pos_d;

    int iter = 0;
    while (iter < itermax) {
        // update half the board
        for (i=1; i<n-1; i+=2) {
            for (j=1; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                pos_l  = getIndex(i-1, j, m);
                pos_r = getIndex(i+1, j, m);
                pos_u    = getIndex(i, j+1, m);
                pos_d  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[pos_l] + T[pos_r] +
                            T[pos_u]   + T[pos_d] )) + b[index];
            }
        }
        for (i=2; i<n-1; i+=2) {
            for (j=2; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                pos_l  = getIndex(i-1, j, m);
                pos_r = getIndex(i+1, j, m);
                pos_u    = getIndex(i, j+1, m);
                pos_d  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[pos_l] + T[pos_r] +
                            T[pos_u]   + T[pos_d] )) + b[index];
            }
        }

        // update second half of the board
        for (i=1; i<n-1; i+=2) {
            for(j=2; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                pos_l  = getIndex(i-1, j, m);
                pos_r = getIndex(i+1, j, m);
                pos_u    = getIndex(i, j+1, m);
                pos_d  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[pos_l] + T[pos_r] +
                            T[pos_u]   + T[pos_d] )) + b[index];
            }
        }
        for (i=2; i<n-1; i+=2) {
            for (j=1; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                pos_l  = getIndex(i-1, j, m);
                pos_r = getIndex(i+1, j, m);
                pos_u    = getIndex(i, j+1, m);
                pos_d  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[pos_l] + T[pos_r] +
                            T[pos_u]   + T[pos_d] )) + b[index];
            }
        }
        iter ++;
    }
}
/* Red-Black 2D algorithm with parallelization

    In: T, b, n, m, itermax (optional)
        n,m are # rows / colummns
*/
void redBlackParallel(double * T, double * b,
                    int n, int m, double alpha, int itermax) {
    int i, j;
    //int index, pos_l, pos_r, pos_u, pos_d, shift;

    int iter = 0;
    while (iter < itermax) {
        // update half the board in parallel
	#pragma omp parallel for collapse(2) shared(T, b)
        for (i=1; i<n-1; i++) {
           // #pragma omp parallel for shared(T, b)
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
        #pragma omp parallel for collapse(2) shared(T,b)
        for (i=1; i<n-1; i++) {
	    //#pragma omp parallel for shared(T, b)
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
        iter ++;
    }

}

void gaussSeidelParallel(double* T, double* b, int n, int m, double alpha, int itermax) {
    int iter=0;
    int i, j;
    while (iter < itermax) {
	#pragma omp parallel for shared(T,b) private(j)
        for (i=1; i<n-1; i++) {
	    //#pragma omp parallel for shared(T, b)
            for (j=1; j<m-1; j++) {
                int index       = getIndex(i, j, m);
                int pos_l  = getIndex(i, j-1, m);
                int pos_r = getIndex(i, j+1, m);
                int pos_u    = getIndex(i+1, j, m);
                int pos_d  = getIndex(i-1, j, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[pos_l] + T[pos_r] +
                            T[pos_u]   + T[pos_d] )) + b[index];
            }
        }
	iter++;
    }
}

int getIndex(int x, int y, int width) { return x*width + y; }

double computeNorm(double * x, int N) {
    double norm = 0.0;
    for (int i=0; i<N; i++) { norm = norm + (x[i]*x[i]); }
    return sqrt(norm);
}

double computeResidual(double * T, double * b, int n, int m, double alpha) {
    int i, j;
    int index, pos_l, pos_r, pos_u, pos_d;

    double residual = 0.0;
    for (i=1; i<n-1; i++) {
        for (j=1; j<m-1; j++) {
            double r;
            index       = getIndex(i, j, m);
            pos_l  = getIndex(i-1, j, m);
            pos_r = getIndex(i+1, j, m);
            pos_u    = getIndex(i, j+1, m);
            pos_d  = getIndex(i, j-1, m);

            r = (-1.0*alpha*T[index]) +
                ((alpha/4.0) *
                (T[pos_l] + T[pos_r] +
                 T[pos_u]   + T[pos_d] )) + b[index];
            residual += r*r;
        }
    }
    return sqrt(residual);

}
