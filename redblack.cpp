#include <iostream>
#include <cmath>
#include "redblack.h"
using namespace std;
#define TOL 1.0e-6

/* Red-Black 2D algorithm without parallelization

    In: T, b, n, m, itermax (optional)
        n,m are # rows / colummns
*/
void redBlackSerial(double * T, double * b,
                    int n, int m, double alpha, int itermax) {
    int i, j;
    int index, index_left, index_right, index_up, index_down;

    int iter = 0;
    while (iter < itermax) {
        // update half the board
        for (i=1; i<n-1; i+=2) {
            for (j=1; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                index_left  = getIndex(i-1, j, m);
                index_right = getIndex(i+1, j, m);
                index_up    = getIndex(i, j+1, m);
                index_down  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[index_left] + T[index_right] +
                            T[index_up]   + T[index_down] )) + b[index];
            }
        }
        for (i=2; i<n-1; i+=2) {
            for (j=2; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                index_left  = getIndex(i-1, j, m);
                index_right = getIndex(i+1, j, m);
                index_up    = getIndex(i, j+1, m);
                index_down  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[index_left] + T[index_right] +
                            T[index_up]   + T[index_down] )) + b[index];
            }
        }

        // update second half of the board
        for (i=1; i<n-1; i+=2) {
            for(j=2; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                index_left  = getIndex(i-1, j, m);
                index_right = getIndex(i+1, j, m);
                index_up    = getIndex(i, j+1, m);
                index_down  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[index_left] + T[index_right] +
                            T[index_up]   + T[index_down] )) + b[index];
            }
        }
        for (i=2; i<n-1; i+=2) {
            for (j=1; j<m-1; j+=2) {
                index       = getIndex(i, j, m);
                index_left  = getIndex(i-1, j, m);
                index_right = getIndex(i+1, j, m);
                index_up    = getIndex(i, j+1, m);
                index_down  = getIndex(i, j-1, m);

                T[index] = (1.0 - alpha) * T[index] +
                           ((alpha/4.0) *
                           (T[index_left] + T[index_right] +
                            T[index_up]   + T[index_down] )) + b[index];
            }
        }
        iter ++;
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
    int index, index_left, index_right, index_up, index_down;

    double residual = 0.0;
    for (i=1; i<n-1; i++) {
        for (j=1; j<m-1; j++) {
            double r;
            index       = getIndex(i, j, m);
            index_left  = getIndex(i-1, j, m);
            index_right = getIndex(i+1, j, m);
            index_up    = getIndex(i, j+1, m);
            index_down  = getIndex(i, j-1, m);

            r = (-1.0*alpha*T[index]) +
                ((alpha/4.0) *
                (T[index_left] + T[index_right] +
                 T[index_up]   + T[index_down] )) + b[index];
            residual += r*r;
        }
    }
    return sqrt(residual);

}
