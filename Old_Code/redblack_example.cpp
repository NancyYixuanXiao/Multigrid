#include "redblack.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <omp.h>
using namespace std;

void initArrays(double * T, double * b, int n, int m) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (i==0 || i==n-1) { T[getIndex(i, j, m)] = 0.0; }
            else if (j==0 || j==m-1) { T[getIndex(i, j, m)] = 20.0; }
            else { T[getIndex(i, j, m)] = 0.0; }
            b[getIndex(i, j, m)] = 0.0;
        }
    }
    b[getIndex(n/2, m/2, m)] = 1.0;
}

int main() {
    const int n = 500;
    const int m = 500;
    double alpha = 0.9;
    double * T = new double[n*m];
    double * b = new double[n*m];
    initArrays(T, b, n, m);

    /* see if threads are utilized properly */
    int tid, nthreads;
    #pragma omp parallel shared(nthreads) private(tid)
    {
	tid = omp_get_thread_num();
	if (tid == 0) { nthreads = omp_get_num_threads(); }
	#pragma omp barrier
	printf("hello from thread %d of %d!\n", tid, nthreads);
    }

    double cc = 99.0;
    int iter = 0;
    cout << "Starting without parallel!\n";
    double t1, t2;
    double duration;

    t1 = omp_get_wtime();
    while (cc > 1e-5) {
        redBlackSerial(T, b, n, m, alpha, 1);
        cc = computeResidual(T, b, n, m, alpha);
        if (iter%1000==0) {
            cout << "Iteration: " << iter << ", Res: " << cc << "\n";
        }
        iter ++;

        //remove this
        cc=1e-6;
    }
    t2 = omp_get_wtime();
    duration = (t2 - t1); /// (double) CLOCKS_PER_SEC;
    cout << "Non-parallel: " << duration << " seconds.\n";

    FILE * output = fopen("T_serial.txt", "w");
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            fprintf(output, "%i %i %0.16g\n", i, j, T[getIndex(i, j, m)]);
        }
    }
    fclose(output);
    initArrays(T, b, n, m);

    cc = 99.0;
    iter = 0;
    cout << "Starting with parallel!\n";
    t1 = omp_get_wtime();
    while (cc > 1e-5) {
        //redBlackParallel(T, b, n, m, alpha, 1);
	gaussSeidelParallel(T, b, n, m, alpha, 1);
        cc = computeResidual(T, b, n, m, alpha);
        if (iter%1000==0) {
            cout << "Iteration: " << iter << ", Res: " << cc << "\n";
        }
        iter ++;
    }
    t2 = omp_get_wtime();
    duration = (t2 - t1);// / (double) CLOCKS_PER_SEC;
    cout << "Parallel: " << duration << " seconds.\n";

    output = fopen("T_parallel.txt", "w");
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            fprintf(output, "%i %i %0.16g\n", i, j, T[getIndex(i, j, m)]);
        }
    }
    fclose(output);



    return 0;
}
