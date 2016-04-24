#include "redblack.h"
#include <iostream>
#include <cmath>
using namespace std;


int main() {
    const int n = 100;
    const int m = 100;
    double alpha = 0.9;
    double * T = new double[n*m];
    double * b = new double[n*m];
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (i==0 || i==n-1) { T[getIndex(i, j, m)] = 0.0; }
			else if (j==0 || j==m-1) { T[getIndex(i, j, m)] = 20.0; }
            else { T[getIndex(i, j, m)] = 0.0; }
            b[getIndex(i, j, m)] = 0.0;
        }
    }
    b[getIndex(n/2, m/2, m)] = 1.0;

    double cc = 99.0;
    int iter = 0;
    while (cc > 1e-5) {
        redBlackSerial(T, b, n, m, alpha, 1);
        cc = computeResidual(T, b, n, m, alpha);
        if (iter%1000==0) {
            cout << "Iteration: " << iter << ", Res: " << cc << "\n";
        }
        iter ++;

    }

    FILE * output = fopen("T.txt", "w");
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            fprintf(output, "%i %i %0.16g\n", i, j, T[getIndex(i, j, m)]);
        }
    }



    return 0;
}
