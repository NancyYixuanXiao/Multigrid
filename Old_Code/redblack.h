void redBlackSerial(double * T, double * b,
                   int n, int m, double alpha, int itermax);
void redBlackParallel(double * T, double * b,
                       int n, int m, double alpha, int itermax);
void gaussSeidelParallel(double* T, double* b, int n, int m, double alpha, int itermax);
int getIndex(int x, int y, int width);
double computeResidual(double * T, double * b, int n, int m, double alpha);
