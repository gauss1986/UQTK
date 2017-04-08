#ifndef DUFFINGISP
#define DUFFINGISP

void comcov(Array2D<double> &cov, int npts, Array1D<double> &xgrid, const double clen, const double sigma, char* cov_type); 
void genGrid(Array1D<double> &xgrid, const int npts, const double t_final); 
double genCovAnl(const double x, const double y, const double clen, const double sigma, const string covtype);
void genKL(Array2D<double> &scaledKLmodes, const int npts, const int nkl, const double clen, const double sigma, const double t_final, const char* cov_type);
Array2D<int> computeBasis(int nvars,int order);
Array1D<int> index1(int dim, int i, int order);


#endif
