#include "kldecompuni.h"
#include "Array1D.h"
#include "Array2D.h"
#include <math.h>
#include "KL.h"
#include "KL_exp.h"

    void comcov(Array2D<double> &cov, const int npts, Array1D<double> &xgrid, const double clen, const double sigma, const char* cov_type) {
        cov.Resize(npts,npts,0.e0);	
        for ( int i = 0; i < npts; i++)
            for ( int j = 0; j < npts; j++)
            cov(i,j)=genCovAnl(xgrid(i),xgrid(j),clen,sigma,cov_type);
    }

    void genGrid(Array1D<double> &xgrid, const int npts, const double t_final) {
        if (npts<=0)
            throw Tantrum("kl_sample::genGrid() : number of grid points needs to be greater than 0") ;
        xgrid.Resize(npts,0.0);
        for ( int i=0; i<npts; i++ ) xgrid(i) = (double) (i / ((double) npts - 1.0))*t_final;
            return;
    }
 
    double genCovAnl(const double x, const double y, const double clen, const double sigma, const string covtype) {
        double cov;
        if ( covtype == "SqExp" )
            cov = exp( - (x-y) * (x-y) / ( clen * clen ) ) * sigma * sigma ;
        else if ( covtype == "Exp" )
            cov = exp( - fabs(x-y) / clen ) * sigma * sigma ;
        else 
            throw Tantrum("kl_sample.cpp::genCovAnl(): covariance type is not recognized!");
        return ( cov );
    }

    void genKL(Array2D<double> &scaledKLmodes, const int npts, const int nkl, const double clen, const double sigma, const double t_final, const char* cov_type){
        Array1D<double> xgrid;
        Array2D<double> cov;
        genGrid(xgrid,npts,t_final);

        Array2D<double> KLmodes(npts,nkl,0.e0);
        string cov_type_s(cov_type);
        Array1D<double> eigs(nkl,0.e0);
        if (cov_type_s == "Exp"){
            eigs = KL_exp(sigma, t_final/2, clen, nkl, KLmodes, xgrid);
        }

        else{
        comcov(cov, npts, xgrid, clen, sigma, cov_type);

        KLDecompUni decomposer(xgrid);
        int n_eig = decomposer.decompose(cov,nkl);

        if(n_eig < nkl){
            printf("There are only %d eigenvalues available (requested %d) \n", n_eig, nkl);
        //    nkl = n_eig;    
        }

        eigs = decomposer.eigenvalues();
        KLmodes = decomposer.KLmodes();
        }

        for ( int i = 0; i < npts; i++ ){
            //scaledKLmodes(i,0) = xgrid(i);
            for ( int j = 0; j < nkl; j++ )
                scaledKLmodes(i,j) = KLmodes(i,j)*sqrt(eigs(j));
        }
    return;
    }
