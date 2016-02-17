#include <math.h>
#include <tgmath.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "uqtktools.h"
#include "uqtkmcmc.h"
#include "PCBasis.h"
#include "PCSet.h"
#include "arraytools.h"
#include "getopt.h"

#include "Utils.h"
#include "Utilsave.h"
#include "Duffing.h"
#include "KL.h"
#include "MCS.h"
#include "nGhanemSpanos.h"
#include "AAPG.h"
#include "ticktock.h"

int main(int argc, char *argv[]){

    int dof=10;
    int ord_GS=2;
    int dim=5;
    int noutput=2;
    string pcType="LU";  //PC type
    Array1D<Array1D<double> > initial(dof); // Initial conditions

    // epsilon
    Array1D<Array1D<double> > epsilon(dof);

    // Time marching info
    double dTym = 0.01;
    double tf = 10;
    // Number of steps
    int nStep=(int) tf / dTym;

    // MCK
    Array1D<Array1D<double> > mck(3);
    // m
    Array1D<double> temp_m(dof,1.e4);
    mck(0) = temp_m;
    // k
    Array1D<double> temp_k(dof,4.e7);
    mck(2) = temp_k;
    // c
    Array1D<double> temp_c(dof,2*0.04*sqrt(1.e4*4.e7));
    mck(1) = temp_c;

    // force
    Array1D<Array2D<double> > f_GS(dof);
    Array2D<double> scaledKLmodes(2*nStep+1,dim,0.e0);
    double clen = 0.05;
    double sigma=0.5;
    char* cov_type = (char *)"Exp";
    genKL(scaledKLmodes, 2*nStep+1, dim, clen, sigma, tf, cov_type);
    Array1D<double> fbar(2*nStep+1,0.e0);//mean of forcing
    double t_temp = 0.0; 
    for (int i=0;i<2*nStep+1;i++){
        fbar(i) = 2.0*(2.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
        t_temp +=dTym/2;
    }
    
    //Array1D<Array2D<double> > mstd_MCS(dof);
    cout << "Starting GS..." << endl;

    //Array2D<double> e_GS(ord_GS,)
    for (int ord=1;ord<=ord_GS;ord++){
        TickTock tt;
        tt.tick();
	    PCSet myPCSet("ISP",ord,dim,pcType,0.0,1.0); 
        tt.tock("Took");
        int nPCTerms = myPCSet.GetNumberPCTerms();
        // Prepare the force in PC format
        Array1D<Array2D<double> > f_GS(dof);
        Array2D<double> temp_f(2*nStep+1,nPCTerms,0.e0);
        for (int i=0;i<dof;i++)
            f_GS(i)=temp_f;
        f_GS(0).replaceCol(fbar,0);
        for (int i=0;i<dim;i++){
            Array1D<double> tempf(2*nStep+1,0.e0);
            getCol(scaledKLmodes,i,tempf);
            f_GS(0).replaceCol(tempf,i+1);
        }
        // initial conditions set to zero for now
        for (int i=0;i<dof;i++){
            Array1D<double> temp_init(2*nPCTerms,0.e0);
            initial(i)= temp_init;   
        }
        Array1D<Array2D<double> > uv_solution(dof);
        for (int i=0;i<dof;i++){
            Array2D<double> temp_solution(nStep,2*nPCTerms,0.e0);
            uv_solution(i) = temp_solution;
        }
        //epsilon
        for (int i=0;i<dof;i++){
            Array1D<double> temp_epsilon(nPCTerms,0.e0);
            epsilon(i) = temp_epsilon;
            epsilon(i)(0) = 0.1;
        }
        cout << "Starting GS order " << ord << endl;
        nGS(dof, myPCSet, epsilon, mck, nStep, initial, dTym, f_GS, uv_solution);
        cout << "Order " << ord << " finished." <<endl; 
        // Post-process the solution
        //Array2D<double> e_GS = postprocess_nGS(dof,noutput,nPCTerms,nStep,uv_solution,myPCSet,dTym,ord,mstd_MCS_u,mstd_MCS_v);
    }

    return 0;
}
