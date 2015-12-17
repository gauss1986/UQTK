/*===================================================================================== 
Lorenz.cpp
by Lin Gao (lingao.gao@mail.utoronto.ca)
Dec 14, 2015
===================================================================================== */

#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include "uqtktools.h"
#include "uqtkmcmc.h"
#include "PCSet.h"
#include "arraytools.h"
#include "getopt.h"

#include "UtilsLorenz.h"
#include "Utils.h"
#include "Utilsave.h"
#include "Duffing.h"
#include "KL.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "AAPG.h"
#include "ticktock.h"

#define DIM 3
#define CLEN 0.1
#define SIG 0.0
#define ORDER_GS 2
#define ORDER_AAPG_GS 2
#define ORDER_AAPG 3
#define TFINAL 10.0
#define DTYM 0.01
#define NSPL 10000
#define AMP 0.0
#define FBAR 8.0
#define x0 1.0 // A point on or nearly on the attractor by Lorenz 2005
#define y0 0.0
#define z0 -0.75
#define FACTOR_OD 1.0
#define ACTD false
#define THRESHOLD 0.99

#define COVTYPE "Exp"
#define PCTYPE "HG"

/// \brief Displays information about this program
int usage(){
  printf("usage: cor_kl [-h] [-c<cov_type>] [-d<dTym>] [-n<dim>] [-p<nspl>] [-s<sigma>] [-l<clen>] [-t<tf>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -c <cov_type>    : define wether using analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -d <dTym>        : define the time marching step (default=%e) \n",DTYM);
  printf(" -n <dim>         : define the number of KL modes retained (default=%d) \n",DIM);
  printf(" -p <nspl>        : define the number of samples (default=%d) \n",NSPL);
  printf(" -s <sigma>       : define standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : define correlation length (default=%e) \n",CLEN);
  printf(" -t <tf>          : define end time (default=%e) \n",TFINAL);
  printf(" -a <AMP>         : define the oscillation amp of F (default=%e) \n",AMP);
  printf(" -F <FBAR>        : define mean of F(default=%e) \n",FBAR);
  printf(" -f <factor_OD>   : define the factor in overdetermined dynimcal-orthogonal calculation (default=%e) \n", FACTOR_OD);
  printf(" -G <ord_GS>      : define the order of Ghanem-Spanos method (default=%d) \n", ORDER_GS);
  printf(" -A <ord_AAPG_GS> : define the order of Ghanem-Spanos method used in subproblems of AAPG (default=%d) \n", ORDER_AAPG_GS);
  printf(" -P <ord_AAPG>    : define the order of AAPG method(default=%d) \n", ORDER_AAPG);
  printf(" -T <p>           : define the threshold used in the adaptive AAPG method(default=%e) \n", THRESHOLD);
  printf("================================================================================================================\n");
  printf("Input:: Nothing hard-coded.\n");
  printf("Output:: \n");
  printf("        - KLmodes.dat:  eigenmodes scaled with sqrt(eig)\n");
  printf("================================================================================================================\n");
  exit(0);
  return 0;
}

/// \brief c++ version of the Matlab code GS_duffing_multi.m.
/// Dynamic response of a Duffing oscillator subject to uncertain excitation force.
/// The unceratain force is generated using KL expansion.
/// Mass, stiffness, damping and nonlinear coefficients are assumed to be deterministic.

/// Main program of uncertainty propagation of the ODE model excitation force via intrusive spectral projection (ISP)
int main(int argc, char *argv[])
{
    // UQ specific info
    // PC type
    string pcType = PCTYPE;
    // Stochastic dimensionality
    int dim = DIM;
    // Covariance type
    char* cov_type = (char *)COVTYPE;
    // Standar deviation
    double sigma = SIG;
    // Correlation length
    double clen = CLEN;
    // MC sample size
    int nspl = NSPL;
    // dynamical-orthogonal info
    double factor_OD = FACTOR_OD;
    // GS order
    int ord_GS = ORDER_GS;
    // AAPG order
    int ord_AAPG = ORDER_AAPG;
    // GS order in AAPG subproblems
    int ord_AAPG_GS = ORDER_AAPG_GS;
    // wheather doing AAPG adaptive or not
    bool act_D = ACTD;
    // Threashold used in the adaptive AAPG scheme
    double p = THRESHOLD;
    // spatial dof
    double dof = 3;

    // Time marching info
    double dTym = DTYM;
    double tf = TFINAL;
  
    /* Read the user input */
    int c;

    while ((c=getopt(argc,(char **)argv,"h:c:d:n:p:s:l:t:f:G:A:P:D:T:"))!=-1){
        switch (c) {
        case 'h':
            usage();
            break;
        case 'c':
            cov_type = optarg;
            break;
        case 'd':
            dTym = strtod(optarg, (char **)NULL);
            break;
        case 'n':
            dim = strtod(optarg, (char **)NULL);
            break;
        case 'p':
            nspl = strtod(optarg, (char **)NULL);
            break;
        case 's':
            sigma = strtod(optarg, (char **)NULL);
            break;
        case 'l':
            clen = strtod(optarg, (char **)NULL);
            break;
        case 't':
            tf = strtod(optarg, (char **)NULL);
            break;
        case 'f':
            factor_OD = strtod(optarg, (char **)NULL);
            break;
        case 'G':
            ord_GS = strtod(optarg, (char **)NULL);
            break;
        case 'A':
            ord_AAPG_GS = strtod(optarg, (char **)NULL);
            break;
        case 'P':
            ord_AAPG = strtod(optarg, (char **)NULL);
            break;
        case 'T':
            p = strtod(optarg,(char **)NULL);
            break;
        case 'D':
            int temp_D = strtod(optarg, (char **)NULL);
            if (temp_D==1)
                act_D=true;
            break;
        }
    }
    
    /* Print the input information on screen */
    cout << "Lorenz model--" << endl<<flush;
    cout << " - Number of KL modes:              " << dim  << endl<<flush;
    cout << " - Monte Carlo sample size:         " << nspl  << endl<<flush;
    cout << " - Will generate covariance of type "<<cov_type<<", with correlation length "<<clen<<" and standard deviation "<<sigma<<endl<<flush;
    cout << " - Time marching step:              " << dTym  << endl<<flush;
    cout << " - Process end time:                " << tf  << endl<<flush;
    cout << " - Dynamical orthogonal AAPG factor:" << factor_OD  << endl<<flush;
    cout << " - The threshold in adaptive AAPG:  " << p  << endl<<flush;
    cout << " - Order of GS:                     " << ord_GS  << endl<<flush;
    cout << " - Order of GS in AAPG subproblems: " << ord_AAPG_GS  << endl<<flush;
    cout << " - Order of AAPG:                   " << ord_AAPG  << endl<<flush;
    if (act_D){
        cout << " - Adaptive AAPG applied"<<endl<<flush;
    }
    else{
        cout << " - Full AAPG applied"<<endl<<flush;
    }
    if (!strcmp(pcType.c_str(),"LU")){
        cout << " - Unifrom random variable used."<<endl<<flush;
    }
    else{
        if (!strcmp(pcType.c_str(),"HG"))
        cout << " - Gaussian random variable used"<<endl<<flush;
    }


    // Number of steps
    int nStep=(int) tf / dTym;

    // Generate the KL expansion of the excitation force
    int nkl = dim;
    Array2D<double> scaledKLmodes(2*nStep+1,nkl,0.0);
    genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
    write_datafile(scaledKLmodes,"Lorenz/KL.dat");
 
    // Monte Carlo simulation
    // Define input parameters
    Array1D<double> inpParams(6,0.e0);
    inpParams(0) = 1;       // code for problem: 0-Duffing, 1-Lorenz
    inpParams(1) = 0.25;    //a, from Lorenz 2005
    inpParams(2) = 4.0;     //b, from Lorenz 2005
    inpParams(3) = 1.23;    //G, from Lorenz 2005
    inpParams(4) = AMP;     //AMP
    inpParams(5) = 2*M_PI/73; //w
    double fbar = FBAR;
    // Draw Monte Carlo samples
    Array2D<double> samPts(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",ord_GS,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts);
    Array1D<double> samPts_1D(nspl,0.e0);
    getCol(samPts,0,samPts_1D);
    Array1D<double> sample_mstd = mStd(samPts_1D,nspl);
    cout << "Mean of sample is "<< sample_mstd(0) << endl;
    cout << "Std of sample is "<< sample_mstd(1) << endl;
    // Open files to write out
    string Solufile_mean = "Lorenz/MCS_mean.dat";
    string Solufile_std = "Lorenz/MCS_std.dat";
    FILE *f_mean = createfile(Solufile_mean);
    FILE *f_std  = createfile(Solufile_std);

    // Initialize solutions
    Array2D<double> x(nStep+1,nspl,0.e0);
    Array2D<double> y(nStep+1,nspl,0.e0);
    Array2D<double> z(nStep+1,nspl,0.e0);
    Array1D<Array2D<double> > result(3);
    result(0) = x;
    result(1) = y;
    result(2) = z;
    // Define and write solutions at initial step
    Array1D<double> initial(dof,0.e0);  // deterministic initial condition
    initial(0) = x0;
    initial(1) = y0;
    initial(2) = z0;
    WriteMeanStdDevToFilePtr_lorenz(0, initial(0), initial(1), initial(2), f_mean);        
    WriteMeanStdDevToFilePtr_lorenz(0, 0, 0, 0, f_std);        
    cout << "\nMCS...\n" << endl;  
    double t_MCS = MCS(dof, nspl, dim, nStep, nkl, dTym, fbar, scaledKLmodes, inpParams, samPts, result, initial);
    Array2D<double> mean(dof,nStep+1,0.e0);
    Array2D<double> std(dof,nStep+1,0.e0);
    Array1D<Array2D<double> > mstd_MCS(dof);
    // post-process the solution
    for (int ivar=0;ivar<dof;ivar++){
        Array2D<double> temp(2,nStep+1,0.e0);
        mstd_MCS(ivar) = temp;
        for (int ix=1;ix<nStep+1;ix++){
            Array1D<double> temp(nspl,0.e0);
            getRow(result(ivar),ix,temp);
            Array1D<double> mstd(2,0.e0);
            mstd = mStd(temp,nspl);
            mstd_MCS(ivar).replaceCol(mstd,ix); 
            mean(ivar,ix) = mstd(0);
            std(ivar,ix) = mstd(1);
        }
    }
    // save results
    for (int ix=1;ix<nStep+1;ix++){
        WriteMeanStdDevToFilePtr_lorenz((ix+1)*dTym, mean(0,ix),mean(1,ix),mean(2,ix),f_mean);         
        WriteMeanStdDevToFilePtr_lorenz((ix+1)*dTym, std(0,ix),std(1,ix),std(2,ix),f_std);         
    }
    cout << "Results of mean..." << endl;
    // Output mean to scren
    for (int ix=0;ix<nStep;ix++){
        if ((ix+1) % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut_lorenz(ix+1, (ix+1)*dTym, mean(0,ix), mean(1,ix), mean(2,ix));
        }
    }
    cout << "Results of std..." << endl;
    // Output std to scren
    for (int ix=0;ix<nStep;ix++){
        if ((ix+1) % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut_lorenz(ix+1, (ix+1)*dTym, std(0,ix), std(1,ix), std(2,ix));
        }
    }
    closefile(f_mean,Solufile_mean);
    closefile(f_std,Solufile_std);
    

    // Ghanem-Spanos method
    printf("\nGhanem-Spanos method...\n");

    Array1D<double> t_GS(ord_GS,0.e0);
    
    ostringstream err_stream;
    if (act_D){
        err_stream << "error_n" << dim <<"_actD.dat";
    }
    else {
        err_stream << "error_n" << dim <<".dat";
    }
    string errstr(err_stream.str());
    FILE *err_dump=createfile(errstr);
    for(int ord=1;ord<ord_GS+1;ord++){
        TickTock tt;
    	tt.tick();
	    PCSet myPCSet("ISP",ord,dim,pcType,0.0,1.0); 
            tt.tock("Took");
	    cout << "Order "<< ord << endl;

	    // The number of PC terms
        const int nPCTerms = myPCSet.GetNumberPCTerms();
        cout << "The number of PC terms in an expansion is " << nPCTerms << endl;
        // Print the multiindices on screen
    
        // Prepare the force in PC format
        Array2D<double> f_GS(2*nStep+1,nPCTerms,0.e0);
        for (int i=0;i<2*nStep+1;i++){
            f_GS(i,0) = fbar+inpParams(4)*cos(inpParams(5)*i);;
        }
        for (int i=0;i<nkl;i++){
            Array1D<double> tempf(2*nStep+1,0.e0);
            getCol(scaledKLmodes,i,tempf);
            f_GS.replaceCol(tempf,i+1);
        }

        // Allocate result variable
        Array1D<Array2D<double> > result(3);

        clock_t start = clock();
    	tt.tick();
        GS(dof, myPCSet, ord, dim, nPCTerms, pcType, nStep, initial, dTym, inpParams, f_GS, result);
        t_GS(ord-1) = tt.silent_tock();
	    cout << "Cost time: "<<(clock()-start)/(double)(CLOCKS_PER_SEC)<<endl;
	
        // output and save the solution
        // Open files to write out solutions
        ostringstream s1;
        s1 << "Lorenz/x_GS_modes_" << ord<<".dat";
        string SoluGSmodes_x(s1.str());
        FILE *GS_dump_x = createfile(SoluGSmodes_x);
        ostringstream s2;
        s2 << "Lorenz/y_GS_modes_" << ord<<".dat";
        string SoluGSmodes_y(s2.str());
        FILE *GS_dump_y = createfile(SoluGSmodes_y);
        ostringstream s3;
        s3 << "Lorenz/z_GS_modes_" << ord<<".dat";
        string SoluGSmodes_z(s3.str());
        FILE *GS_dump_z = createfile(SoluGSmodes_z);
        // Open files to write out statistics of the solutions
        ostringstream sx;
        sx << "Lorenz/x_GS_stat_" << ord<<".dat";
        string SoluGSstat_x(sx.str());
        FILE *GSstat_dump_x = createfile(SoluGSstat_x);
        ostringstream sy;
        sy << "Lorenz/y_GS_stat_" << ord<<".dat";
        string SoluGSstat_y(sy.str());
        FILE *GSstat_dump_y = createfile(SoluGSstat_y);
        ostringstream sz;
        sz << "Lorenz/z_GS_stat_" << ord<<".dat";
        string SoluGSstat_z(sz.str());
        FILE *GSstat_dump_z = createfile(SoluGSstat_z);
        
        Array1D<double> e_GS_ord_x = postprocess_GS(nPCTerms, nStep, x0, result(0), myPCSet, dTym, GS_dump_x, GSstat_dump_x, mstd_MCS(0));
        Array1D<double> e_GS_ord_y = postprocess_GS(nPCTerms, nStep, y0, result(1), myPCSet, dTym, GS_dump_y, GSstat_dump_y, mstd_MCS(1));
        Array1D<double> e_GS_ord_z = postprocess_GS(nPCTerms, nStep, z0, result(2), myPCSet, dTym, GS_dump_z, GSstat_dump_z, mstd_MCS(2));

        fprintf(err_dump, "%lg %lg", e_GS_ord_x(0),e_GS_ord_x(1));
        fprintf(err_dump, "\n");

        closefile(GS_dump_x,SoluGSmodes_x);
        closefile(GS_dump_y,SoluGSmodes_y);
        closefile(GS_dump_z,SoluGSmodes_z);
        closefile(GSstat_dump_x,SoluGSstat_x);
        closefile(GSstat_dump_y,SoluGSstat_x);
        closefile(GSstat_dump_z,SoluGSstat_x);

    }
    return 0;
}
