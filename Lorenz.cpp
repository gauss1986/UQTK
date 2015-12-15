/*===================================================================================== 
DuffingISP.cpp
Modified based on the surf_rxn example in UQTK v.2.1.1
by Lin Gao (lingao.gao@mail.utoronto.ca)
July 25, 2015
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
#include "Duffing.h"
#include "KL.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "AAPG.h"
#include "ticktock.h"

#define DIM 10
#define CLEN 0.1
#define SIG 0.2
#define ORDER_GS 2
#define ORDER_AAPG_GS 2
#define ORDER_AAPG 3
#define TFINAL 10.0
#define DTYM 0.01
#define NSPL 100000
#define AMP 2.0
#define FBAR 7.0
#define x0 1.0
#define y0 1.0
#define z0 1.0
#define FACTOR_OD 1.0
#define ACTD false
#define THRESHOLD 0.99

#define COVTYPE "Exp"
#define PCTYPE "LU"

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

    // Number of steps
    int nStep=(int) tf / dTym;

    // Generate the KL expansion of the excitation force
    int nkl = dim;
    Array2D<double> scaledKLmodes(nStep+1,nkl,0.0);
    genKL(scaledKLmodes, nStep+1, nkl, clen, sigma, tf, cov_type);
    write_datafile(scaledKLmodes,"./Lorenz/KL.dat");
 
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
    // Open files to write out
    string Solufile_mean = "MCS_mean.dat";
    string Solufile_std = "MCS_std.dat";
    FILE *f_mean;
    FILE *f_std;
    if(!(f_mean=fopen(Solufile_mean.c_str(),"w"))){
        printf("Could not open file '%s'\n",Solufile_mean.c_str());
        exit(1);    
    }
    if(!(f_std=fopen(Solufile_std.c_str(),"w"))){
        printf("Could not open file '%s'\n",Solufile_std.c_str());
        exit(1);    
    }
    // Initialize solutions
    Array2D<double> x(nStep+1,nspl,0.e0);
    Array2D<double> y(nStep+1,nspl,0.e0);
    Array2D<double> z(nStep+1,nspl,0.e0);
    Array1D<Array2D<double> > result(3);
    result(0) = x;
    result(1) = y;
    result(2) = z;
    // Define and write solutions at initial step
    Array1D<double> initial(dof,0.e0);  // deterministic and zero initial condition
    WriteMeanStdDevToFilePtr_lorenz(0, initial(0), initial(1), initial(2), f_mean);        
    cout << "\nMCS...\n" << endl;  
    double t_MCS = MCS(dof, nspl, dim, nStep, nkl, dTym, fbar, scaledKLmodes, inpParams, samPts, result, initial);
    Array2D<double> mean(dof,nStep+1,0.e0);
    Array2D<double> std(dof,nStep+1,0.e0);
    // post-process the solution
    for (int ivar=0;ivar<dof;ivar++){
        for (int ix=0;ix<nStep;ix++){
            Array1D<double> temp(nspl,0.e0);
            getRow(result(ivar),ix,temp);
            Array1D<double> mstd(2,0.e0);
            mstd = mStd(temp,nspl);
            mean(ivar,ix) = mstd(0);
            std(ivar,ix) = mstd(1);
        }
    }
    // save results
    for (int ix=0;ix<nStep;ix++){
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
    if(fclose(f_mean)){
        printf("Could not close file '%s'\n",Solufile_mean.c_str());
        exit(1);    
    }
    if(fclose(f_std)){
        printf("Could not close file '%s'\n",Solufile_std.c_str());
        exit(1);    
    }
    
    return 0;
}

    // Ghanem-Spanos method
    printf("\nGhanem-Spanos method...\n");

    double dis0 = DIS0;
    double vel0 = VEL0;
    Array1D<double> t_GS(ord_GS,0.e0);
    
    ostringstream err_stream;
    if (act_D){
        err_stream << "error_n" << dim << "_e"<<epsilon<<"_s"<<sigma<<"_actD.dat";
    }
    else {
        err_stream << "error_n" << dim << "_e"<<epsilon<<"_s"<<sigma<<".dat";
    }
    string errstr(err_stream.str());
    FILE *err_dump;
    if(!(err_dump=fopen(errstr.c_str(),"w"))){
        printf("Could not open file '%s'\n",errstr.c_str());
        exit(1);    
    }
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
            f_GS(i,0) = fbar+inpParams(4)*cos(inpParams(5)*it);;
        }
        for (int i=0;i<nkl;i++){
            Array1D<double> tempf(2*nStep+1,0.e0);
            getCol(scaledKLmodes,i,tempf);
            f_GS.replaceCol(tempf,i+1);
        }

        // Assumed deterministic initial conditions
        Array2D<double> dis_GS(nStep+1,nPCTerms);

        clock_t start = clock();
    	tt.tick();
        dis_GS = GS(myPCSet, ord, dim, nPCTerms, pcType, nStep, dis0, vel0, dTym, inpParams, f_GS);
        t_GS(ord-1) = tt.silent_tock();
	    cout << "Cost time: "<<(clock()-start)/(double)(CLOCKS_PER_SEC)<<endl;
	
        // output and save the solution
        // Open files to write out solutions
        ostringstream s;
        s << "GS_modes_" << ord<<".dat";
        string SoluGSmodes(s.str());
        FILE *GS_dump;
        if(!(GS_dump=fopen(SoluGSmodes.c_str(),"w"))){
            printf("Could not open file '%s'\n",SoluGSmodes.c_str());
            exit(1);    
        }
        // Open files to write out statistics of the solutions
        ostringstream s1;
        s1 << "GS_stat_" << ord<<".dat";
        string SoluGSstat(s1.str());
        //string SoluGSstat = "GS_stat.dat";
        FILE *GSstat_dump;
        if(!(GSstat_dump=fopen(SoluGSstat.c_str(),"w"))){
            printf("Could not open file '%s'\n",SoluGSstat.c_str());
            exit(1);    
        }
        
        Array1D<double> e_GS_ord = postprocess_GS(nPCTerms, nStep, dis0, dis_GS, myPCSet, dTym, GS_dump, GSstat_dump, mstd_MCS);

        fprintf(err_dump, "%lg %lg", e_GS_ord(0),e_GS_ord(1));
        fprintf(err_dump, "\n");

        // close files
        if(fclose(GS_dump)){
            printf("Could not close file '%s'\n",SoluGSmodes.c_str());
            exit(1);    
        }
        
        if(fclose(GSstat_dump)){
            printf("Could not close file '%s'\n",SoluGSstat.c_str());
            exit(1);    
        }
    }
