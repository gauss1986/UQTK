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
#include "uqtktools.h"
#include "uqtkmcmc.h"
#include "PCSet.h"
#include "arraytools.h"
#include "getopt.h"

#include "UtilsDuffing.h"
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
  printf(" -a <AMP>           : define the oscillation amp of F (default=%e) \n",A);
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
    double fbar = FBAR;
    double Amp = AMP;
    Array2D<double> dis_MC(nStep+1,nspl,0.e0);
    Array2D<double> samPts(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",ord_GS,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts);
    // Open files to write out
    string Solufile = "MCS.dat";
    FILE *f_dump;
    if(!(f_dump=fopen(Solufile.c_str(),"w"))){
        printf("Could not open file '%s'\n",Solufile.c_str());
        exit(1);    
    }
    // Write solutions at initial step
    WriteMeanStdDevToFilePtr(0, 0, 0, f_dump);         
    cout << "\nMCS...\n" << endl;  
    double t_MCS = MCS(nspl, dim, nStep, nkl, dTym, fbar, scaledKLmodes, samPts, dis_MC);
    Array2D<double> mstd_MCS(2,nStep,0.e0);
    // post-process the solution
    for (int ix=0;ix<nStep;ix++){
        Array1D<double> tempdis(nspl,0.e0);
        getRow(dis_MC,ix,tempdis);
        Array1D<double> mstd(2,0.e0);
        mstd = mStd(tempdis,nspl);
        mstd_MCS.replaceCol(mstd,ix); 
        WriteMeanStdDevToFilePtr((ix+1)*dTym, mstd(0),mstd(1),f_dump);         
        if ((ix+1) % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut(ix+1, (ix+1)*dTym, mstd(0), mstd(1));
        }
    }
    if(fclose(f_dump)){
        printf("Could not close file '%s'\n",Solufile.c_str());
        exit(1);    
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
        Array2D<double> f_GS(nStep+1,nPCTerms,0.e0);
        for (int i=0;i<nStep+1;i++){
            f_GS(i,0) = fbar;
        }
        for (int i=0;i<nkl;i++){
            Array1D<double> tempf(nStep+1,0.e0);
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

    // AAPG
    printf("\nAAPG...\n");
    TickTock tt;
    tt.tick();
    PCSet myPCSet("ISP",ord_AAPG_GS,dim,pcType,0.0,1.0,false); 
    tt.tock("Took");
    
    Array1D<double> t_AAPG = AAPG(inpParams, fbar, dTym, ord_AAPG_GS, pcType, dim, nStep, scaledKLmodes, dis0, vel0, myPCSet, factor_OD, ord_AAPG, act_D, p, mstd_MCS, err_dump);
    
    // output the timing
    Array1D<double> t(3+ord_GS+ord_AAPG,0.e0);
    t(0) = t_MCS;
    for (int i=0;i<ord_GS;i++)
	    t(i+1)=t_GS(i);
    for (int i=0;i<ord_AAPG+2;i++)
        t(i+1+ord_GS)=t_AAPG(i);

    ostringstream time_stream;
    if (act_D){
        time_stream << "time_n" << dim << "_e"<<epsilon<<"_s"<<sigma<<"_actD.dat";
    }
    else{
        time_stream << "time_n" << dim << "_e"<<epsilon<<"_s"<<sigma<<".dat";
    }
    string timestr(time_stream.str());
    WriteToFile(t, timestr.c_str());
        
    if(fclose(err_dump)){
        printf("Could not close file '%s'\n",errstr.c_str());
        exit(1);    
    }

    return 0;
}
