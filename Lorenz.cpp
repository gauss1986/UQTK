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
#include <omp.h>
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
#define SIG 0.01*sqrt(3)
#define ORDER_GS 2
#define ORDER_AAPG_GS 2
#define ORDER_AAPG 3
#define TFINAL 10.0
#define DTYM 0.01
#define NSPL 1000000
#define AMP 0.0
#define FBAR 8.0
#define x0 0.0 // A point on or nearly on the attractor by Lorenz 2005
#define y0 0.0
#define z0 0.0
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
    //cout << " - Will generate covariance of type "<<cov_type<<", with correlation length "<<clen<<" and standard deviation "<<sigma<<endl<<flush;
    cout << " - Standard deviation:              " << sigma <<endl <<flush;
    cout << " - Initial condition:               " << "x0="<<x0 <<",y0="<< y0<<",z0="<<z0 <<endl <<flush;
    cout << " - Mean of force:                   " << FBAR <<endl <<flush;
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
    //Array2D<double> scaledKLmodes(2*nStep+1,nkl,0.0);
    //genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
    //write_datafile(scaledKLmodes,"Lorenz/KL.dat");
 
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
    int Rand_force = 1;
    // Draw Monte Carlo samples
    Array2D<double> samPts(nspl,dim,0.e0);
    PCSet MCPCSet("NISP",ord_GS,dim,pcType,0.0,sigma);//alpha and beta coeff is useless in LU and HG
    MCPCSet.DrawSampleVar(samPts);
    for (int i=0;i<nspl+1;i++){
        samPts(i,0)=samPts(i,0)*sigma+x0;
        samPts(i,1)=samPts(i,1)*sigma+y0;
        samPts(i,2)=samPts(i,2)*sigma+z0;
    }
    Array2D<double> sample_mstd_2D(dof,2);
    for (int i=0;i<dof;i++){
        Array1D<double> samPts_1D(nspl,0.e0);
        getCol(samPts,i,samPts_1D);
        Array1D<double> sample_mstd = mStd(samPts_1D,nspl);
        cout << "Mean of sample on dof " << i << " is "<< sample_mstd(0) << endl;
        cout << "Std of sample on dof "<< i << " is " <<sample_mstd(1) << endl;
        sample_mstd_2D.replaceRow(sample_mstd,i); 
        ostringstream name;
        name << "sample_dof" << i << ".dat";
        string name_str = name.str();
        write_datafile_1d(samPts_1D,name_str.c_str());
    }
    // Open files to write out
    string Solufile_mean = "Lorenz/MCS_mean.dat";
    string Solufile_std = "Lorenz/MCS_std.dat";
    FILE *f_mean = createfile(Solufile_mean);
    FILE *f_std  = createfile(Solufile_std);

    // Initialize solutions
    Array2D<double> x(nStep+1,nspl,0.e0);
    //Array2D<double> y(nStep+1,nspl,0.e0);
    //Array2D<double> z(nStep+1,nspl,0.e0);
    //Array1D<Array2D<double> > result(1);
    //result(0) = x;
    //result(1) = y;
    //result(2) = z;
    // Define and write solutions at initial step
    //WriteMeanStdDevToFilePtr_lorenz(0, x0, y0, z0, f_mean);        
    //WriteMeanStdDevToFilePtr_lorenz(0, sigma, sigma, sigma, f_std);        
    Array1D<double> totalforce(2*nStep,fbar);
    cout << "\nMCS...\n" << endl; 
    int nthreads;
    // Time marching steps
    TickTock tt;
    tt.tick();
    #pragma omp parallel default(none) shared(x,dof,nStep,nspl,samPts,nkl,dTym,inpParams,nthreads,totalforce) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl+1;iq++){
        // sample initial condition
        Array1D<double> initial(dof,0.e0);
        getRow(samPts,iq,initial);
        // compute the solution with different initial conditions
        Array2D<double> temp = det(dof, nspl, nStep, nkl, dTym, totalforce, inpParams, initial);
        //for (int idof=0;idof<dof;idof++){
            for (int it=0;it<nStep+1;it++){
                x(it,iq) = temp(0,it);
            }
        //}
    }
    //output the number of threads
    #pragma omp single
    nthreads = omp_get_num_threads();
    }
    //report the time cost
    tt.tock("Took");
    double t_MCS = tt.silent_tock();
    cout << "Number of threads in OMP:" << nthreads << endl;
    Array2D<double> mean(dof,nStep+1,0.e0);
    Array2D<double> std(dof,nStep+1,0.e0);
    Array1D<Array2D<double> > mstd_MCS(dof);
    for (int ivar=0;ivar<dof;ivar++){
        Array2D<double> temp2(2,nStep+1,0.e0);
        mstd_MCS(ivar) = temp2;
    }
    // post-process the solution
    #pragma omp parallel default(none) shared(x,nStep,dof,nspl,mstd_MCS,mean,std) 
    {
    #pragma omp for
    for (int ix=0;ix<nStep+1;ix++){
        //for (int ivar=0;ivar<1;ivar++){
            Array1D<double> temp(nspl,0.e0);
            getRow(x,ix,temp);
            Array1D<double> mstd(2,0.e0);
            mstd = mStd(temp,nspl);
            mstd_MCS(0)(0,ix)=mstd(0); 
            mstd_MCS(0)(1,ix)=mstd(1); 
            mean(0,ix) = mstd(0);
            std(0,ix) = mstd(1);
        //}
    }
    }
    // save results
    for (int ix=0;ix<nStep+1;ix++){
        WriteMeanStdDevToFilePtr_lorenz(ix*dTym, mean(0,ix),mean(1,ix),mean(2,ix),f_mean);         
        WriteMeanStdDevToFilePtr_lorenz(ix*dTym, std(0,ix),std(1,ix),std(2,ix),f_std);         
    }
    cout << "Results of mean..." << endl;
    // Output mean to scren
    for (int ix=0;ix<nStep+1;ix++){
        if (ix % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut_lorenz(ix, ix*dTym, mean(0,ix), mean(1,ix), mean(2,ix));
        }
    }
    cout << "Results of std..." << endl;
    // Output std to scren
    for (int ix=0;ix<nStep+1;ix++){
        if (ix % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut_lorenz(ix, ix*dTym, std(0,ix), std(1,ix), std(2,ix));
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
	    // The number of PC terms
        const int nPCTerms = myPCSet.GetNumberPCTerms();
        Array1D<double> normsq(nPCTerms,0.e0); 
        myPCSet.EvalNormSq(normsq);
        tt.tock("Took");
	    cout << "Order "<< ord << endl;

        cout << "The number of PC terms in an expansion is " << nPCTerms << endl;
        // Print the multiindices on screen
        for (int i=0;i<nPCTerms;i++)
            cout << "Normsq of term No." << i << " is " << normsq(i) << endl;
    
        // Prepare the force in PC format
        Array2D<double> f_GS(2*nStep+1,nPCTerms,0.e0);
        for (int i=0;i<2*nStep+1;i++){
            f_GS(i,0) = fbar;
        }
        //for (int i=0;i<nkl;i++){
        //    Array1D<double> tempf(2*nStep+1,0.e0);
        //    getCol(scaledKLmodes,i,tempf);
        //    f_GS.replaceCol(tempf,i+1);
        //}

        // Allocate result variable
        Array1D<Array2D<double> > result(3);
        // initial conditions
        Array1D<Array1D<double> > initial_GS(dof);
        Array1D<double> temp(nPCTerms,0.e0);
        for (int i=0;i<dof;i++){
            initial_GS(i)=temp;
            myPCSet.InitMeanStDv(sample_mstd_2D(i,0),sample_mstd_2D(i,1),i+1,initial_GS(i));
        }

        clock_t start = clock();
    	tt.tick();
        GS(dof, myPCSet, ord, dim, nPCTerms, pcType, nStep, initial_GS, dTym, inpParams, f_GS, result);
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
        
        Array1D<double> e_GS_ord_x = postprocess_GS(nPCTerms, nStep, result(0), myPCSet, dTym, GS_dump_x, GSstat_dump_x, mstd_MCS(0));
        Array1D<double> e_GS_ord_y = postprocess_GS(nPCTerms, nStep, result(1), myPCSet, dTym, GS_dump_y, GSstat_dump_y, mstd_MCS(1));
        Array1D<double> e_GS_ord_z = postprocess_GS(nPCTerms, nStep, result(2), myPCSet, dTym, GS_dump_z, GSstat_dump_z, mstd_MCS(2));

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
