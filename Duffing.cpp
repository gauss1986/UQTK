/*===================================================================================== 
DuffingISP.cpp
Modified based on the surf_rxn example in UQTK v.2.1.1
by Lin Gao (lingao.gao@mail.utoronto.ca)
July 25, 2015
===================================================================================== */

#include <math.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "uqtktools.h"
#include "uqtkmcmc.h"
#include "PCSet.h"
#include "arraytools.h"
#include "getopt.h"

#include "Utils.h"
#include "Duffing.h"
#include "KL.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "AAPG.h"
#include "ticktock.h"

#define DIM 2
#define CLEN 0.1
#define SIG 1.0
#define ORDER_GS 2
#define ORDER_AAPG_GS 2
#define ORDER_AAPG 2
#define TFINAL 10.0
#define DTYM 0.01
#define NSPL 100000
#define ZETA 0.1
#define EPSILON 1.0
#define FBAR 2.0
#define DIS0 0.0
#define VEL0 0.0
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
  printf(" -f <factor_OD>   : define the factor in overdetermined dynimcal-orthogonal calculation (default=%e) \n", FACTOR_OD);
  printf(" -e <epsilon>     : define the nonlinearity parameter (default=%e) \n", EPSILON);
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
    // Nonlinearity parameter
    double epsilon = EPSILON;
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
    // Spatial dof
    double dof = 2;
    // Define input parameters
    Array1D<double> inpParams(3,0.e0);
    inpParams(0) = 3;  // code for problem: 0-Duffing, 1-Lorenz, 3-Duffing with stochastic initial condition
    inpParams(1) = ZETA;
    inpParams(2) = epsilon;

    // Time marching info
    double dTym = DTYM;
    double tf = TFINAL;
    // Number of steps
    int nStep=(int) tf / dTym;
    // mean of forcing
    double fbar = FBAR;
  
    /* Read the user input */
    int c;

    while ((c=getopt(argc,(char **)argv,"h:c:d:n:p:s:l:t:f:e:G:A:P:D:T:"))!=-1){
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
        case 'e':
            epsilon = strtod(optarg, (char **)NULL);
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
    cout << " - Number of KL modes:              " << dim  << endl<<flush;
    cout << " - Monte Carlo sample size:         " << nspl  << endl<<flush;
    if (abs(inpParams(0)-1)<1e-10){
        cout << " - Will generate covariance of type "<<cov_type<<", with correlation length "<<clen<<" and standard deviation "<<sigma<<endl<<flush;
    }
    if (abs(inpParams(0)-3)<1e-10){
        cout << "- Will generate random variable of type "<< pcType << ", with standard deviation " << sigma << endl<< flush;
    }
    cout << " - Time marching step:              " << dTym  << endl<<flush;
    cout << " - Process end time:                " << tf  << endl<<flush;
    cout << " - Dynamical orthogonal AAPG factor:" << factor_OD  << endl<<flush;
    cout << " - Nonlinearity parameter:          " << epsilon  << endl<<flush;
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

    // Generate the KL expansion of the excitation force
    int nkl = dim;
    Array2D<double> scaledKLmodes(2*nStep+1,nkl,0.0);
    if (abs(inpParams(0)-1)<1e-10){
        genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
    }
    write_datafile(scaledKLmodes,"KL.dat");
 
    // Monte Carlo simulation
    int Rand_force = 1;
    Array2D<double> samPts(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",ord_GS,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts);
    for (int i=0;i<nspl;i++){
        samPts(i,0)=samPts(i,0)*sigma+VEL0;
        samPts(i,1)=samPts(i,1)*sigma+DIS0;
    }
    // Examine the mean/std of the sample
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
    string Solufile = "MCS.dat";
    FILE *f_dump;
    if(!(f_dump=fopen(Solufile.c_str(),"w"))){
        printf("Could not open file '%s'\n",Solufile.c_str());
        exit(1);    
    }
    // Write solutions at initial step
    Array2D<double>dis_MC(nStep+1,nspl,0.e0);
    Array2D<double>vel_MC(nStep+1,nspl,0.e0);
    Array1D<Array2D<double> > result(dof);
    result(0) = dis_MC;
    result(1) = vel_MC;
    cout << "\nMCS...\n" << endl; 
    int nthreads;
    Array1D<double> temp_init(nspl,0);
    Array1D<double> temp_init2(nspl,0);
    // Time marching steps
    TickTock tt;
    tt.tick();
    #pragma omp parallel default(none) shared(temp_init,temp_init2,result,dof,nStep,nspl,samPts,nkl,dTym,inpParams,nthreads,fbar,scaledKLmodes) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array1D<double> initial(2,0.e0);
        if (abs(inpParams(0)-1)<1e-10){
            initial(0) = VEL0;
            initial(1) = DIS0;
        }
        if (abs(inpParams(0)-3)<1e-10){
            initial(0) = samPts(iq,0);
            initial(1) = samPts(iq,1);
        }
        Array1D<double> totalforce=sample_force(samPts,iq,2*nStep,fbar,nkl,scaledKLmodes,inpParams);
        Array2D<double> temp = det(dof, nspl, nStep, nkl, dTym, totalforce, inpParams, initial);
        for (int it=0;it<nStep+1;it++){
            for (int idof=0;idof<dof;idof++){
                result(idof)(it,iq) = temp(idof,it);
            }
        }
        temp_init2(iq)=initial(1);
        temp_init(iq)=temp(1,0);
        //test if the initial condition is the same as passed in
        //if (abs(initial(1)-temp(1,0))>10e-5 || abs(initial(0)-temp(0,0))>10e-5){
        //    temp_init(iq)=1;
        //}
        
    }
    //output the number of threads
    #pragma omp single
    nthreads = omp_get_num_threads();
    }
    write_datafile_1d(temp_init,"initial.dat");
    write_datafile_1d(temp_init2,"initial2.dat");
    // output the test result
    //for (int iq=0;iq<nspl;iq++){
    //    if (temp_init(iq)==1)
    //        cout << "The initial state of the solution is different from the initial condition passd in on dof No." << iq << endl;
    //}
    //report the time cost
    tt.tock("Took");
    cout << "Size of result(1) is:" <<result(1).XSize() << "X"<<result(1).YSize()<< endl;
    double t_MCS = tt.silent_tock();
    cout << "Number of threads in OMP:" << nthreads << endl;
    Array2D<double> mstd_MCS(2,nStep+1,0.e0);
    //Array1D<double> tempdis(nspl,0.e0);
    //getRow(result(1),0,tempdis);
    //write_datafile_1d(tempdis,"MCS_dis.dat");
    // post-process the solution
    for (int ix=0;ix<nStep+1;ix++){
        Array1D<double> tempdis(nspl,0.e0);
        getRow(result(1),ix,tempdis);
        Array1D<double> mstd  = mStd(tempdis,nspl);
        mstd_MCS.replaceCol(mstd,ix); 
        WriteMeanStdDevToFilePtr(ix*dTym, mstd(0),mstd(1),f_dump);         
        if (ix % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut(ix, ix*dTym, mstd(0), mstd(1));
        }
    }
    if(fclose(f_dump)){
        printf("Could not close file '%s'\n",Solufile.c_str());
        exit(1);    
    }
    
    // Ghanem-Spanos method
    printf("\nGhanem-Spanos method...\n");

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
            f_GS(i,0) = fbar;
        }
        for (int i=0;i<nkl;i++){
            Array1D<double> tempf(2*nStep+1,0.e0);
            getCol(scaledKLmodes,i,tempf);
            f_GS.replaceCol(tempf,i+1);
        }

        // Initial conditions
        Array1D<Array2D<double> > result(2);
        Array1D<Array1D<double> > initial_GS(2);
        Array1D<double> temp(nPCTerms,0.e0);
        initial_GS(0)=temp;
        initial_GS(1) = temp;
        if (abs(inpParams(0)-1)<1e-10){
            initial_GS(0)(0) = VEL0;
            initial_GS(1)(0)=DIS0;
        }
        if (abs(inpParams(0)-3)<1e-10){
            myPCSet.InitMeanStDv(sample_mstd_2D(0,0),sample_mstd_2D(0,1),1,initial_GS(0));
            myPCSet.InitMeanStDv(sample_mstd_2D(1,0),sample_mstd_2D(1,1),2,initial_GS(1));
        }

        clock_t start = clock();
    	tt.tick();
        GS(dof, myPCSet, ord, dim, nPCTerms, pcType, nStep, initial_GS, dTym, inpParams, f_GS, result);
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
        
        Array1D<double> e_GS_ord = postprocess_GS(nPCTerms, nStep, result(1), myPCSet, dTym, GS_dump, GSstat_dump, mstd_MCS);

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
    tt.tick();
    PCSet myPCSet("ISP",ord_AAPG_GS,dim,pcType,0.0,1.0,false); 
    tt.tock("Took");
    
    Array1D<double> t_AAPG = AAPG(dof, inpParams, fbar, dTym, ord_AAPG_GS, pcType, dim, nStep, scaledKLmodes, myPCSet, factor_OD, ord_AAPG, act_D, p, mstd_MCS, err_dump, sample_mstd_2D);
    
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

