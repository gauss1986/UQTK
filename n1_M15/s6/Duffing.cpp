/*===================================================================================== 
DuffingISP.cpp
Modified based on the surf_rxn example in UQTK v.2.1.1
by Lin Gao (lingao.gao@mail.utoronto.ca)
July 25, 2015
===================================================================================== */

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
#include "GhanemSpanos.h"
#include "AAPG.h"
#include "nAAPG.h"
#include "ticktock.h"

#define DIS0 0.0
#define VEL0 0.0

/// \brief c++ version of the Matlab code GS_duffing_multi.m.
/// Dynamic response of a Duffing oscillator subject to uncertain excitation force.
/// The unceratain force is generated using KL expansion.
/// Mass, stiffness, damping and nonlinear coefficients are assumed to be deterministic.

/// Main program of uncertainty propagation of the ODE model excitation force via intrusive spectral projection (ISP)
int main(int argc, char *argv[])
{   int CASE=2;
    int nspl=1000000;       //MCS sample size
    string pcType;  //PC type
    double epsilon = 1;
    double sigma = 2.5;   //Standard deviation
    double dTym = 0.01;
    int ord_AAPG_GS = 2;//GS order in AAPG subproblems
    int dim;        //Stochastic dimensionality
    int nkl;        //Number of terms retained in KL expansion
    char* cov_type; //Covariance type
    bool act_D;     //wheather doing AAPG adaptive or not
    double dof;     //Spatial dof
    //int noutput;    //Number of output points
    Array1D<int> lout(4,0);
    lout(0) = 2/dTym;
    lout(1) = 4/dTym;
    lout(2) = 6.5/dTym;
    lout(3) = 9/dTym;
    unsigned int noutput = lout.XSize();
    double clen;    // Correlation length of the random process
    Array1D<double> inpParams(5,0.e0);//Input parameters
    Array1D<int> init_D(2,0);
    Array1D<int> coeff_D(2,0);
    bool PDF=true;       //control coeff for outputing PDF

    /* Read the user input */
    int c;
    while ((c=getopt(argc,(char **)argv,"C:N:s:e:n:g:P:p"))!=-1){
        switch (c) {
        case 'C':
            CASE = strtod(optarg, (char **)NULL);
            break;
        case 'N':
            nspl= strtod(optarg, (char **)NULL);
            break;
        case 's':
            sigma= strtod(optarg, (char **)NULL);
            break;
        case 'e':
            epsilon= strtod(optarg, (char **)NULL);
            break;
        case 'n':
            dTym = strtod(optarg, (char **)NULL);
            break;
        case 'g':
            ord_AAPG_GS = strtod(optarg, (char **)NULL);
            break;
        case 'P':
            int temp_PDF = strtod(optarg, (char **)NULL);
            if (temp_PDF == 0)
                PDF = false;
            else
                PDF = true;
            break;
        }
    }

    // Time marching info
    double tf = 10;
    // Number of steps
    int nStep=(int) tf / dTym;

    Array1D<double> fbar(2*nStep+1,0.e0);//mean of forcing

    double factor_OD;//dynamical-orthogonal info
    double p;       //Threashold used in the adaptive AAPG
    int ord_GS;     //GS order
    int ord_AAPG;   //AAPG order
    
    if (CASE==1){//Stochastic initial conditions and deterministic forcing
        pcType = "HG";
        dim = 2;
        nkl = 2;
        cov_type = (char *)"Exp";
        //sigma = 0.5;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 2;
        //ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        //noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = epsilon;//epsilon
        double t_temp = 0.0; 
        for (int i=0;i<2*nStep+1;i++){
            fbar(i) = 2.0*(1.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
            t_temp +=dTym/2;
        }
        init_D(0)=0;
        init_D(1)=1;
        coeff_D(0)=0;
        coeff_D(1)=1;
    }
    if (CASE==2){//Stochastic forcing and deterministic initial conditions
        pcType = "LU";
        clen = 0.05;
        dim = 15;
        nkl = 15;
        cov_type = (char *)"Exp";
        //sigma = 2.5;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 3;
        //ord_AAPG_GS = 3;
        act_D = false;
        p = 0.99;
        dof = 2;
        //noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = epsilon;//epsilon
        double t_temp = 0.0; 
        for (int i=0;i<2*nStep+1;i++){
            fbar(i) = 2.0*(2.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
            t_temp +=dTym/2;
        }
        init_D(0)=10000;
        init_D(1)=10001;
        coeff_D(0)=0;
        coeff_D(1)=1;
    }
    if (CASE==3){//stochastic initial conditions and stochastic forcing
        clen = 0.05;
        dim = 12;
        nkl = 10;
        dof = 2;
        cov_type = (char *)"Exp";
        sigma=0.5;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 3;
        //ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        //noutput = 2;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = epsilon;//epsilon
        double t_temp = 0.0; 
        for (int i=0;i<2*nStep+1;i++){
            fbar(i) = 2.0*(1-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
            t_temp +=dTym/2;
        }
        if ((dof+nkl-dim)>10e-5){
            cout << "This test case is configured so that total stochastic dim should equal the number of modes in KL exapansion and dof. Now this is not true!!" << endl<<flush;
            return 1;
        }
        init_D(0)=10;
        init_D(1)=11;
        coeff_D(0)=0;
        coeff_D(1)=1;
    }
    if (CASE==4){//Stochastic zeta and epsilon
        dim = 2;
        nkl = 2;
        cov_type = (char *)"Exp";
        sigma = 0.5;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 2;
        //ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        //noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = epsilon;//epsilon
        inpParams(3) = 0.02;//std for zeta
        inpParams(4) = 0.2;//std for epsilon
        double t_temp = 0.0; 
        for (int i=0;i<2*nStep+1;i++){
            fbar(i) = 2.0*(1.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
            t_temp +=dTym/2;
        }
        init_D(0)=10000;
        init_D(1)=10001;
        coeff_D(0)=0;
        coeff_D(1)=1;
    }
    if (CASE==5){//Stochastic zeta and epsilon and stochastic forcing
        pcType = "LU";
        clen = 0.05;
        dim = 100;
        nkl = 98;
        cov_type = (char *)"Exp";
        sigma = 0.5;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 2;
        //ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        //noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = epsilon;//epsilon
        inpParams(3) = 0.05;//std for zeta
        inpParams(4) = 0.35;//std for epsilon
        double t_temp = 0.0; 
        for (int i=0;i<2*nStep+1;i++){
            fbar(i) = 2.0*(1.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
            t_temp +=dTym/2;
        }
        init_D(0)=10000;
        init_D(1)=10001;
        coeff_D(0)=dim-2;
        coeff_D(1)=dim-1;
    }

    // Save the force
    write_datafile_1d(fbar,"fbar.dat");
  
    
    /* Print the input information on screen */
    cout << " - Number of KL modes:              " << nkl  << endl<<flush;
    cout << " - Monte Carlo sample size:         " << nspl  << endl<<flush;
    if (CASE == 1){
        cout << " - Case 1 Will generate random variable of type "<< pcType << ", with standard deviation " << sigma << endl<< flush;
    }
    if (CASE == 2){
        cout << " - Case 2 with correlation length "<<clen<<", standard deviation "<<sigma<<". Random variable is of type "<<pcType<<"."<<endl<<flush;
    }
    if (CASE == 3){
        cout << " - Case3 where the initial condition and the forcing are both stochastic, with standard deviation " << sigma << ". Correlation length of the process is " << clen << ". Variance type " << pcType  << endl << flush;
    }
    if (CASE == 4){
        cout << " - Case4 where the damping coeff zeta and epsilon are stochastic."<< endl<< flush;
    }
    if (CASE == 5){
        cout << " - Case5 where zeta and epsilon are stochastic." << endl << flush;
    }
    cout << " - Time marching step:              " << dTym  << endl<<flush;
    cout << " - Process end time:                " << tf  << endl<<flush;
    cout << " - Dynamical orthogonal AAPG factor:" << factor_OD  << endl<<flush;
    cout << " - Nonlinearity parameter:          " << inpParams(2)  << endl<<flush;
    cout << " - The threshold in adaptive AAPG:  " << p  << endl<<flush;
    cout << " - Order of GS:                     " << ord_GS  << endl<<flush;
    cout << " - Order of GS in AAPG subproblems: " << ord_AAPG_GS  << endl<<flush;
    cout << " - Order of AAPG:                   " << ord_AAPG  << endl<<flush;
    if (PDF){
        cout << " - Sampling GS/AAPG:true"<<endl<<flush;
    }
    if (act_D){
        cout << " - Adaptive AAPG applied"<<endl<<flush;
    }
    else{
        cout << " - Full AAPG applied"<<endl<<flush;
    }

    // Generate random samples
    Array2D<double> samPts_ori(nspl,dim,0.e0);
    Array2D<double> samPts_norm(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",0,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts_ori);
    if (strcmp(pcType.c_str(),"HG")==0){
        for (int i=0;i<nspl;i++){
            for (int j=0;j<dim;j++){
                // normalize samPts so that it has variance 1/3
                samPts_norm(i,j)=samPts_ori(i,j)*sqrt(1.0/3.0);
            }  
        }
    }
    else if (strcmp(pcType.c_str(),"LU")==0){
        MCPCSet.DrawSampleVar(samPts_norm);
    }

    // Generate the KL expansion of the excitation force, the first nkl terms are stochastic
    Array2D<double> scaledKLmodes(2*nStep+1,dim,0.e0);
    if (CASE==2 || CASE ==3 || CASE==5){
        Array2D<double> scaledKLmodes_small(2*nStep+1,nkl,0.e0);
        genKL(scaledKLmodes_small, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
        Array1D<double> KLmodes_temp(2*nStep+1,0.e0);
        for (int i=0;i<nkl;i++){
            getCol(scaledKLmodes_small,i,KLmodes_temp);
            scaledKLmodes.replaceCol(KLmodes_temp,i);
        }
    }
    write_datafile(scaledKLmodes,"KL.dat");
 
    // Sample initial conditions for MCS 
    Array2D<double> initial_samples(nspl,dim,0.e0);
    Array2D<double> stat_init(dof,2,0.e0);
    stat_init(0,0)=VEL0;
    stat_init(1,0)=DIS0;
    for (int i=0;i<nspl;i++){
        initial_samples(i,0) = VEL0;
        initial_samples(i,1) = DIS0;
    }
    if (CASE==1 || CASE == 3){
        for (int i=0;i<nspl;i++){
            for (int j=0;j<dof;j++){
                if (strcmp(pcType.c_str(),"HG")==0){
                    initial_samples(i,j) += samPts_ori(i,j)*sigma;
                }
                if (strcmp(pcType.c_str(),"LU")==0){
                    initial_samples(i,j) += samPts_ori(i,j)/sqrt(1.0/3.0)*sigma;
                }
            }
        }
        for (int i=0;i<dof;i++){
            Array1D<double> samPts_1D(nspl,0.e0);
            getCol(initial_samples,i,samPts_1D);
            Array1D<double> sample_mstd = mStd(samPts_1D,nspl);
            cout << "Mean of sample on dim " << i << " is "<< sample_mstd(0) << endl;
            cout << "Std of sample on dim "<< i << " is " <<sample_mstd(1) << endl;
            stat_init.replaceRow(sample_mstd,i); 
            ostringstream name;
            name << "sample_dof" << i << ".dat";
            string name_str = name.str();
            write_datafile_1d(samPts_1D,name_str.c_str());
        }
    }

    // Sample epsilon and zeta for MCS
    Array2D<double> temp_inp(nspl,2,0.e0); 
    Array1D<Array1D<double> > inpParams_MCS(nspl);
    for (int i=0;i<nspl;i++){
        temp_inp(i,0) = inpParams(1);
        temp_inp(i,1) = inpParams(2);
        inpParams_MCS(i) = inpParams;
        if (CASE==4 || CASE==5){
            for (int j=0;j<2;j++){
                if (strcmp(pcType.c_str(),"HG")==0){
                    temp_inp(i,j) += samPts_ori(i,j)*inpParams(j+3);
                }
                if (strcmp(pcType.c_str(),"LU")==0){
                    temp_inp(i,j) += samPts_ori(i,j)/sqrt(1.0/3.0)*inpParams(j+3);
                }
                inpParams_MCS(i)(j+1)=temp_inp(i,j);
            }
        }
    }
    Array2D<double> stat_inp(2,2,0.e0);
    for (int i=0;i<2;i++){
        Array1D<double> samPts_1D(nspl,0.e0);
        getCol(temp_inp,i,samPts_1D);
        Array1D<double> sample_mstd = mStd(samPts_1D,nspl);
        cout << "Mean of sample on inpParams " << i << " is "<< sample_mstd(0) << endl;
        cout << "Std of sample on inpParams "<< i << " is " <<sample_mstd(1) << endl;
        //stat_inp.replaceRow(sample_mstd,i); 
        inpParams(i+1)=sample_mstd(0);
        inpParams(i+3)=sample_mstd(1);
        ostringstream name;
        name << "sample_dof" << i << ".dat";
        string name_str = name.str();
        write_datafile_1d(samPts_1D,name_str.c_str());
    }

    // MCS
    // Open files to write out
    string f_MCS_dis = "MCS_dis.dat"; 
    string f_MCS_vel = "MCS_vel.dat"; 
    string force_MCS = "KL_sample.dat"; 
    FILE *f_dump=createfile(f_MCS_dis);
    FILE *f_dump2=createfile(f_MCS_vel);
    FILE *force_dump=createfile(force_MCS);
    // Write solutions at initial step
    Array2D<double>dis_MC(nStep+1,nspl,0.e0);
    Array2D<double>vel_MC(nStep+1,nspl,0.e0);
    Array2D<double>totalforce(2*nStep+1,nspl,0.e0);
    Array1D<Array2D<double> > result(dof);
    result(0) = dis_MC;
    result(1) = vel_MC;
    cout << "\nMCS...\n" << endl; 
    int nthreads;
    // Time marching steps
    TickTock tt;
    tt.tick();
    #pragma omp parallel default(none) shared(result,dof,nStep,nspl,samPts_norm,initial_samples,nkl,dTym,inpParams_MCS,nthreads,fbar,scaledKLmodes,totalforce) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array1D<double> initial(2,0.e0);
        getRow(initial_samples,iq,initial);
        Array1D<double> sampleforce=sample_force(samPts_norm,iq,2*nStep,fbar,nkl,scaledKLmodes,inpParams_MCS(iq));
        for (int i=0;i<2*nStep+1;i++)
            totalforce(i,iq)=sampleforce(i);
        Array2D<double> temp = det(dof, nspl, nStep, dTym, sampleforce, inpParams_MCS(iq), initial);
        for (int it=0;it<nStep+1;it++){
            for (int idof=0;idof<dof;idof++){
                result(idof)(it,iq) = temp(idof,it);
            }
        }
    }
    //output the number of threads
    #pragma omp single
    nthreads = omp_get_num_threads();
    }
    //report the time cost
    tt.tock("Took");
    cout << "Size of result(1) is:" <<result(1).XSize() << "X"<<result(1).YSize()<< endl;
    double t_MCS = tt.silent_tock();
    cout << "Number of threads in OMP:" << nthreads << endl;
    Array2D<double> MCS_s_dis(2,nStep+1,0.e0);
    Array2D<double> MCS_s_vel(2,nStep+1,0.e0);
    Array2D<double> dis_MCS(nspl,noutput+1,0.e0);
    Array2D<double> vel_MCS(nspl,noutput+1,0.e0);
    unsigned int ind = 0;
    Array1D<double> mstd;
    for (int ix=0;ix<nStep+1;ix++){
        Array1D<double> tempforce(nspl,0.e0);
        // save mean/std of force
        getRow(totalforce,2*ix,tempforce);
        mstd  = mStd(tempforce,nspl);
        WriteMeanStdDevToFilePtr(ix*dTym, mstd(0),mstd(1),force_dump);         
        // save mean/std of vel
        Array1D<double> tempvel(nspl,0.e0);
        getRow(result(0),ix,tempvel);
        mstd  = mStd(tempvel,nspl);
        MCS_s_vel.replaceCol(mstd,ix); 
        WriteMeanStdDevToFilePtr(ix*dTym, mstd(0),mstd(1),f_dump2);         
        // save mean/std of dis
        Array1D<double> tempdis(nspl,0.e0);
        getRow(result(1),ix,tempdis);
        mstd  = mStd(tempdis,nspl);
        MCS_s_dis.replaceCol(mstd,ix); 
        WriteMeanStdDevToFilePtr(ix*dTym, mstd(0),mstd(1),f_dump);         
        // report to screen
        //if (ix % ((int) nStep/noutput) == 0){
        if (ind<noutput && ix == lout(ind)){
            WriteMeanStdDevToStdOut(ix, ix*dTym, mstd(0), mstd(1));
            dis_MCS.replaceCol(tempdis,ind); 
            Array1D<double> tempvel(nspl,0.e0);
            getRow(result(0),ix,tempvel);
            vel_MCS.replaceCol(tempvel,ind); 
            ind++;
        }
    }
    write_datafile(dis_MCS,"dis_MCS.dat");
    write_datafile(vel_MCS,"vel_MCS.dat");
    closefile(f_dump, f_MCS_dis); 
    closefile(f_dump2, f_MCS_vel); 
    closefile(force_dump, force_MCS); 
    
    // Ghanem-Spanos method
    printf("\nGhanem-Spanos method...\n");

    Array1D<double> t_GS(ord_GS,0.e0);
    
    ostringstream err_stream;
    err_stream << "error_s" << sigma << "e"<< epsilon<< "dt"<<dTym << ".dat";
    //if (act_D){
        //err_stream << "error_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<"_actD.dat";
    //}
    //else {
    //    err_stream << "error_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<".dat";
    //}
    string errstr(err_stream.str());
    FILE *err_dump;
    if(!(err_dump=fopen(errstr.c_str(),"w"))){
        printf("Could not open file '%s'\n",errstr.c_str());
        exit(1);    
    }
    Array2D<double> e_GS(ord_GS,2);
    Array1D<Array1D<double> > e_GS_sample(ord_GS);
    Array1D<Array1D<double> > e_GS_sample_vel(ord_GS);
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
        f_GS.replaceCol(fbar,0);

        for (int i=0;i<nkl;i++){
            Array1D<double> tempf(2*nStep+1,0.e0);
            getCol(scaledKLmodes,i,tempf);
            f_GS.replaceCol(tempf,i+1);
        }
        ostringstream s_fGS;
        s_fGS << "f_GS" << ord<<".dat";
        string n_fGS(s_fGS.str());
        write_datafile(f_GS,n_fGS.c_str());

       // Initial conditions
        Array1D<Array2D<double> > result(2);
        Array1D<Array1D<double> > initial_GS(2);
        Array1D<double> temp(nPCTerms,0.e0);
        initial_GS(0)=temp;
        initial_GS(1) = temp;
        myPCSet.InitMeanStDv(stat_init(0,0),stat_init(0,1),1,initial_GS(0));
        myPCSet.InitMeanStDv(stat_init(1,0),stat_init(1,1),2,initial_GS(1));

        clock_t start = clock();
    	tt.tick();
        GS(dof, myPCSet, coeff_D, nPCTerms, nStep, initial_GS, dTym, inpParams, f_GS, result);
        t_GS(ord-1) = tt.silent_tock();
	    cout << "Cost time: "<<(clock()-start)/(double)(CLOCKS_PER_SEC)<<endl;
	
        // output and save the solution
        // Open files to write out solutions
        ostringstream s;
        s << "GS_dis_modes_" << ord<<".dat";
        string SoluGSmodes(s.str());
        FILE *GS_dump=createfile(SoluGSmodes);
        ostringstream s0;
        s0 << "GS_vel_modes_" << ord<<".dat";
        string SolvGSmodes(s0.str());
        FILE *GS_dump_v=createfile(SolvGSmodes);
        // Open files to write out statistics of the solutions
        ostringstream s1;
        s1 << "GS_dis_stat_" << ord<<".dat";
        string SoluGSstat(s1.str());
        FILE *GSstat_dump=createfile(SoluGSstat);
        ostringstream s4;
        s4 << "GS_vel_stat_" << ord<<".dat";
        string SolvGSstat(s4.str());
        FILE *GSstat_dump_v=createfile(SolvGSstat);
       
        //assemble dis and output dis_sample
        Array2D<double> stat_GS(2,nStep+1,0.e0);
        Array2D<double> stat_GS_vel(2,nStep+1,0.e0);
        Array2D<double> et_dis(nStep+1,2,0.e0);
        Array2D<double> et_vel(nStep+1,2,0.e0);
        Array1D<double> e_GS_ord = postprocess_GS(lout, nPCTerms, nStep, result(1), myPCSet, dTym,GS_dump, GSstat_dump, MCS_s_dis, stat_GS, et_dis);
        Array1D<double> e_GS_ord_vel = postprocess_GS(lout, nPCTerms, nStep, result(0), myPCSet, dTym, GS_dump_v, GSstat_dump_v, MCS_s_vel, stat_GS_vel, et_vel);
        Array1D<double> normsq(nPCTerms,0.e0); 
        ostringstream s5;
        s5 << "et_dis_GS" << ord<<".dat";
        string et(s5.str());
        write_datafile(et_dis,et.c_str()); 
        //myPCSet.OutputNormSquare(normsq);
        //cout << "Normsq 1=" <<normsq(1)<<endl;
        //cout << "Normsq 2=" <<normsq(2)<<endl;

        if (PDF){ 
        cout << "Sampling dis..."<< endl;
        Array2D<double> GS_dis_sampt=sampleGS(lout,dim, nStep, nPCTerms, myPCSet, result(1), samPts_norm, stat_GS, e_GS_sample(ord-1));
        ostringstream s2;
        s2 << "GS_dis_sample" << ord<<".dat";
        string SoluGSsample(s2.str());
        write_datafile(GS_dis_sampt,SoluGSsample.c_str());
        //ouput vel_sample
        cout << "Sampling vel..."<< endl;
        Array2D<double> GS_vel_sampt=sampleGS(lout, dim, nStep, nPCTerms, myPCSet, result(0), samPts_norm, stat_GS_vel, e_GS_sample_vel(ord-1));
        ostringstream s3;
        s3 << "GS_vel_sample" << ord<<".dat";
        SoluGSsample=s3.str();
        write_datafile(GS_vel_sampt,SoluGSsample.c_str());
        }

        fprintf(err_dump, "%lg %lg", e_GS_ord(0),e_GS_ord(1));
        fprintf(err_dump, "\n");
        e_GS.replaceRow(e_GS_ord,ord-1);

        closefile(GS_dump, SoluGSmodes); 
        closefile(GS_dump_v, SolvGSmodes); 
        closefile(GSstat_dump, SoluGSstat); 
        closefile(GSstat_dump_v, SolvGSstat); 
    }

    // AAPG
    printf("AAPG...\n");

   
    Array2D<double> e_AAPG(ord_AAPG,2,0.e0); 
    Array1D<Array1D<double> > e_sample_AAPG_dis(ord_AAPG); 
    Array1D<Array1D<double> > e_sample_AAPG_vel(ord_AAPG); 

    Array1D<double> t_AAPG = AAPG(dof, inpParams, fbar, dTym, ord_AAPG_GS, pcType, lout, dim, nStep, scaledKLmodes, factor_OD, ord_AAPG, act_D, p, MCS_s_dis, err_dump, stat_init, samPts_norm, e_AAPG, e_sample_AAPG_dis, e_sample_AAPG_vel, init_D, coeff_D, PDF);
    
    // output the timing
    Array1D<double> t(1+ord_GS+ord_AAPG,0.e0);
    t(0) = t_MCS;
    for (int i=0;i<ord_GS;i++)
        t(i+1)=t_GS(i);
    t(ord_GS+1)=t_AAPG(0)+t_AAPG(1);
    t(ord_GS+2)=t_AAPG(0)+t_AAPG(1)+t_AAPG(2);
    if (ord_AAPG>=3)
        t(ord_GS+3)=t_AAPG(0)+t_AAPG(1)+t_AAPG(2)+t_AAPG(3);
    // output the error info
    printf("\n");
    cout << "Printing the error of GS in displacement..." << endl;
    for(int i=0;i<ord_GS;i++){
        cout << "GS ord " << i+1<< " is em="<<e_GS(i,0) <<", es= "<< e_GS(i,1)<< endl;
    }
    cout << "Printing the error of AAPG in displacement..." << endl;
    for(int i=0;i<ord_AAPG;i++){
        cout << "AAPG ord " << i+1<< " is em="<<e_AAPG(i,0) <<", es= "<< e_AAPG(i,1)<<endl;
    }
    if (PDF){
    cout << "Printing the error of MCS assembled GS in displacement..." << endl;
    for(int i=0;i<ord_GS;i++){
        cout << "GS ord " << i+1<< " is em="<<e_GS_sample(i)(0) <<", es= "<< e_GS_sample(i)(1)<< endl;
    }
    cout << "Printing the error of MCS assembled AAPG in displacement..." << endl;
    for(int i=0;i<ord_AAPG;i++){
        cout << "AAGP ord " << i+1<< " is em="<<e_sample_AAPG_dis(i)(0) <<", es= "<< e_sample_AAPG_dis(i)(1)<< endl;
    }
    }

    ostringstream time_stream;
    if (act_D){
        time_stream << "time_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<"_actD.dat";
    }
    else{
        time_stream << "time_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<".dat";
    }
    string timestr(time_stream.str());
    WriteToFile(t, timestr.c_str());
        
    if(fclose(err_dump)){
        printf("Could not close file '%s'\n",errstr.c_str());
        exit(1);    
    }

    return 0;
}

