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
#include "ticktock.h"

#define CASE 1

#define DIS0 0.0
#define VEL0 0.0

/// \brief c++ version of the Matlab code GS_duffing_multi.m.
/// Dynamic response of a Duffing oscillator subject to uncertain excitation force.
/// The unceratain force is generated using KL expansion.
/// Mass, stiffness, damping and nonlinear coefficients are assumed to be deterministic.

/// Main program of uncertainty propagation of the ODE model excitation force via intrusive spectral projection (ISP)
int main(int argc, char *argv[])
{
    string pcType;  //PC type
    int dim;        //Stochastic dimensionality
    char* cov_type; //Covariance type
    double sigma;   //Standard deviation
    int nspl;       //MCS sample size
    double factor_OD;//dynamical-orthogonal info
    int ord_GS;     //GS order
    int ord_AAPG;   //AAPG order
    int ord_AAPG_GS;//GS order in AAPG subproblems
    bool act_D;     //wheather doing AAPG adaptive or not
    double p;       //Threashold used in the adaptive AAPG
    double dof;     //Spatial dof
    int noutput;    //Number of output points
    Array1D<double> inpParams(3,0.e0);//Input parameters
    double clen;
    double FBAR;    //Scalar defining mean of f
    
    if (CASE==1){//Stochastic forcing and deterministic initial conditions
        // Correlation length
        clen = 0.05;
        pcType = "HG";
        dim = 10;
        cov_type = (char *)"Exp";
        sigma = 0.8;
        nspl = 100000;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 3;
        ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = 1.0;//epsilon
        FBAR = 2.0;
    }
    if (CASE==2){//Stochastic initial conditions and deterministic forcing
        pcType = "HG";
        dim = 2;
        cov_type = (char *)"Exp";
        sigma = 0.5;
        nspl = 100000;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 2;
        ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = 1.0;//epsilon
        FBAR = 2.0;
    }
    if (CASE==3){//stochastic initial conditions and stochastic forcing
        clen = 0.05;
        pcType = "LU";
        dim = 10;
        cov_type = (char *)"Exp";
        sigma=0.5;
        nspl = 100000;
        factor_OD = 1.0;
        ord_GS = 2;
        ord_AAPG = 3;
        ord_AAPG_GS = 2;
        act_D = false;
        p = 0.99;
        dof = 2;
        noutput = 10;
        inpParams(0) = 0.0;//Problem to solve 
        inpParams(1) = 0.1;//zeta
        inpParams(2) = 1.0;//epsilon
        FBAR = 2.0;
    }

    // Time marching info
    double dTym = 0.01;
    double tf = 10;
    // Number of steps
    int nStep=(int) tf / dTym;
    // mean of forcing
    Array1D<double> fbar(2*nStep+1,0.e0);
    double t_temp = 0.0; 
    for (int i=0;i<2*nStep+1;i++){
        fbar(i) = FBAR*(1-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
        t_temp +=dTym/2;
    }
    write_datafile_1d(fbar,"fbar.dat");
  
    /* Read the user input */
    int c;

    while ((c=getopt(argc,(char **)argv,"h:c:d:n:p:s:l:t:f:e:G:A:P:D:T:"))!=-1){
        switch (c) {
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
            inpParams(2) = strtod(optarg, (char **)NULL);
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
    if (CASE == 1){
        cout << " - Will generate KL expansion with covariance of type "<<cov_type<<", with correlation length "<<clen<<" and standard deviation "<<sigma<<endl<<flush;
    }
    if (CASE == 2){
        cout << " - Will generate random variable of type "<< pcType << ", with standard deviation " << sigma << endl<< flush;
    }
    if (CASE == 3){
        cout << " - Will generate ranodm variable of type "<< pcType <<", with standard deviation and  KL expansion with covariance of type "<<cov_type<<", with correlation length "<<clen<<" and standard deviation "<<sigma<<endl<<flush; 
    }
    cout << " - Time marching step:              " << dTym  << endl<<flush;
    cout << " - Process end time:                " << tf  << endl<<flush;
    cout << " - Dynamical orthogonal AAPG factor:" << factor_OD  << endl<<flush;
    cout << " - Nonlinearity parameter:          " << inpParams(2)  << endl<<flush;
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
    if (CASE==1 || CASE ==3){
        genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
    }
    write_datafile(scaledKLmodes,"KL.dat");
 
    // Generate random samples
    Array2D<double> samPts(nspl,dim,0.e0);
    Array2D<double> samPts_ori(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",ord_GS,dim,pcType,0.0,1.0);
    Array2D<double> sample_mstd_2D(dof,2,0.e0);
    sample_mstd_2D(0,0)=VEL0;
    sample_mstd_2D(1,0)=DIS0;
    MCPCSet.DrawSampleVar(samPts);
    if (CASE==1){
        if (strcmp(pcType.c_str(),"HG")==0){
            for (int i=0;i<nspl;i++){
                for (int j=0;j<dim;j++){
                    samPts(i,j)=samPts(i,j)*sqrt(0.3333);
                }
            }
        }
    }
    if (CASE==2){
        MCPCSet.DrawSampleVar(samPts);
        for (int i=0;i<nspl;i++){
            if (strcmp(pcType.c_str(),"HG")==0){
                samPts(i,0)=samPts(i,0)*sigma+VEL0;
                samPts(i,1)=samPts(i,1)*sigma+DIS0;
            }
            if (strcmp(pcType.c_str(),"LU")==0){
                samPts(i,0)=samPts(i,0)/sqrt(1.0/3.0)*sigma+VEL0;
                samPts(i,1)=samPts(i,1)/sqrt(1.0/3.0)*sigma+DIS0;
            }
        }
    }
    // Examine the mean/std of the sample
    for (int i=0;i<dim;i++){
        Array1D<double> samPts_1D(nspl,0.e0);
        getCol(samPts,i,samPts_1D);
        Array1D<double> sample_mstd = mStd(samPts_1D,nspl);
        cout << "Mean of sample on dim " << i << " is "<< sample_mstd(0) << endl;
        cout << "Std of sample on dim "<< i << " is " <<sample_mstd(1) << endl;
        if (CASE==2 || CASE==3){
            sample_mstd_2D.replaceRow(sample_mstd,i); 
        }
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
    #pragma omp parallel default(none) shared(result,dof,nStep,nspl,samPts,nkl,dTym,inpParams,nthreads,fbar,scaledKLmodes,totalforce) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array1D<double> initial(2,0.e0);
        if (CASE == 1){
            initial(0) = VEL0;
            initial(1) = DIS0;
        }
        if (CASE == 2){
            initial(0) = samPts(iq,0);
            initial(1) = samPts(iq,1);
        }
        Array1D<double> sampleforce=sample_force(samPts,iq,2*nStep,fbar,nkl,scaledKLmodes,inpParams);
        for (int i=0;i<2*nStep+1;i++)
            totalforce(i,iq)=sampleforce(i);
        Array2D<double> temp = det(dof, nspl, nStep, nkl, dTym, sampleforce, inpParams, initial);
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
    int ind = 0;
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
        if (ix % ((int) nStep/noutput) == 0){
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
    if (act_D){
        err_stream << "error_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<"_actD.dat";
    }
    else {
        err_stream << "error_n" << dim << "_e"<<inpParams(2)<<"_s"<<sigma<<".dat";
    }
    string errstr(err_stream.str());
    FILE *err_dump;
    if(!(err_dump=fopen(errstr.c_str(),"w"))){
        printf("Could not open file '%s'\n",errstr.c_str());
        exit(1);    
    }
    Array2D<double> samPts_norm(nspl,2,0.e0);
    Array2D<double> e_GS(ord_GS,2);
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

       // Initial conditions
        Array1D<Array2D<double> > result(2);
        Array1D<Array1D<double> > initial_GS(2);
        Array1D<double> temp(nPCTerms,0.e0);
        initial_GS(0)=temp;
        initial_GS(1) = temp;
        myPCSet.InitMeanStDv(sample_mstd_2D(0,0),sample_mstd_2D(0,1),1,initial_GS(0));
        myPCSet.InitMeanStDv(sample_mstd_2D(1,0),sample_mstd_2D(1,1),2,initial_GS(1));

        clock_t start = clock();
    	tt.tick();
        GS(dof, myPCSet, ord, dim, nPCTerms, nStep, initial_GS, dTym, inpParams, f_GS, result);
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
        Array1D<double> e_GS_ord = postprocess_GS(noutput, nPCTerms, nStep, result(1), myPCSet, dTym,GS_dump, GSstat_dump, MCS_s_dis);
        Array1D<double> e_GS_ord_vel = postprocess_GS(noutput, nPCTerms, nStep, result(0), myPCSet, dTym, GS_dump_v, GSstat_dump_v, MCS_s_vel);
        Array1D<double> normsq(nPCTerms,0.e0); 
        myPCSet.OutputNormSquare(normsq);
        cout << "Normsq 1=" <<normsq(1)<<endl;
        cout << "Normsq 2=" <<normsq(2)<<endl;

        if (ord == 1){
            for (int i=0;i<nspl;i++){
                if (strcmp(pcType.c_str(),"HG")==0){
                    samPts_norm(i,0)=samPts_ori(i,0)*sqrt(normsq(1));
                    samPts_norm(i,1)=samPts_ori(i,1)*sqrt(normsq(2));
                }
                else if (strcmp(pcType.c_str(),"LU")==0){
                    samPts_norm(i,0)=samPts_ori(i,0);
                    samPts_norm(i,1)=samPts_ori(i,1);
                }
                if (CASE==2){
                    samPts_norm(i,0)=samPts_norm(i,0)+VEL0;
                    samPts_norm(i,1)=samPts_norm(i,1)+DIS0;
                }
            }
        }
        cout << "Sampling dis..."<< endl;
        Array2D<double> GS_dis_sampt=sampleGS(noutput,dim, nStep, nPCTerms, myPCSet, result(1), samPts_norm);
        ostringstream s2;
        s2 << "GS_dis_sample" << ord<<".dat";
        string SoluGSsample(s2.str());
        write_datafile(GS_dis_sampt,SoluGSsample.c_str());
        //ouput vel_sample
        cout << "Sampling vel..."<< endl;
        Array2D<double> GS_vel_sampt=sampleGS(noutput, dim, nStep, nPCTerms, myPCSet, result(0), samPts_norm);
        ostringstream s3;
        s3 << "GS_vel_sample" << ord<<".dat";
        SoluGSsample=s3.str();
        write_datafile(GS_vel_sampt,SoluGSsample.c_str());

        fprintf(err_dump, "%lg %lg", e_GS_ord(0),e_GS_ord(1));
        fprintf(err_dump, "\n");
        e_GS.replaceRow(e_GS_ord,ord-1);

        closefile(GS_dump, SoluGSmodes); 
        closefile(GS_dump_v, SolvGSmodes); 
        closefile(GSstat_dump, SoluGSstat); 
        closefile(GSstat_dump_v, SolvGSstat); 
    }

    // AAPG
    printf("\nAAPG...\n");
    tt.tick();
    PCSet myPCSet("ISP",ord_AAPG_GS,dim,pcType,0.0,1.0,false); 
    tt.tock("Took");
   
    Array2D<double> e_AAPG(ord_AAPG,2,0.e0); 
    Array1D<double> t_AAPG = AAPG(dof, inpParams, fbar, dTym, ord_AAPG_GS, pcType, noutput, dim, nStep, scaledKLmodes, myPCSet, factor_OD, ord_AAPG, act_D, p, MCS_s_dis, err_dump, sample_mstd_2D, samPts_norm, e_AAPG);
    
    // output the timing
    Array1D<double> t(3+ord_GS+ord_AAPG,0.e0);
    t(0) = t_MCS;
    for (int i=0;i<ord_GS;i++)
	    t(i+1)=t_GS(i);
    for (int i=0;i<ord_AAPG+2;i++)
        t(i+1+ord_GS)=t_AAPG(i);

    // output the error info
    cout << "Printing the error of GS in displacement..." << endl;
    for(int i=0;i<ord_GS;i++){
        cout << "GS ord " << i+1<< " is em="<<e_GS(i,0) <<", es= "<< e_GS(i,1)<< endl;
    }
    cout << "Printing the error of AAPG in displacement..." << endl;
    for(int i=0;i<ord_AAPG;i++){
        cout << "AAPG ord " << i+1<< " is em="<<e_AAPG(i,0) <<", es= "<< e_AAPG(i,1)<<endl;
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

