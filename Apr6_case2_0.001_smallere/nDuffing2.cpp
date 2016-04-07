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
#include "nMCS.h"
#include "nGhanemSpanos.h"
#include "AAPG.h"
#include "nAAPG.h"
#include "ticktock.h"

int main(int argc, char *argv[]){

    int dof=2;
    int ord_GS=2;
    int ord_AAPG=2;
    int ord_AAPG_GS=2;
    int refine = 1;
    bool act_D = false;
    Array1D<double> time(1+ord_GS,0.e0);
    int nkl=8;
    int dim=nkl+6*dof;// set epsilon to be stochastic coeffs on each dof
    int noutput=2;
    int nspl =10000;
    int factor_OD = 0.99;
    string pcType="LU";  //PC type
    double dTym = 0.01;
    Array1D<double> initial(2*dof,0.e0); // initial condition
    Array1D<double> initial_sigma(2*dof,0.e0);
    for (int i=0;i<dof;i++){
        initial_sigma(i)=0.5;
        initial_sigma(dof+i)=0.1;
    } 

    // epsilon
    //Array1D<double>  epsilon_mean(dof,1e4);
    double e1 = 1.0;
    double e2 = 0.1;
    /* Read the user input */
    int c;
    while ((c=getopt(argc,(char **)argv,"r:G:d:e:m:N:"))!=-1){
        switch (c) {
        case 'r':
            refine = (strtod(optarg, (char **)NULL));
            break;
        case 'G':
            ord_AAPG_GS = (strtod(optarg, (char **)NULL));
            break;
        case 'd':
            dTym = 1/(strtod(optarg, (char **)NULL));
            break;
        case 'e':
            e1 = strtod(optarg, (char **)NULL);
            break;
        case 'm':
            e2 = strtod(optarg, (char **)NULL);
            break;
        case 'N':
            nspl = strtod(optarg, (char **)NULL);
            break;
        }
    }
    e1=e1;
    cout << "epsilon_mean=" << e1 << endl;
    cout << "e_sigma=" << e2 << endl;
    cout << "nspl=" << nspl << endl;
    cout << "dTym=" << dTym << endl;
    cout << "ord_AAPG_GS=" << ord_AAPG_GS << endl;
    cout << "AAPG1 is refined by "<< refine << " times" << endl;
    Array1D<double>  epsilon_mean(dof,e1);
    Array1D<double>  e_sigma(dof,e2);

    // Time marching info
    double tf = 10;
    // Number of steps
    int nStep=(int) tf / dTym;
    int nStep_fine=nStep*refine;

    // MCK
    Array1D<Array2D<double> > mck(3);
    // m
    //Array1D<double> temp_m(dof,1e4);
    Array2D<double> temp_m(dof,2,0.0);
    for (int i=0;i<dof;i++){
        temp_m(i,0)=1.0;
        temp_m(i,1)=0.1;
    }
    mck(0) = temp_m;
    // k
    //Array1D<double> temp_k(dof,4e7);
    Array2D<double> temp_k(dof,2,0.0);
    for (int i=0;i<dof;i++){
        //temp_k(i,0)=1.0-(i-i%(dof/3))/(dof/3)*0.1;
        temp_k(i,0)=1.0;
        temp_k(i,1)=0.1;
    }
    write_datafile(temp_k,"temp_k.dat");
    //temp_k(4)=0.9;
    //temp_k(5)=0.9;
    //temp_k(6)=0.9;
    //temp_k(7)=0.8;
    //temp_k(8)=0.8;
    //temp_k(9)=0.8;
    mck(2) = temp_k;
    // c
    Array2D<double> temp_c(dof,2,0.e0);
    for (int i=0;i<dof;i++){
        temp_c(i,0)=2*0.1*sqrt(temp_m(i,0)*temp_k(i,0));
        temp_c(i,1)=0.1*temp_c(i,0);
    }
    mck(1) = temp_c;

    // sample point
    Array2D<double> samPts_norm(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",0,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts_norm);
    //write_datafile(samPts_norm,"samPts_norm.dat");

    // force
    Array2D<double> scaledKLmodes(2*nStep+1,nkl,0.e0);
    Array2D<double> scaledKLmodes_fine(2*nStep_fine+1,nkl,0.e0);
    if (nkl>0){
        cout << "Generating KL..." << endl;
        double clen = 0.1;
        double sigma=0.8;
        char* cov_type = (char *)"Exp";
        //genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
        genKL(scaledKLmodes_fine, 2*nStep_fine+1, nkl, clen, sigma, tf, cov_type);
    }
    Array1D<double> fbar(2*nStep+1,0.e0);//mean of forcing
    Array1D<double> fbar_fine(2*nStep_fine+1,0.e0);//mean of forcing
    double t_temp = 0.0; 
    int i_temp = 0;
    cout << "Generating fbar..." << endl;
    for (int i=0;i<2*nStep_fine+1;i++){
        fbar_fine(i) = 2.0-2.0*sin(2*3.1415926*t_temp)*exp(-0.3*t_temp);
        //fbar_fine(i) = 2.0-sin(40*3.1415926*t_temp)*exp(-0.3*t_temp);
        //fbar_fine(i)=2.0;
        if (i%refine == 0){
            fbar(i_temp)=fbar_fine(i);
            Array1D<double> temp(nkl,0.e0);
            getRow(scaledKLmodes_fine,i,temp); 
            scaledKLmodes.replaceRow(temp,i_temp);
            i_temp++;
        }
        t_temp +=dTym/refine/2;
    }
    write_datafile_1d(fbar,"nfbar.dat");
    write_datafile_1d(fbar_fine,"nfbar_fine.dat");
    write_datafile(scaledKLmodes,"nKL.dat");
    write_datafile(scaledKLmodes_fine,"nKL_fine.dat");

    /////////////---MCS--/////////////
    cout << "Starting MCS..." << endl;
    Array1D<Array2D<double> > result_MCS(2*dof);
    for (int idof=0;idof<2*dof;idof++){
        Array2D<double> temp_result(nStep+1,nspl);
        result_MCS(idof)=temp_result;
    }
    Array2D<double> epsilon_MCS_samples(nspl,dof,0.e0);
    Array2D<double> initial_MCS_samples(nspl,2*dof,0.e0);
    Array2D<double> mck_MCS_samples(nspl,3*dof,0.e0);
    Array1D<Array2D<double> >f_MCS_samples(dof);
    for (int i=0;i<dof;i++){
        Array2D<double> f_MCS(nspl,2*nStep+1,0.e0);
        f_MCS_samples(i)=f_MCS;
    }
    TickTock tt;
    tt.tick();
    #pragma omp parallel default(none) shared(mck,result_MCS,dof,nStep,nspl,samPts_norm,initial,initial_sigma,initial_MCS_samples,dTym,fbar,nkl,epsilon_mean,scaledKLmodes,e_sigma,epsilon_MCS_samples,mck_MCS_samples,f_MCS_samples) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        // initialize epsilon, initial and mck
        Array1D<double> epsilon_MCS(epsilon_mean);
        Array1D<double> initial_MCS(initial);
        // sample epsilon
        for (int i=0;i<dof;i++){
            epsilon_MCS(i)+=samPts_norm(iq,nkl+i)*e_sigma(i); 
            epsilon_MCS_samples(iq,i)=epsilon_MCS(i);
        }
        // sample initial conditions
        for (int i=0;i<2*dof;i++){
            initial_MCS(i)+=samPts_norm(iq,nkl+dof+i)*initial_sigma(i);
            initial_MCS_samples(iq,i)=initial_MCS(i);
        }
        Array1D<double> m(dof,0.e0);
        getCol(mck(0),0,m);    
        Array1D<double> c(dof,0.e0);
        getCol(mck(1),0,c);    
        Array1D<double> k(dof,0.e0);
        getCol(mck(2),0,k);    
        // sample force using the mean of the mass
        Array2D<double> force_MCS=nsample_force(dof,samPts_norm,iq,fbar,nkl,scaledKLmodes,m); 
        // reorganize force for output
        for (int i=0;i<dof;i++)
            for (int ix=0;ix<2*nStep+1;ix++)
                f_MCS_samples(i)(iq,ix)=force_MCS(ix,i);
        // sample mck
        Array1D<Array1D<double> > mck_MCS(3);
        for (int i=0;i<dof;i++){
            m(i)+=samPts_norm(iq,nkl+3*dof+i)*mck(0)(i,1);
            c(i)+=samPts_norm(iq,nkl+4*dof+i)*mck(1)(i,1);
            k(i)+=samPts_norm(iq,nkl+5*dof+i)*mck(2)(i,1);
            mck_MCS_samples(iq,i)=m(i);
            mck_MCS_samples(iq,i+dof)=c(i);
            mck_MCS_samples(iq,i+2*dof)=k(i);
        }
        mck_MCS(0)=m;
        mck_MCS(1)=c;
        mck_MCS(2)=k;
        Array1D<Array1D<double> > temp=ndet(dof,nStep,dTym,force_MCS,epsilon_MCS,mck_MCS,initial_MCS); 
        for (int idof=0;idof<2*dof;idof++){
            for (int ix=0;ix<nStep+1;ix++){
                //result_MCS(idof).replaceCol(temp(idof),iq);
                result_MCS(idof)(ix,iq)=temp(idof)(ix);
            }
        }
    }
    }
    time(0) = tt.silent_tock();
    tt.tock("Took");
    //write_datafile(epsilon_MCS_samples,"epsilon_samples.dat");
    //write_datafile(initial_MCS_samples,"initial_samples.dat");
    // examine statistics of epsilon & initial conditions
    Array2D<double> stat_e(dof,2,0.e0);
    Array2D<double> stat_m(dof,2,0.e0);
    Array2D<double> stat_c(dof,2,0.e0);
    Array2D<double> stat_k(dof,2,0.e0);
    for (int i=0;i<dof;i++){
        Array1D<double> epsilon_sample(nspl,0.e0);
        getCol(epsilon_MCS_samples,i,epsilon_sample); 
        mstd_MCS(epsilon_sample,stat_e(i,0),stat_e(i,1));
        cout << "Mean of epsilon on dof " << i << " is " << stat_e(i,0) << ". Std is " << stat_e(i,1) << "." << endl;
        Array1D<double> m_sample(nspl,0.e0);
        getCol(mck_MCS_samples,i,m_sample); 
        mstd_MCS(m_sample,stat_m(i,0),stat_m(i,1));
        cout << "Mean of m on dof " << i << " is " << stat_m(i,0) << ". Std is " << stat_m(i,1) << "." << endl;
        Array1D<double> c_sample(nspl,0.e0);
        getCol(mck_MCS_samples,i+dof,c_sample); 
        mstd_MCS(c_sample,stat_c(i,0),stat_c(i,1));
        cout << "Mean of c on dof " << i << " is " << stat_c(i,0) << ". Std is " << stat_c(i,1) << "." << endl;
        Array1D<double> k_sample(nspl,0.e0);
        getCol(mck_MCS_samples,i+dof*2,k_sample); 
        mstd_MCS(k_sample,stat_k(i,0),stat_k(i,1));
        cout << "Mean of k on dof " << i << " is " << stat_k(i,0) << ". Std is " << stat_k(i,1) << "." << endl;
    }
    Array2D<double> stat_i(2*dof,2,0.e0);
    for (int i=0;i<2*dof;i++){
        Array1D<double> initial_sample(nspl,0.e0);
        getCol(initial_MCS_samples,i,initial_sample); 
        mstd_MCS(initial_sample,stat_i(i,0),stat_i(i,1));
        cout << "Mean on initial condition on dof " << i << " is " << stat_i(i,0) << ". Std is " << stat_i(i,1) << "." << endl;
    }
    Array2D<double> f_mean(dof,2*nStep+1,0.e0);
    Array2D<double> f_std(dof,2*nStep+1,0.e0);
    for (int i=0;i<dof;i++){
        for (int ix=0;ix<2*nStep+1;ix++){
            Array1D<double> f_sample(nspl,0.e0);
            getCol(f_MCS_samples(i),ix,f_sample);
            mstd_MCS(f_sample,f_mean(i,ix),f_std(i,ix));
        }
    }
    write_datafile(f_mean,"f_mean.dat");
    write_datafile(f_std,"f_std.dat");

    // post-process result
    Array2D<double> mean_MCS(nStep+1,2*dof,0.e0);
    Array2D<double> std_MCS(nStep+1,2*dof,0.e0);
    #pragma omp parallel default(none) shared(nStep,dof,nspl,result_MCS,mean_MCS,std_MCS)
    {
        #pragma omp for 
        for (int ix=0;ix<nStep+1;ix++){
            for (int i=0;i<2*dof;i++){
                Array1D<double> temp_result(nspl,0.e0);
                getRow(result_MCS(i),ix,temp_result);        
                mstd_MCS(temp_result,mean_MCS(ix,i),std_MCS(ix,i));
            }
        }
    }
    for (int ix=0;ix<nStep+1;ix++){
        // report to screen
        if (ix % ((int) nStep/noutput) == 0){
            for (int i=0;i<dof;i++){
                cout << "dof " << i << endl;
                WriteMeanStdDevToStdOut(ix, ix*dTym, mean_MCS(ix,i), std_MCS(ix,i));
            }
        }
    }
    ostringstream nameMCS;
    nameMCS << "mean_MCS_eps" << e1 << ".dat";
    string nameMCS_str = nameMCS.str();
    write_datafile(mean_MCS,nameMCS_str.c_str());
    //write_datafile(mean_MCS,"mean_MCS_n.dat");
    write_datafile(std_MCS,"std_MCS_n.dat");
   
    /////////////---GS---///////////// 
    //Array1D<Array2D<double> > mstd_MCS(dof);
    cout << "Starting GS..." << endl;
    //Array2D<double> e_GS(ord_GS,)
    Array1D<Array1D<double> > initial_GS(dof); // Initial conditions
    Array1D<Array1D<double> > epsilon_GS(dof);
    Array1D<Array1D<Array1D<double> > > mck_GS(3);
    for (int ord=1;ord<=ord_GS;ord++){
        // Generate PCSet
        TickTock tt;
        tt.tick();
	    PCSet myPCSet("ISP",ord,dim,pcType,0.0,1.0); 
        tt.tock("Took");
        int nPCTerms = myPCSet.GetNumberPCTerms();
        // Prepare the force in PC format
        Array1D<Array2D<double> > f_GS(dof);
        for (int idof=0;idof<dof;idof++){
            Array2D<double> temp_f(2*nStep+1,nPCTerms,0.e0);
            f_GS(idof)=temp_f;
            for (int ix=0;ix<2*nStep+1;ix++){
                f_GS(idof)(ix,0)=fbar(ix)*mck(0)(idof,0);
                for (int i=0;i<nkl;i++){
                    f_GS(idof)(ix,i+1)=scaledKLmodes(ix,i)*mck(0)(idof,0);
                }
            }
        }
        //epsilon
        for (int i=0;i<dof;i++){
            Array1D<double> temp_epsilon(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_e(i,0),stat_e(i,1),nkl+i+1,temp_epsilon);
            epsilon_GS(i)=temp_epsilon;
        }
        // initial conditions
        for (int i=0;i<dof;i++){
            Array1D<double> temp_init(2*nPCTerms,0.e0);
            Array1D<double> temp_init2(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_i(i,0),stat_i(i,1),1+nkl+dof+i,temp_init2);
            Array1D<double> temp_init3(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_i(i+dof,0),stat_i(i+dof,1),1+nkl+2*dof+i,temp_init3);
            merge(temp_init2,temp_init3,temp_init); 
            initial_GS(i)= temp_init;
        }
        // mck
        Array1D<Array1D<double> > m_GS(dof);
        Array1D<Array1D<double> > c_GS(dof);
        Array1D<Array1D<double> > k_GS(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> temp_m(nPCTerms,0.e0);
            Array1D<double> temp_c(nPCTerms,0.e0);
            Array1D<double> temp_k(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_m(i,0),stat_m(i,1),1+nkl+3*dof+i,temp_m);
            myPCSet.InitMeanStDv(stat_c(i,0),stat_c(i,1),1+nkl+4*dof+i,temp_c);
            myPCSet.InitMeanStDv(stat_k(i,0),stat_k(i,1),1+nkl+5*dof+i,temp_k);
            m_GS(i)=temp_m;
            c_GS(i)=temp_c;
            k_GS(i)=temp_k;
        }
        mck_GS(0)=m_GS;
        mck_GS(1)=c_GS;
        mck_GS(2)=k_GS;
        // initialize solution
        Array1D<Array2D<double> > uv_solution(dof);
        for (int i=0;i<dof;i++){
            Array2D<double> temp_solution(nStep,2*nPCTerms,0.e0);
            uv_solution(i) = temp_solution;
        }
        cout << "Starting GS order " << ord << endl;
        tt.tick();
        nGS(dof, myPCSet, epsilon_GS, mck_GS, nStep, initial_GS, dTym, f_GS, uv_solution);
        time(ord)=tt.silent_tock();
        cout << "Order " << ord << " finished." <<endl; 
        // Post-process the solution
        Array2D<double> e2(dof,4,0.e0);
        Array2D<double> e_GS = postprocess_nGS(dof,nStep,uv_solution,myPCSet,dTym,ord,mean_MCS,std_MCS,e2);
        // print out the error
        cout << "Error kind 1 is" << endl;
        for (int i=0;i<dof;i++){
            cout << "Dof " << i << ", m_v:" << e_GS(i,0) << ",s_v:" << e_GS(i,1) << "," <<",m_u:" << e_GS(i,2) << ",s_u:"<<e_GS(i,3) << "." << endl;
        }
        cout << "Error kind 2 is" << endl;
        for (int i=0;i<dof;i++){
            cout << "Dof " << i << ", m_v:" << e2(i,0) << ",s_v:" << e2(i,1) << "," <<",m_u:" << e2(i,2) << ",s_u:"<<e2(i,3) << "." << endl;
        }
        ostringstream name;
        name << "e_GS_" << ord << "_d_"<<dTym<<".dat";
        string name_str = name.str();
        write_datafile(e_GS,name_str.c_str());
    }

    write_datafile_1d(time,"t_nMCSGS.dat");

    ////////---AAPG----///////
    // Compute multiindex
    Array2D<int> Pbtot;
    computeMultiIndex(dim,ord_AAPG_GS,Pbtot);
    // Initialize pcbasis
    PCBasis p_basis(pcType, 0.0, 1.0, ord_AAPG_GS);
    // Get the 1d norms-squared
    Array1D<double> norms1d;
    p_basis.Get1dNormsSq(norms1d);
    // Compute normsq
    Array1D<double> normsq(Pbtot.XSize(),1.e0);
    // For each term, multiply appropriate 1d norms-squared
    for(unsigned int ipc=0; ipc<Pbtot.XSize(); ipc++)
        for(int id=0; id<dim; id++)
            normsq(ipc) *= norms1d(Pbtot(ipc,id));
   
    cout << "AAPG..." << endl;
    nAAPG(refine,dof,nkl,dim,nStep,ord_AAPG_GS,noutput,factor_OD,ord_AAPG,act_D,fbar,fbar_fine,dTym,epsilon_mean,pcType,scaledKLmodes,scaledKLmodes_fine,stat_e,stat_i,stat_m,stat_c,stat_k,normsq,mean_MCS,std_MCS,mck);

    return 0;
}
