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
#include "ticktock.h"

int main(int argc, char *argv[]){

    int dof=5;
    int ord_GS=2;
    int nkl=5;
    int dim=nkl+3*dof;// set epsilon to be stochastic coeffs on each dof
    int noutput=2;
    int nspl =1000;
    string pcType="LU";  //PC type
    Array1D<double> initial(2*dof,0.e0); // initial condition
    Array1D<double> initial_sigma(2*dof,0.e0);
    for (int i=0;i<dof;i++)
        initial_sigma(i)=0.1; 

    // epsilon
    Array1D<double>  epsilon_mean(dof,1e4);
    Array1D<double>  e_sigma(dof,1e4);

    // Time marching info
    double dTym = 0.01;
    double tf = 10;
    // Number of steps
    int nStep=(int) tf / dTym;

    // MCK
    Array1D<Array1D<double> > mck(3);
    // m
    Array1D<double> temp_m(dof,1e4);
    mck(0) = temp_m;
    // k
    Array1D<double> temp_k(dof,4e7);
    temp_k(4)=temp_k(4)*0.9;
    temp_k(5)=temp_k(5)*0.9;
    temp_k(6)=temp_k(6)*0.9;
    temp_k(7)=temp_k(7)*0.8;
    temp_k(8)=temp_k(8)*0.8;
    temp_k(9)=temp_k(9)*0.8;
    mck(2) = temp_k;
    // c
    Array1D<double> temp_c(dof,0.e0);
    for (int i=0;i<dof;i++)
        temp_c(i)=2*0.04*sqrt(temp_m(i)*temp_k(i));
    mck(1) = temp_c;

    // sample point
    Array2D<double> samPts_norm(nspl,dim,0.e0);
    PCSet MCPCSet("NISPnoq",ord_GS,dim,pcType,0.0,1.0);
    MCPCSet.DrawSampleVar(samPts_norm);
    write_datafile(samPts_norm,"samPts_norm.dat");

    // force
    Array1D<Array2D<double> > f_GS(dof);
    Array2D<double> scaledKLmodes(2*nStep+1,nkl,0.e0);
    double clen = 0.05;
    double sigma=0.5;
    char* cov_type = (char *)"Exp";
    genKL(scaledKLmodes, 2*nStep+1, nkl, clen, sigma, tf, cov_type);
    Array1D<double> fbar(2*nStep+1,0.e0);//mean of forcing
    double t_temp = 0.0; 
    for (int i=0;i<2*nStep+1;i++){
        fbar(i) = 0.2*(2.0-sin(2*3.1415926*t_temp)*exp(-0.3*t_temp));
        t_temp +=dTym/2;
    }
    write_datafile_1d(fbar,"nfbar.dat");
    write_datafile(scaledKLmodes,"nKL.dat");

    /////////////---MCS--/////////////
    cout << "Starting MCS..." << endl;
    Array1D<Array2D<double> > result_MCS(2*dof);
    for (int idof=0;idof<2*dof;idof++){
        Array2D<double> temp_result(nStep+1,nspl);
        result_MCS(idof)=temp_result;
    }
    Array2D<double> epsilon_MCS_samples(nspl,dof,0.e0);
    Array2D<double> initial_MCS_samples(nspl,2*dof,0.e0);
    TickTock tt;
    tt.tick();
    #pragma omp parallel default(none) shared(mck,result_MCS,dof,nStep,nspl,samPts_norm,initial,initial_sigma,initial_MCS_samples,dTym,fbar,nkl,epsilon_mean,scaledKLmodes,e_sigma,epsilon_MCS_samples) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array2D<double> force_MCS=nsample_force(dof,samPts_norm,iq,fbar,nkl,scaledKLmodes,mck); 
        Array1D<double> epsilon_MCS(epsilon_mean);
        Array1D<double> initial_MCS(initial);
        for (int i=0;i<dof;i++){
            epsilon_MCS(i)+=samPts_norm(iq,nkl+i)*e_sigma(i); 
            epsilon_MCS_samples(iq,i)=epsilon_MCS(i);
        }
        for (int i=0;i<2*dof;i++){
            initial_MCS(i)+=samPts_norm(iq,nkl+dof+i)*initial_sigma(i);
            initial_MCS_samples(iq,i)=initial_MCS(i);
        }
        Array1D<Array1D<double> > temp=ndet(dof,nStep,dTym,force_MCS,epsilon_MCS,mck,initial_MCS); 
        for (int idof=0;idof<2*dof;idof++){
            for (int ix=0;ix<nStep+1;ix++){
                //result_MCS(idof).replaceCol(temp(idof),iq);
                result_MCS(idof)(ix,iq)=temp(idof)(ix);
            }
        }
    }
    }
    write_datafile(epsilon_MCS_samples,"epsilon_samples.dat");
    write_datafile(initial_MCS_samples,"initial_samples.dat");
    tt.tock("Took");
    // examine statistics of epsilon & initial conditions
    Array2D<double> stat_e(dof,2,0.e0);
    for (int i=0;i<dof;i++){
        Array1D<double> epsilon_sample(nspl,0.e0);
        getCol(epsilon_MCS_samples,i,epsilon_sample); 
        mstd_MCS(epsilon_sample,stat_e(i,0),stat_e(i,1));
        cout << "Mean of epsilon on dof " << i << " is " << stat_e(i,0) << ". Std is " << stat_e(i,1) << "." << endl;
    }
    Array2D<double> stat_i(2*dof,2,0.e0);
    for (int i=0;i<2*dof;i++){
        Array1D<double> initial_sample(nspl,0.e0);
        getCol(initial_MCS_samples,i,initial_sample); 
        mstd_MCS(initial_sample,stat_i(i,0),stat_i(i,1));
        cout << "Mean on initial condition on dof " << i << " is " << stat_i(i,0) << ". Std is " << stat_i(i,1) << "." << endl;
    }

    // post-process result
    Array2D<double> mean_MCS(nStep+1,2*dof,0.e0);
    Array2D<double> std_MCS(nStep+1,2*dof,0.e0);
    for (int ix=0;ix<nStep+1;ix++){
        for (int i=0;i<2*dof;i++){
            Array1D<double> temp_result(nspl,0.e0);
            getRow(result_MCS(i),ix,temp_result);        
            mstd_MCS(temp_result,mean_MCS(ix,i),std_MCS(ix,i));
        }
        // report to screen
        if (ix % ((int) nStep/noutput) == 0){
            for (int i=0;i<dof;i++){
                cout << "dof " << i << endl;
                WriteMeanStdDevToStdOut(ix, ix*dTym, mean_MCS(ix,i), std_MCS(ix,i));
            }
        }
    }
    write_datafile(mean_MCS,"mean_MCS_n.dat");
    write_datafile(std_MCS,"std_MCS_n.dat");
   
    /////////////---GS---///////////// 
    //Array1D<Array2D<double> > mstd_MCS(dof);
    cout << "Starting GS..." << endl;
    //Array2D<double> e_GS(ord_GS,)
    Array1D<Array1D<double> > initial_GS(dof); // Initial conditions
    Array1D<Array1D<double> > epsilon_GS(dof);
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
                f_GS(idof)(ix,0)=fbar(ix)*mck(0)(idof);
                for (int i=0;i<nkl;i++){
                    //Array1D<double> tempf(2*nStep+1,0.e0);
                    //getCol(scaledKLmodes,i,tempf);
                    f_GS(idof)(ix,i+1)=scaledKLmodes(ix,i)*mck(0)(idof);
                }
            }
        }
        //epsilon
        for (int i=0;i<dof;i++){
            Array1D<double> temp_epsilon(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_e(i,0),stat_e(i,1),nkl+i+1,temp_epsilon);
            epsilon_GS(i)=temp_epsilon;
        }
        // initial conditions set to zero for now
        for (int i=0;i<dof;i++){
            Array1D<double> temp_init(2*nPCTerms,0.e0);
            //initial_GS(i)(0)=initial(i);   
            //initial_GS(i)(nPCTerms)=initial(dof+i);   
            Array1D<double> temp_init2(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_i(i,0),stat_i(i,1),nkl+dof+i,temp_init2);
            Array1D<double> temp_init3(nPCTerms,0.e0);
            myPCSet.InitMeanStDv(stat_i(i+dof,0),stat_i(i+dof,1),nkl+2*dof+i,temp_init3);
            merge(temp_init2,temp_init3,temp_init); 
            initial_GS(i)= temp_init;
        }
        // initialize solution
        Array1D<Array2D<double> > uv_solution(dof);
        for (int i=0;i<dof;i++){
            Array2D<double> temp_solution(nStep,2*nPCTerms,0.e0);
            uv_solution(i) = temp_solution;
        }
        cout << "Starting GS order " << ord << endl;
        nGS(dof, myPCSet, epsilon_GS, mck, nStep, initial_GS, dTym, f_GS, uv_solution);
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
        //for (int i=0;i<dof;i++){
        //    ostringstream name2;
        //    name2 << "sol_GS_" << ord << "dof_"<< i <<".dat";
        //    string name2_str = name2.str();
        //    write_datafile(uv_solution(i),name2_str.c_str());
        //}
        ostringstream name;
        name << "e_GS_" << ord << ".dat";
        string name_str = name.str();
        write_datafile(e_GS,name_str.c_str());
    }

    return 0;
}
