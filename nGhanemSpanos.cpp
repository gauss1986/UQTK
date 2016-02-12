#include <math.h>
#include <omp.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "Utils.h"

void nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& initial,  Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck, Array1D<Array1D<double> >& f, double dTym, int nStep, Array1D<Array1D<double> >& solution){
    // n-dimensional Duffing oscillator

    // deterministic M and C
    Array2D<double> M(dof,dof,0.e0);
    Array2D<double> C(dof,dof,0.e0);
    for (int i=0;i<dof;i++){
        M(i,i)=mck(0)(i);
        C(i,i)=mck(1)(i);
    }
    for (int i=0;i<dof-1;i++){
        C(i,i+1)=-mck(1)(i+1);
        C(i+1,i)=-mck(1)(i+1);
        C(i,i)+=mck(1)(i+1);
    }

     
}

void forward_duffing_GS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& epsilon, Array1D<double>& m, Array2D<double>& C){
    int nPC = myPCSet.GetNumberPCTerms();

    // kbar 
    Array1D<Array1D<double> > kbar(dof);
    // kbar0
    Array1D<double> temp_u(nPC,0.e0);
    Array1D<double> temp_k(nPC,0.e0);
    myPCSet.IPow(u(0),2,temp_u);
    myPCSet.Prod(temp_u,epsilon(0),temp_k);
    temp_k(0)=temp_k(0)+1; 
    prodVal(temp_k,mck(2)(0));
    kbar(0) = temp_k;
    // kbar rest 
    for (int i=1;i<dof;i++){
        Array1D<double> temp_k1(nPC,0.e0);
        Array1D<double> temp_u1(u(i));
        subtractVec(u(i-1),temp_u1); 
        myPCSet.IPow(temp_u1,2,temp_u);
        myPCSet.Prod(temp_u,epsilon(i),temp_k1);
        temp_k1(0)=temp_k1(0)+1; 
        prodVal(temp_k1,mck(2)(i));
        kbar(i) = temp_k1;
    }

    // K
    Array2D<double> K(dof,dof);
    Array1D<double> temp_K(nPC,0.e0);
    for (int i=0;i<dof;i++)
        for (int j=0;j<dof;j++)
            K(i,j)=temp_K;
    for int i=0;i<dof-1;i++){
        addVec(kbar(i),K(i,i));
        addVec(kbar(i+1),K(i,i));
        prodVal(kbar(i+1),-1);
        K(i,i+1)=kbar(i+1);
        K(i+1,i)=kbar(i+1);
    }
    K(dof,dof)=kbar(dof);

    // Compute u at next step
    for (int i=1;i<dof-1;i++){
        Array1D<double> temp_u2(ndof,0.e0);
        Array1D<double> u_minus(u(i-1));
        Array1D<double> u_current(u(i));
        Array1D<double> u_plus(u(i+1));
        prodVal(u_minus,m(i));
                
    }

    return;
}

void GS(int dof, PCSet& myPCSet, Array1D<int>& coeff_D, int nPCTerms, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS, Array1D<Array2D<double> >& solution){

    // Initialize working variable
    Array1D<Array1D<double> > result(dof);
    for (int i=0;i<dof;i++){
        Array2D<double> temp(nStep+1,nPCTerms,0.e0);
        solution(i) = temp;
        solution(i).replaceRow(initial(i),0);
        Array1D<double> temp2(initial(i));
        result(i)=temp2;
    }
    
    // Forward run
    Array1D<double> force_current(nPCTerms,0.e0);
    Array1D<double> force_mid(nPCTerms,0.e0);
    Array1D<double> force_plus(nPCTerms,0.e0);
    for (int ix=0;ix<nStep;ix++){
        // Subtract the force at currect, mid and next step
        getRow(f_GS,2*ix,force_current);
        getRow(f_GS,2*ix+1,force_mid);
        getRow(f_GS,2*ix+2,force_plus);
        // Step forward 
        forward_duffing_GS(dof, myPCSet, inpParams, force_current, force_mid, force_plus, dTym, result, coeff_D);
        // Update solution
        for (int i=0;i<dof;i++){
            solution(i).replaceRow(result(i),ix+1);
        }
    }    
}

void forward_duffing_GS(int dof, PCSet& myPCSet, Array1D<double>& inpParams, Array1D<double>& force_current, Array1D<double>& force_mid, Array1D<double>& force_plus,  double dTym, Array1D<Array1D<double> >& x, Array1D<int>& coeff_D){
        //Save solution at current time step
        Array1D<Array1D<double> > x0(x.XSize());
        for (unsigned int i=0;i<x.XSize();i++){
            Array1D<double> temp(x(i));
            x0(i)=temp;
        }
        
        // Integrate with classical 4th order Runge-Kutta
        // k1
        Array1D<Array1D<double> > dev1(dof);
        RHS_GS(dof, myPCSet, force_current, x, inpParams, dev1, coeff_D);
        for (unsigned int i=0;i<x.XSize();i++)
            addVecAlphaVecPow(x(i),0.5*dTym,dev1(i),1);

        // k2
        Array1D<Array1D<double> > dev2(dof);
        RHS_GS(dof, myPCSet, force_mid, x, inpParams,dev2, coeff_D);
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i);
            addVecAlphaVecPow(x(i),0.5*dTym,dev2(i),1);
        }

        // k3
        Array1D<Array1D<double> > dev3(dof);
        RHS_GS(dof, myPCSet, force_mid, x, inpParams,dev3,coeff_D);
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i);
            addVecAlphaVecPow(x(i),dTym,dev3(i),1);
        }

        // k4
        Array1D<Array1D<double> > dev4(dof);
        RHS_GS(dof, myPCSet, force_plus, x, inpParams,dev4,coeff_D);

        // Advance to next time step
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i);
            prodVal(dev2(i),2);
            prodVal(dev3(i),2);
            addVec(dev3(i),dev4(i));
            addVec(dev2(i),dev4(i));
            addVec(dev1(i),dev4(i));
            addVecAlphaVecPow(x(i),dTym/6,dev4(i),1);
        }        
   }
   
void RHS_GS(int dof, PCSet& myPCSet, Array1D<double>& force, Array1D<Array1D<double> >& x, Array1D<double>& inpParams, Array1D<Array1D<double> >& dev, Array1D<int>& coeff_D){
        for (int i=0;i<dof;i++){
            Array1D<double> temp(x(i));
            dev(i) = temp;
        }
        int nPCTerms = x(0).Length();
        
        if ((abs(inpParams(0))<1e-10)){ //Duffing
            // parse input parameters
            Array1D<double> epsilon(nPCTerms,0.e0); 
            Array1D<double> zeta(nPCTerms,0.e0); 
            if (coeff_D(0)+1<nPCTerms){
                myPCSet.InitMeanStDv(inpParams(1),inpParams(3),coeff_D(0)+1,zeta);
            }
            if (coeff_D(1)+1<nPCTerms){
                myPCSet.InitMeanStDv(inpParams(2),inpParams(4),coeff_D(1)+1,epsilon);
            }
            // buff to store temperary results 
            Array1D<double> temp(nPCTerms,0.e0);                     
            Array1D<double> temp2(nPCTerms,0.e0);                     
            // nonlinear term
            myPCSet.IPow(x(1),temp,3);
            myPCSet.Prod(epsilon,temp,temp2);
            myPCSet.Prod(zeta,x(0),temp);
            for (int ind=0;ind<nPCTerms;ind++){
                dev(1)(ind) = x(0)(ind);//velocity term
                dev(0)(ind) = force(ind)-temp2(ind)-2*temp(ind)-x(1)(ind);
            }
        }
}

Array1D<double> postprocess_GS(int noutput, int nPCTerms, int nStep,  Array2D<double>& solution, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS, Array2D<double>& stat){
    // Output solution (mean and std) 
    Array1D<double> temp(nPCTerms,0.e0);
    Array1D<double> StDv(nStep+1,0.e0);
    for (int i=0;i<nStep+1;i++){
        getRow(solution,i,temp);
        StDv(i) = myPCSet.StDv(temp);
        // Write time and solution to file
        WriteModesToFilePtr(i, temp.GetArrayPointer(), nPCTerms, GS_dump);
        // Write dis (mean and std) to file and screen
        WriteMeanStdDevToFilePtr(i,temp(0),StDv(i),GSstat_dump);
        if (i % ((int) nStep/noutput) == 0){
            WriteMeanStdDevToStdOut(i,i*dTym,temp(0),StDv(i));
        }
    }

    Array1D<double> mean;
    getCol(solution,0,mean);
    Array1D<double> e =  error(mean, StDv, mstd_MCS);

    stat.replaceRow(mean,0);
    stat.replaceRow(StDv,1);    
    return e;
}
