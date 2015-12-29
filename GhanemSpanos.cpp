#include <math.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "Utils.h"

void GS(int dof, PCSet& myPCSet, int order, int dim, int nPCTerms, string pcType, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS, Array1D<Array2D<double> >& solution){

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
        forward_duffing_GS(dof, myPCSet, inpParams, force_current, force_mid, force_plus, dTym, result);
        // Update solution
        for (int i=0;i<dof;i++){
            solution(i).replaceRow(result(i),ix+1);
        }
    }    
}

void forward_duffing_GS(int dof, PCSet& myPCSet, Array1D<double>& inpParams, Array1D<double>& force_current, Array1D<double>& force_mid, Array1D<double>& force_plus,  double dTym, Array1D<Array1D<double> >& x){
        //Save solution at current time step
        Array1D<Array1D<double> > x0(x.XSize());
        for (unsigned int i=0;i<x.XSize();i++){
            Array1D<double> temp(x(i));
            x0(i)=temp;
        }
        
        // Integrate with classical 4th order Runge-Kutta
        // k1
        Array1D<Array1D<double> > dev1(dof);
        RHS_GS(dof, myPCSet, force_current, x, inpParams, dev1);
        for (unsigned int i=0;i<x.XSize();i++)
            addVecAlphaVecPow(x(i),0.5*dTym,dev1(i),1);

        // k2
        Array1D<Array1D<double> > dev2(dof);
        RHS_GS(dof, myPCSet, force_mid, x, inpParams,dev2);
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i);
            addVecAlphaVecPow(x(i),0.5*dTym,dev2(i),1);
        }

        // k3
        Array1D<Array1D<double> > dev3(dof);
        RHS_GS(dof, myPCSet, force_mid, x, inpParams,dev3);
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i);
            addVecAlphaVecPow(x(i),dTym,dev3(i),1);
        }

        // k4
        Array1D<Array1D<double> > dev4(dof);
        RHS_GS(dof, myPCSet, force_plus, x, inpParams,dev4);

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
   
void RHS_GS(int dof, PCSet& myPCSet, Array1D<double>& force, Array1D<Array1D<double> >& x, Array1D<double>& inpParams, Array1D<Array1D<double> >& dev){
        for (int i=0;i<dof;i++){
            Array1D<double> temp(x(i));
            dev(i) = temp;
        }
        int nPCTerms = x(0).Length();
        
        if ((abs(inpParams(0))<1e-10)||(abs(inpParams(0)-3)<1e-10)){ //Duffing
            // parse input parameters
            const double zeta = inpParams(1);
            const double epsilon = inpParams(2);
            // buff to store temperary results 
            Array1D<double> temp(nPCTerms,0.e0);                     
            // nonlinear term
            myPCSet.IPow(x(1),temp,3);
            for (int ind=0;ind<nPCTerms;ind++){
                dev(1)(ind) = x(0)(ind);//velocity term
                dev(0)(ind) = force(ind)-temp(ind)*epsilon-2*zeta*x(0)(ind)-x(1)(ind);
            }
        }
        if (abs(inpParams(0)-1)<1e-10){ //Lorenz
            // parse input parameters
            const double a = inpParams(1);
            const double b = inpParams(2);
            const double G = inpParams(3);
            // update x
            Array1D<double> temp(nPCTerms,0.e0);
            addVecAlphaVecPow(temp,a,force,1); //aF
            addVecAlphaVecPow(temp,-a,x(0),1);//aF-ax
            Array1D<double> temp2(nPCTerms,0.e0);
            myPCSet.IPow(x(1),temp2,2); //y^2
            addVecAlphaVecPow(temp,-1,temp2,1);//-y^2-ax+aF
            Array1D<double> temp3(nPCTerms,0.e0);
            myPCSet.IPow(x(2),temp3,2); //z^2
            addVecAlphaVecPow(temp,-1,temp3,1);//-y^2-z^2-ax+aF
            dev(0) = temp;
            // update y
            Array1D<double> tempxy(nPCTerms,0.e0);
            Array1D<double> tempxz(nPCTerms,0.e0);
            myPCSet.Prod(x(0),x(1),tempxy); //xy
            myPCSet.Prod(x(0),x(2),tempxz); //xz
            copy(tempxy,temp);
            copy(tempxz,temp2);
            addVecAlphaVecPow(temp,-b,temp2,1);//xy-bxz
            subtractVec(x(1),temp);//xy-bxz-y
            temp2.Resize(nPCTerms);
            temp2(0) = G;
            addVec(temp2,temp);//xy-bxz-y+G
            dev(1) = temp;
            // update z
            copy(tempxz,temp);
            copy(tempxy,temp2);
            addVecAlphaVecPow(temp,b,temp2,1);//bxy+xz
            subtractVec(x(2),temp);
            dev(2) = temp;
        }
}

Array1D<double> postprocess_GS(int nPCTerms, int nStep,  Array2D<double>& solution, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS){
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
        if (i % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut(i,i*dTym,temp(0),StDv(i));
        }
    }

    Array1D<double> mean;
    getCol(solution,0,mean);
    Array1D<double> e =  error(mean, StDv, mstd_MCS);    
    return e;
}
