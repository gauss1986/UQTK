#include <math.h>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "UtilsDuffing.h"

void GS(int dof, PCSet& myPCSet, int order, int dim, int nPCTerms, string pcType, int nStep, Array1D<double>& initial, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS, Array1D<Array2D<double> >& solution){

    // Working array of the solutions
    Array1D<double> vel_temp(nPCTerms,0.e0);
    Array1D<double> dis_temp(nPCTerms,0.e0);
    Array1D<Array1D<double> > result(dof);  //variable to store results at each iteration step

    // initial step
    vel_temp(0) = initial(0);
    dis_temp(0) = initial(1);
    result(0) = vel_temp;
    result(1) = dis_temp;
    for (int i=0;i<dof;i++){
        Array2D<double> temp(nStep+1,nPCTerms,0.e0);
        solution(i) = temp;
        solution(i).replaceRow(result(i),0);
    }
    
    // Forward run
    Array1D<double> force_current(nPCTerms,0.e0);
    Array1D<double> force_mid(nPCTerms,0.e0);
    Array1D<double> force_plus(nPCTerms,0.e0);
    for (int ix=0;ix<nStep;ix++){
        // Subtract the current force
        getRow(f_GS,2*ix,force_current);
        getRow(f_GS,2*ix+1,force_mid);
        getRow(f_GS,2*ix+2,force_plus);
        // Step forward 
        forward_duffing_GS(dof, myPCSet, inpParams, force_current, force_mid, force_plus, dTym, result);
        // Update dis/vel
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

Array1D<double> postprocess_GS(int nPCTerms, int nStep, double init, Array2D<double>& solution, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS){
    // Output solution (mean and std) 
    Array1D<double> temp(nPCTerms,0.e0);
    temp(0) = init;
    Array1D<double> StDv(nStep+1,0.e0);
    for (int i=0;i<nStep+1;i++){
        getRow(solution,i,temp); 
        double StDv_temp = myPCSet.StDv(temp);
        StDv(i) = StDv_temp; 
        // Write time and solution to file
        WriteModesToFilePtr(i, temp.GetArrayPointer(), nPCTerms, GS_dump);
        // Write dis (mean and std) to file and screen
        WriteMeanStdDevToFilePtr(i,temp(0),StDv_temp,GSstat_dump);
        if (i % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut(i,i*dTym,temp(0),StDv_temp);
        }
    }

    Array1D<double> mean;
    getCol(solution,0,mean);
    Array1D<double> e =  error(mean, StDv, mstd_MCS);    
    return e;
}
