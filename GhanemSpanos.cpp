#include <math.h>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include "GhanemSpanos.h"
#include "UtilsDuffing.h"
#include <ctime>

Array2D<double> GS(PCSet& myPCSet, int order, int dim, int nPCTerms, string pcType, int nStep, double dis0, double vel0, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS){

    // Working array of the solutions
    Array1D<double> dis_temp(nPCTerms,0.e0);
    Array1D<double> vel_temp(nPCTerms,0.e0);
    Array2D<double> dis_GS(nStep+1,nPCTerms,0.e0);
    Array2D<double> vel_GS(nStep+1,nPCTerms,0.e0);
    dis_temp(0) = dis0;
    vel_temp(0) = vel0;
    dis_GS.replaceRow(dis_temp,0);
    vel_GS.replaceRow(vel_temp,0);
    
    // Forward run
    for (int ix=0;ix<nStep;ix++){
        // Subtract the current force
        Array1D<double> force_current(nPCTerms,0.e0);
        Array1D<double> force_mid(nPCTerms,0.e0);
        Array1D<double> force_plus(nPCTerms,0.e0);
        getRow(f_GS,2*ix,force_current);
        getRow(f_GS,2*ix+1,force_mid);
        getRow(f_GS,2*ix+2,force_plus);
        // Step forward 
        forward_duffing_GS(myPCSet, inpParams, force_current, force_mid, force_plus, dTym, dis_temp, vel_temp);
        // Update dis/vel
        dis_GS.replaceRow(dis_temp,ix+1);
        vel_GS.replaceRow(vel_temp,ix+1);
        for (int j=0;j<nPCTerms;j++){
            if ((dis_temp(j) != dis_temp(j) )|| (vel_temp(j)!= vel_temp(j))){
                cout << "GS is producing NaN values at time step " << ix << " on term "<< j<< " !\n"<< endl;
                exit(1);
            }
        }
    }
    
    return (dis_GS);
}

void forward_duffing_GS(PCSet& myPCSet, Array1D<double>& inpParams, Array1D<double>& force_current, Array1D<double>& force_mid, Array1D<double>& force_plus,  double dTym, Array1D<double>& dis, Array1D<double>& vel){
        // Prepare force term for PC calculation
        int nPCTerms = force_current.Length();
        
        // variable to store right hand side
        Array1D<double> dudt1(nPCTerms,0.e0);
        Array1D<double> dudt2(nPCTerms,0.e0);
        Array1D<double> dudt3(nPCTerms,0.e0);
        Array1D<double> dudt4(nPCTerms,0.e0);
        Array1D<double> dvdt1(nPCTerms,0.e0);
        Array1D<double> dvdt2(nPCTerms,0.e0);
        Array1D<double> dvdt3(nPCTerms,0.e0);
        Array1D<double> dvdt4(nPCTerms,0.e0);
        
        // parse input parameters
        const double zeta = inpParams(1);
        const double epsilon = inpParams(2);

        //Save solution at current time step
        Array1D<double> dis0(dis);
        Array1D<double> vel0(vel);
        
        // Integrate with classical 4th order Runge-Kutta
        // k1
        RHS_GS(myPCSet, force_current, dis, vel, dudt1, dvdt1, epsilon, zeta);
        
        // Advance to mid-point
        addVecAlphaVecPow(dis,0.5*dTym,dudt1,1);
        addVecAlphaVecPow(vel,0.5*dTym,dvdt1,1);
    
        // k2
        RHS_GS(myPCSet, force_mid, dis, vel, dudt2, dvdt2, epsilon, zeta);

        // Advance to mid-point
        dis = dis0;
        vel = vel0;
        addVecAlphaVecPow(dis,0.5*dTym,dudt2,1);
        addVecAlphaVecPow(vel,0.5*dTym,dvdt2,1);

        // k3
        RHS_GS(myPCSet, force_mid, dis, vel, dudt3, dvdt3, epsilon, zeta);

        // Advance
        dis = dis0;
        vel = vel0;
        addVecAlphaVecPow(dis,dTym,dudt3,1);
        addVecAlphaVecPow(vel,dTym,dvdt3,1);

        // k4
        RHS_GS(myPCSet, force_plus, dis, vel, dudt4, dvdt4, epsilon, zeta);

        // Advance to next time step
        dis = dis0;
        prodVal(dudt2,2);
        prodVal(dudt3,2);
        addVec(dudt3,dudt4);
        addVec(dudt2,dudt4);
        addVec(dudt1,dudt4);
        addVecAlphaVecPow(dis,dTym/6,dudt4,1);
        
        vel = vel0;
        prodVal(dvdt2,2);
        prodVal(dvdt3,2);
        addVec(dvdt3,dvdt4);
        addVec(dvdt2,dvdt4);
        addVec(dvdt1,dvdt4);
        addVecAlphaVecPow(vel,dTym/6,dvdt4,1);
   }
   
void RHS_GS(PCSet& myPCSet, Array1D<double>& force, Array1D<double>& dis, Array1D<double>& vel, Array1D<double>& dudt, Array1D<double>& dvdt, double epsilon, double zeta){
        int nPCTerms = dis.Length();
        
        // dvdt
        // buff to store temperary results 
        Array1D<double> temp(nPCTerms,0.e0);                     
        
        // nonlinear term
        myPCSet.IPow(dis,temp,3);
        prodVal(temp,epsilon);
        myPCSet.Subtract(force,temp,dvdt);
        
        // damping term
        temp = vel;
        prodVal(temp,2*zeta);
        myPCSet.SubtractInPlace(dvdt,temp);
        
        // stiffness term
        myPCSet.SubtractInPlace(dvdt,dis);
            
        // dudt
        dudt = vel;
}

Array1D<double> postprocess_GS(int nPCTerms, int nStep, double dis0, Array2D<double>& dis_GS, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS){
    // Output solution (mean and std) 
    Array1D<double> dis_temp(nPCTerms,0.e0);
    dis_temp(0) = dis0;
    Array1D<double> disStDv(nStep+1,0.e0);
    double StDv_temp = 0.e0;
    // Initial conditions
    WriteModesToFilePtr(0, dis_temp.GetArrayPointer(), nPCTerms, GS_dump);
    WriteMeanStdDevToFilePtr(0,dis_temp(0),StDv_temp,GSstat_dump);
    WriteMeanStdDevToStdOut(0,0,dis_temp(0),StDv_temp);
    for (int i=1;i<nStep+1;i++){
        getRow(dis_GS,i,dis_temp); 
        StDv_temp = myPCSet.StDv(dis_temp);
        disStDv(i) = StDv_temp; 
        // Write time and solution to file
        WriteModesToFilePtr(i, dis_temp.GetArrayPointer(), nPCTerms, GS_dump);
        // Write dis (mean and std) to file and screen
        WriteMeanStdDevToFilePtr(i,dis_temp(0),StDv_temp,GSstat_dump);
        if (i % ((int) nStep/10) == 0){
            WriteMeanStdDevToStdOut(i,i*dTym,dis_temp(0),StDv_temp);
        }
    }

    Array1D<double> dis_mean;
    getCol(dis_GS,0,dis_mean);
    Array1D<double> e =  error(dis_mean, disStDv, mstd_MCS);    
    return e;
}
