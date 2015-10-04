#include <math.h>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "Utils.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include <ctime>
#include <omp.h>
#include "ticktock.h"

void MCS(int nspl, int dim, int nStep, int nkl, double dTym, double fbar, Array2D<double>& scaledKLmodes, Array1D<double>& inpParams, Array2D<double>& samPts, Array2D<double>& dis_MC){
    // Time marching steps
    TickTock tt;
    tt.tick();

    int nthreads;
    #pragma omp parallel default(none) shared(dis_MC,nStep,nspl,fbar,dim,samPts,nkl,scaledKLmodes,dTym,inpParams,nthreads) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        // sample force
        Array1D<double> temp(nStep+1,fbar);
        Array1D<double> sampt(dim,0.e0);
        getRow(samPts,iq,sampt);
        for (int iy=0;iy<nkl;iy++){
            Array1D<double> tempf(nStep+1,0.e0);
            getCol(scaledKLmodes,iy,tempf);
            prodVal(tempf,sampt(iy));
            addVec(tempf,temp);
            }
        Array2D<double> totalforce(2,nStep+1,0.e0);
    	totalforce.replaceRow(temp,0);
        
        // initialize the solution
        Array2D<double> x(2,nStep+1,0.e0);
        // zero initial condition
        Array1D<double> tempx(2,0.e0);

        // Time marching steps
        for (int ix=0;ix<nStep;ix++){
            //cout << "Time step No." << ix << endl;
            Array1D<double> tempf(2,0.e0);
            getCol(totalforce,ix,tempf);
            forward_duffing_dt(inpParams, tempf, dTym, tempx);
            x.replaceCol(tempx,ix+1);
        }
       
        // store solution for this sample force
        Array1D<double> dis(nStep+1,0.e0);
        getRow(x,1,dis);
        for (int ix=0;ix<nStep+1;ix++)
            dis_MC(ix,iq) = dis(ix);
    }
    #pragma omp single
    nthreads = omp_get_num_threads();
    }
    cout << "Number of threads in OMP:" << nthreads << endl;
        
    tt.tock("Took");
    return;
}
   
    void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force, double dTym, Array1D<double>& x){
        // variable to store right hand side
        Array1D<double> dxdt1(2,0.e0);    
        Array1D<double> dxdt2(2,0.e0);    
        Array1D<double> dxdt3(2,0.e0);    
        Array1D<double> dxdt4(2,0.e0);    
        
        // parse input parameters
        const double zeta = inpParams(0);
        const double epsilon = inpParams(1);

        //Save solution at current time step
        Array1D<double> x0(2,0.e0);
        x0(0) = x(0);
        x0(1) = x(1);
        
        // Integrate with classical 4th order Runge-Kutta
        // k1
        RHS(force,x,dxdt1,epsilon,zeta);
        
        // Advance to mid-point
        x(0) = x0(0) + 0.5*dTym*dxdt1(0);
        x(1) = x0(1) + 0.5*dTym*dxdt1(1);

        // k2
        RHS(force,x,dxdt2,epsilon,zeta);

        // Advance to mid-point
        x(0) = x0(0) + 0.5*dTym*dxdt2(0);
        x(1) = x0(1) + 0.5*dTym*dxdt2(1);

        // k3
        RHS(force,x,dxdt3,epsilon,zeta);

        // Advance
        x(0) = x0(0) + dTym*dxdt3(0);
        x(1) = x0(1) + dTym*dxdt3(1);

        // k4
        RHS(force,x,dxdt4,epsilon,zeta);

        // Advance to next time step
        x(0) = x0(0) + dTym*(dxdt1(0)+2*dxdt2(0)+2*dxdt3(0)+dxdt4(0))/6;
        x(1) = x0(1) + dTym*(dxdt1(1)+2*dxdt2(1)+2*dxdt3(1)+dxdt4(1))/6;
   }
   
   void RHS(Array1D<double>& force, Array1D<double>& x, Array1D<double>& dxdt, double epsilon, double zeta){
        double temp = pow(x(1),3);
        if (temp != temp){
            cout << "Pow 3 in MCS is generating NaN!\n" << endl;
            exit(1);
        }
        dxdt(0) = force(0)-epsilon*temp-2*zeta*x(0)-x(1);
        dxdt(1) = x(0); 
   }
   
   Array1D<double> mStd(Array1D<double>& x,int nspl){
        Array1D<double> mstd(2,0.e0);
        double sum = 0.e0;
        double sum2 = 0.e0;
        for (int i=0;i<nspl;i++){
            sum += x(i);
            sum2 += pow(x(i),2);
        }
        mstd(0) = sum/nspl;
        mstd(1) = sqrt((sum2-nspl*pow(mstd(0),2))/nspl); 
        return(mstd);   
   }
