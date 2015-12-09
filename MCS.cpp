#include <math.h>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "UtilsDuffing.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include <ctime>
#include <omp.h>
#include "ticktock.h"

double MCS(int nspl, int dim, int nStep, int nkl, double dTym, double fbar, Array2D<double>& scaledKLmodes, Array1D<double>& inpParams, Array2D<double>& samPts, Array2D<double>& dis_MC){
    // Time marching steps
    TickTock tt;
    tt.tick();

    int nthreads;
    #pragma omp parallel default(none) shared(dis_MC,nStep,nspl,fbar,dim,samPts,nkl,scaledKLmodes,dTym,inpParams,nthreads) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array2D<double> totalforce(2,2*nStep+1,0.e0);
        sample_duffing(samPts,iq,2*nStep,fbar,nkl,scaledKLmodes,totalforce);
        
        // initialize the solution
        Array2D<double> x(2,nStep+1,0.e0);
        // zero initial condition
        Array1D<double> tempx(2,0.e0);

        // Time marching steps
        for (int ix=0;ix<nStep;ix++){
            //cout << "Time step No." << ix << endl;
            Array1D<double> tempf1(2,0.e0);
            Array1D<double> tempf2(2,0.e0);
            Array1D<double> tempf3(2,0.e0);
            getCol(totalforce,ix*2,tempf1);
            getCol(totalforce,ix*2+1,tempf2);
            getCol(totalforce,ix*2+2,tempf3);
            forward_duffing_dt(inpParams, tempf1, tempf2, tempf3, dTym, tempx);
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
    double t = tt.silent_tock();
    return (t);
}
   
    void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force1, Array1D<double>& force2, Array1D<double>& force3, double dTym, Array1D<double>& x){
        // Integrate with classical 4th order Runge-Kutta

        //Save solution at current time step
        Array1D<double> x0(x);
        
        // k1
        Array1D<double> dxdt1 = RHS(force1,x,inpParams);
        // Advance to mid-point
        addVecAlphaVecPow(x,0.5*dTym,dxdt1,1);

        // k2
        Array1D<double> dxdt2 = RHS(force2,x,inpParams);
        // Advance to mid-point
        copy(x0,x);
        addVecAlphaVecPow(x,0.5*dTym,dxdt2,1);

        // k3
        Array1D<double> dxdt3 = RHS(force2,x,inpParams);
        // Advance
        copy(x0,x);
        addVecAlphaVecPow(x,dTym,dxdt3,1);

        // k4
        Array1D<double> dxdt4 = RHS(force3,x,inpParams);

        // Advance to next time step
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i)+dTym/6*(dxdt1(i)+2*dxdt2(i)+2*dxdt3(i)+dxdt4(i));
        }
   }
   
   Array1D<double> RHS(Array1D<double>& force, Array1D<double>& x,Array1D<double>& inpParams){
        Array1D<double> dxdt(x);
        // parse input parameters
        const double zeta = inpParams(0);
        const double epsilon = inpParams(1);

        double temp = pow(x(1),3);
        if (temp != temp){
            cout << "Pow 3 in MCS is generating NaN!\n" << endl;
            exit(1);
        }
        dxdt(0) = force(0)-epsilon*temp-2*zeta*x(0)-x(1);
        dxdt(1) = x(0); 

        return dxdt;
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

   Array1D<double> error(Array1D<double>& dis, Array1D<double>& StDv, Array2D<double>& mstd_MCS){
   // Return the integrated error in mean/std
        Array1D<double> e(2,0.e0);
        for (unsigned int i=0;i<dis.XSize();i++){
            e(0) = e(0) + fabs(dis(i) - mstd_MCS(0,i));
            e(1) = e(1) + fabs(StDv(i) - mstd_MCS(1,i));  
        }
        Array1D<double> temp_m;
        getRow(mstd_MCS,0,temp_m);
        Array1D<double> temp_s;
        getRow(mstd_MCS,1,temp_s);

        // normalized relative error
        e(0) = e(0)/sum(temp_m)*100;
        e(1) = e(1)/sum(temp_s)*100;
    
        return e;
    }

void sample_duffing(Array2D<double>& samPts,int iq,int nStep,double fbar,int nkl,Array2D<double>& scaledKLmodes,Array2D<double>& totalforce){    
    for (int it=0;it<nStep+1;it++){
        totalforce(0,it) = fbar;
	    for (int iy=0;iy<nkl;iy++){
            totalforce(0,it) = totalforce(0,it)+samPts(iq,iy)*scaledKLmodes(it,iy);
        }
    }
}

void sample_lorenz(Array2D<double>& samPts,int iq,int nStep,double fbar,double Amp,double w,int nkl,Array2D<double>& scaledKLmodes,Array2D<double>& totalforce){    
    for (int it=0;it<nStep+1;it++){
        totalforce(0,it) = fbar+Amp*cos(w*it);
	    for (int iy=0;iy<nkl;iy++){
            totalforce(0,it) = totalforce(0,it)+samPts(iq,iy)*scaledKLmodes(it,iy);
        }
    }
}
