#include <math.h>
#include <cmath>
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

double MCS(int dof, int nspl, int dim, int nStep, int nkl, double dTym, double fbar, Array2D<double>& scaledKLmodes, Array1D<double>& inpParams, Array2D<double>& samPts, Array1D<Array2D<double> >& result, Array1D<double>& initial){
    // Time marching steps
    TickTock tt;
    tt.tick();

    int nthreads;
    #pragma omp parallel default(none) shared(result,dof,nStep,nspl,fbar,samPts,nkl,scaledKLmodes,dTym,inpParams,nthreads,initial) 
    {
    #pragma omp for 
    for (int iq=0;iq<nspl;iq++){
        Array1D<double> totalforce(2*nStep+1,0.e0);
        sample(samPts,iq,2*nStep,fbar,nkl,scaledKLmodes,totalforce,inpParams);
        
        // zero initial condition
        Array1D<double> tempx(initial);
        // initialize the forcing terms
        Array1D<double> tempf(3,0.e0);
        tempf(2) = totalforce(0);
        for (int idof=0;idof<dof;idof++){
            result(idof)(0,iq) = tempx(idof);
        }
        for (int ix=0;ix<nStep;ix++){
            tempf(0) = tempf(2);
            tempf(1) = totalforce(ix*2+1);
            tempf(2) = totalforce(ix*2+2);
            forward_duffing_dt(inpParams, tempf, dTym, tempx);
            for (int idof=0;idof<dof;idof++){
                result(idof)(ix+1,iq) = tempx(idof);
            }
        }
       
    }
    //output the number of threads
    #pragma omp single
    nthreads = omp_get_num_threads();
    }
    cout << "Number of threads in OMP:" << nthreads << endl;
    
    //report the time cost
    tt.tock("Took");
    double t = tt.silent_tock();
    return (t);
}
   
    void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force,  double dTym, Array1D<double>& x){
        // Integrate with classical 4th order Runge-Kutta

        //Save solution at current time step
        Array1D<double> x0(x);
        
        // k1
        Array1D<double> dxdt1 = RHS(force(0),x,inpParams);
        // Advance to mid-point
        addVecAlphaVecPow(x,0.5*dTym,dxdt1,1);

        // k2
        Array1D<double> dxdt2 = RHS(force(1),x,inpParams);
        // Advance to mid-point
        copy(x0,x);
        addVecAlphaVecPow(x,0.5*dTym,dxdt2,1);

        // k3
        Array1D<double> dxdt3 = RHS(force(1),x,inpParams);
        // Advance
        copy(x0,x);
        addVecAlphaVecPow(x,dTym,dxdt3,1);

        // k4
        Array1D<double> dxdt4 = RHS(force(2),x,inpParams);

        // Advance to next time step
        for (unsigned int i=0;i<x.XSize();i++){
            x(i) = x0(i)+dTym/6*(dxdt1(i)+2*dxdt2(i)+2*dxdt3(i)+dxdt4(i));
        }
   }
   
   Array1D<double> RHS(double force, Array1D<double>& x,Array1D<double>& inpParams){
        Array1D<double> dxdt(x);
        if (abs(inpParams(0))<1e-10){ //Duffing
            // parse input parameters
            const double zeta = inpParams(1);
            const double epsilon = inpParams(2);

            double temp = pow(x(1),3);
            if (temp != temp){
                cout << "Pow 3 in MCS is generating NaN!\n" << endl;
                exit(1);
            }
            dxdt(0) = force-epsilon*temp-2*zeta*x(0)-x(1);
            dxdt(1) = x(0);
            } 
        if (abs(inpParams(0)-1)<1e-10){ //Lorenz
            // parse input parameters
            const double a = inpParams(1);
            const double b = inpParams(2);
            const double G = inpParams(3);

            dxdt(0)=-pow(x(1),2)-pow(x(2),2)-a*x(0)+a*force;
            dxdt(1)=x(0)*x(1)-b*x(0)*x(2)-x(1)+G;
            dxdt(2)=b*x(0)*x(1)+x(0)*x(2)-x(2);
            }
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

void sample(Array2D<double>& samPts,int iq,int nStep,double fbar,int nkl,Array2D<double>& scaledKLmodes,Array1D<double>& totalforce, Array1D<double>& inpParams){    
    for (int it=0;it<nStep+1;it++){
        if (abs(inpParams(0))<1e-10){ //Duffing
            totalforce(it) = fbar;
        }
        if (abs(inpParams(0)-1)<1e-10){ //Lorenz
            totalforce(it) = fbar+inpParams(4)*cos(inpParams(5)*it);
	    }
        for (int iy=0;iy<nkl;iy++){
            totalforce(it) = totalforce(it)+samPts(iq,iy)*scaledKLmodes(it,iy);
        }
    }
}

