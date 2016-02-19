#include <math.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "Utils.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "nMCS.h"
#include <ctime>
#include "ticktock.h"

Array1D<Array1D<double> > ndet(int dof, int nStep, double dTym, Array2D<double>& sampleforce, Array1D<double>& epsilon, Array1D<Array1D<double> >& mck, Array1D<double>& initial){

    // initialize result where solution from each step is stored
    Array1D<Array1D<double> > result(2*dof);
    Array1D<double> tempu(dof,0.e0);
    Array1D<double> tempv(dof,0.e0);
    for (int i=0;i<dof;i++){
        Array1D<double> temp_result(nStep+1,0.e0);
        result(i) = temp_result;
        result(i)(0)=initial(i);
        result(dof+i) = temp_result;
        result(dof+i)(0)=initial(dof+i);
        tempv(i)=initial(i);
        tempu(i)=initial(dof+i);
    }

    // time marching steps
    // initialize the forcing terms
    Array2D<double> force(3,dof,0.e0);
    Array1D<double> tempf(dof,0.e0);
    getRow(sampleforce,0,tempf);
    for (int ix=0;ix<nStep;ix++){
        Array1D<double> tempf2(dof,0.e0);

        force.replaceRow(tempf,0);
        getRow(sampleforce,ix*2+1,tempf2);
        force.replaceRow(tempf2,1);
        getRow(sampleforce,ix*2+2,tempf);
        force.replaceRow(tempf,2);

        nforward_duffing_dt(epsilon,dof,force,mck,dTym,tempu,tempv);
        for (int i=0;i<dof;i++){
            result(i)(ix+1)=tempv(i);
            result(i+dof)(ix+1)=tempu(i);
        }
    }

    return result;
}
   
void nforward_duffing_dt(Array1D<double>& epsilon, int dof, Array2D<double>& force, Array1D<Array1D<double> >& mck, double dTym, Array1D<double>& u, Array1D<double>& v){
    // Integrate with classical 4th order Runge-Kutta
    //Save solution at current time step
    Array1D<double> u0(u);
    Array1D<double> v0(v);
        
    // k1
    Array1D<double> temp_f(dof,0.e0);
    getRow(force,0,temp_f);
    Array1D<double> a1(dof,0.e0);
    nRHS(a1,dof,epsilon,mck,temp_f,u,v);
    Array1D<double> v1(v);
    // Advance to mid-point
    addVecAlphaVecPow(u,0.5*dTym,v,1);
    addVecAlphaVecPow(v,0.5*dTym,a1,1);

    // k2
    getRow(force,1,temp_f);
    Array1D<double> a2(dof,0.e0);
    nRHS(a2,dof,epsilon,mck,temp_f,u,v);
    Array1D<double> v2(v);
    // Advance to mid-point
    copy(u0,u);
    addVecAlphaVecPow(u,0.5*dTym,v,1);
    copy(v0,v);
    addVecAlphaVecPow(v,0.5*dTym,a2,1);

    // k3
    //getRow(force,1,temp_f);
    Array1D<double> a3(dof,0.e0);
    nRHS(a3,dof,epsilon,mck,temp_f,u,v);
    Array1D<double> v3(v);
    // Advance
    copy(u0,u);
    addVecAlphaVecPow(u,dTym,v,1);
    copy(v0,v);
    addVecAlphaVecPow(v,dTym,a3,1);

    // k4
    getRow(force,2,temp_f);
    Array1D<double> a4(dof,0.e0);
    nRHS(a4,dof,epsilon,mck,temp_f,u,v);
    Array1D<double> v4(v);

    // Advance to next time step
    for (int i=0;i<dof;i++){
        u(i) = u0(i)+dTym/6*(v1(i)+2*v2(i)+2*v3(i)+v4(i));
        v(i) = v0(i)+dTym/6*(a1(i)+2*a2(i)+2*a3(i)+a4(i));
    }
}
   
void nRHS(Array1D<double>& acc, int dof, Array1D<double>& epsilon, Array1D<Array1D<double> >& mck, Array1D<double>& force, Array1D<double>& u, Array1D<double>& v){
    Array1D<double> kbar(dof,0.e0);
    kbar(0)=mck(2)(0)*(1+epsilon(0)*u(0)*u(0));
    for (int i=1;i<dof;i++){
        kbar(i)=mck(2)(i)*(1+epsilon(i)*pow((u(i)-u(i-1)),2));
    }

    Array1D<double> kbaru(dof,0.e0);
    Array1D<double> kbaru_plus(dof-1,0.e0);
    Array1D<double> cv(dof,0.e0);
    Array1D<double> cv_plus(dof-1,0.e0);
    kbaru(0) = kbar(0)*u(0);
    cv(0) = mck(1)(0)*v(0);
    for (int i=1;i<dof;i++){
        kbaru(i) = kbar(i)*u(i);
        kbaru_plus(i-1) = kbar(i)*u(i-1);
        cv(i) = mck(1)(i)*v(i);
        cv_plus(i-1) = mck(1)(i)*v(i-1);
    }

    acc(0)=(force(0)-kbaru(0)-kbaru_plus(0)+kbaru(1)-cv(0)-cv_plus(0)+cv(1))/mck(0)(0);
    for (int i=1;i<dof-1;i++){
        acc(i)=(force(i)+kbaru_plus(i-1)-kbaru(i)-kbaru_plus(i)+kbaru(i+1)+cv_plus(i-1)-cv(i)-cv_plus(i)+cv(i+1))/mck(0)(i);
    }
    acc(dof-1)=(force(dof-1)-kbaru(dof-1)+kbaru_plus(dof-2)-cv(dof-1)+cv_plus(dof-2))/mck(0)(dof-1);

    return;
}

Array1D<double> error(Array2D<double>& et, Array1D<double>& mean, Array1D<double>& StDv, Array2D<double>& mstd_MCS){
   // Return the integrated error in mean/std
        Array1D<double> e(2,0.e0);
        //Array2D<double> et(mean.XSize(),2,0.e0);
        for (unsigned int i=0;i<mean.XSize();i++){
            et(i,0) = fabs(mean(i) - mstd_MCS(0,i));
            et(i,1) = fabs(StDv(i) - mstd_MCS(1,i));  
        }
        Array1D<double> temp_m;
        getRow(mstd_MCS,0,temp_m);
        Array1D<double> temp_m2;
        getCol(et,0,temp_m2);
        Array1D<double> temp_s;
        getRow(mstd_MCS,1,temp_s);
        Array1D<double> temp_s2;
        getCol(et,1,temp_s2);

        // normalized relative error
        e(0) = sum(temp_m2)/sum(temp_m)*100;
        e(1) = sum(temp_s2)/sum(temp_s)*100;
    
        for (unsigned int i=0;i<mean.XSize();i++){
            et(i,0) = et(i,0)/fabs(mstd_MCS(0,i))*100;
            et(i,1) = et(i,1)/fabs(mstd_MCS(1,i))*100;  
        }

        return e;
}

Array2D<double>  nsample_force(int dof, Array2D<double>& samPts,int iq, Array1D<double>& fbar,int nkl,Array2D<double>& scaledKLmodes){    
    // force is applied to dof 0 only 
    unsigned int nStep = fbar.XSize();
    Array2D<double> totalforce(nStep+1,dof,0.e0);
    for (unsigned int it=0;it<nStep;it++){
        totalforce(it,0) = fbar(it);
        for (int iy=0;iy<nkl;iy++){
            totalforce(it,0) = totalforce(it,0)+samPts(iq,iy)*scaledKLmodes(it,iy);
        }
    }
    return (totalforce);
}

void mstd_MCS(Array1D<double>& result, double& mean, double& std){
    unsigned int nspl = result.XSize();
    mean = 0;
    double var =0;
    for (unsigned int i=0;i<nspl;i++){
        mean += result(i);
        var += result(i)*result(i);
    }
    mean = mean/nspl;
    var = var/nspl-mean*mean;
    std = sqrt(var);
}
