#include <math.h>
#include <omp.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "MCS.h"
#include "nGhanemSpanos.h"
#include "Utils.h"

void nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<Array2D<double> >& f_GS, Array1D<Array2D<double> >& solution){
    int nPCTerms = myPCSet.GetNumberPCTerms();

    // Initialize
    Array1D<Array1D<double> > result(dof);
    for (int i=0;i<dof;i++){
        // save to Solution
        Array2D<double> temp(nStep+1,2*nPCTerms,0.e0);
        solution(i) = temp;
        solution(i).replaceRow(initial(i),0);
        // result is working variable
        Array1D<double> temp2(initial(i));
        result(i)=temp2;
    }
    
    // Forward run
    Array1D<Array1D<double> > force_current(dof);
    Array1D<Array1D<double> > force_mid(dof);
    Array1D<Array1D<double> > force_plus(dof);
    cout << "Starting time marching step..." << endl;
    for (int ix=0;ix<nStep;ix++){
        // Subtract the force at currect, mid and next step
        for (int i=0;i<dof;i++){
            Array1D<double> temp_force(nPCTerms,0.e0);
            getRow(f_GS(i),2*ix,temp_force);
            force_current(i)=temp_force;
            Array1D<double> temp_force1(nPCTerms,0.e0);
            getRow(f_GS(i),2*ix+1,temp_force1);
            force_mid(i)=temp_force1;
            Array1D<double> temp_force2(nPCTerms,0.e0);
            getRow(f_GS(i),2*ix+2,temp_force2);
            force_plus(i)=temp_force1;
        }
        // Step forward 
        forward_duffing_nGS(dof, myPCSet, epsilon, mck, force_current, force_mid, force_plus, dTym, result);
        // Update solution
        for (int i=0;i<dof;i++){
            solution(i).replaceRow(result(i),ix+1);
        }
    }   

    return; 
}

void forward_duffing_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck, Array1D<Array1D<double> >& force_current, Array1D<Array1D<double> >& force_mid, Array1D<Array1D<double> >& force_plus,  double dTym, Array1D<Array1D<double> >& x){
        int nPCTerms = myPCSet.GetNumberPCTerms();

        Array1D<Array1D<double> > u(dof);
        Array1D<Array1D<double> > v(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> temp_u(nPCTerms,0.e0);
            Array1D<double> temp_v(nPCTerms,0.e0);
            u(i) = temp_u;
            v(i) = temp_v;
            for (int j=0;j<nPCTerms;j++){
                v(i)(j)=x(i)(j);
                u(i)(j)=x(i)(j+nPCTerms);
            }
        }

        //cout << "Save solution at current time step..." << endl;
        //Save solution at current time step
        Array1D<Array1D<double> > u0(dof);
        Array1D<Array1D<double> > v0(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> temp(u(i));
            u0(i)=temp;
            Array1D<double> temp2(v(i));
            v0(i)=temp2;
        }
        
        // Integrate with classical 4th order Runge-Kutta
        // k1
        Array1D<Array1D<double> > du1(dof);
        Array1D<Array1D<double> > dv1(dof);
        Array1D<Array1D<double> > kbar1 = kbar(dof, myPCSet, u, epsilon, mck);
        dev_nGS(dof, myPCSet, force_current, mck, u, v, du1, dv1, kbar1);
        for (int i=0;i<dof;i++){
            addVecAlphaVecPow(u(i),0.5*dTym,du1(i),1);
            addVecAlphaVecPow(v(i),0.5*dTym,dv1(i),1);
        }

        // k2
        Array1D<Array1D<double> > du2(dof);
        Array1D<Array1D<double> > dv2(dof);
        Array1D<Array1D<double> > kbar2 = kbar(dof, myPCSet, u, epsilon, mck);
        dev_nGS(dof, myPCSet, force_mid, mck, u, v, du2, dv2, kbar2);
        for (int i=0;i<dof;i++){
            u(i) = u0(i);
            v(i) = v0(i);
            addVecAlphaVecPow(u(i),0.5*dTym,du2(i),1);
            addVecAlphaVecPow(v(i),0.5*dTym,dv2(i),1);
        }

        // k3
        Array1D<Array1D<double> > du3(dof);
        Array1D<Array1D<double> > dv3(dof);
        Array1D<Array1D<double> > kbar3 = kbar(dof, myPCSet, u, epsilon, mck);
        dev_nGS(dof, myPCSet, force_mid, mck, u, v, du3, dv3, kbar3);
        for (int i=0;i<dof;i++){
            u(i) = u0(i);
            v(i) = v0(i);
            addVecAlphaVecPow(u(i),dTym,du3(i),1);
            addVecAlphaVecPow(v(i),dTym,dv3(i),1);
        }

        // k4
        Array1D<Array1D<double> > du4(dof);
        Array1D<Array1D<double> > dv4(dof);
        Array1D<Array1D<double> > kbar4 = kbar(dof, myPCSet, u, epsilon, mck);
        dev_nGS(dof, myPCSet, force_plus, mck, u, v, du4, dv4, kbar4);

        // Advance to next time step
        for (int i=0;i<dof;i++){
            //cout << "dof " << i << endl;
            //u
            u(i) = u0(i);
            prodVal(du2(i),2);
            prodVal(du3(i),2);
            addVec(du3(i),du4(i));
            addVec(du2(i),du4(i));
            addVec(du1(i),du4(i));
            addVecAlphaVecPow(u(i),dTym/6,du4(i),1);
            //v
            v(i) = v0(i);
            prodVal(dv2(i),2);
            prodVal(dv3(i),2);
            addVec(dv3(i),dv4(i));
            addVec(dv2(i),dv4(i));
            addVec(dv1(i),dv4(i));
            addVecAlphaVecPow(v(i),dTym/6,dv4(i),1);
            //x
            for (int j=0;j<nPCTerms;j++){
                x(i)(j)=v(i)(j);
                x(i)(j+nPCTerms)=u(i)(j);    
            }
        }     
           
        return; 
   }
   
Array1D<Array1D<double> > kbar(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck){
    // kbar 
    int nPC = myPCSet.GetNumberPCTerms();
    Array1D<Array1D<double> > kbar(dof);
    // kbar0
    Array1D<double> temp_u(nPC,0.e0);
    Array1D<double> temp_k(nPC,0.e0);
    myPCSet.IPow(u(0),temp_u,2);
    myPCSet.Prod(temp_u,epsilon(0),temp_k);
    temp_k(0)=temp_k(0)+1; 
    prodVal(temp_k,mck(2)(0));
    kbar(0) = temp_k;
    // kbar rest 
    for (int i=1;i<dof;i++){
        Array1D<double> temp_k1(nPC,0.e0);
        Array1D<double> temp_u1(u(i));
        subtractVec(u(i-1),temp_u1); 
        myPCSet.IPow(temp_u1,temp_u,2);
        myPCSet.Prod(temp_u,epsilon(i),temp_k1);
        temp_k1(0)=temp_k1(0)+1; 
        prodVal(temp_k1,mck(2)(i));
        kbar(i) = temp_k1;
    }

    return(kbar);
}

void dev_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& force, Array1D<Array1D<double> >& mck, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& v, Array1D<Array1D<double> >& du, Array1D<Array1D<double> >& dv, Array1D<Array1D<double> >& kbar){
    int nPCTerms = myPCSet.GetNumberPCTerms();
    // du = v
    for (int i=0;i<dof;i++){
        Array1D<double> temp(v(i));
        du(i) = temp;
    }
    //int nPCTerms = x(0).Length();
       
    // compute kbaru and kbaru_plus 
    Array1D<Array1D<double> > kbaru(dof);
    Array1D<Array1D<double> > kbaru_plus(dof-1);
    Array1D<double> temp_k(nPCTerms,0.e0);
    for (int i=0;i<dof-1;i++){
        //Array1D<double> temp_k(nPCTerms,0.e0);
        kbaru(i)=temp_k; 
        kbaru_plus(i)=temp_k; 
        myPCSet.Prod(kbar(i),u(i),kbaru(i));   
        myPCSet.Prod(kbar(i+1),u(i),kbaru_plus(i));   
    }
    kbaru(dof-1) = temp_k;
    myPCSet.Prod(kbar(dof-1),u(dof-1),kbaru(dof-1));   

    // compute cv and cv_plus 
    Array1D<Array1D<double> > cv(dof);
    Array1D<Array1D<double> > cv_plus(dof-1);
    for (int i=0;i<dof-1;i++){
        cv(i) = v(i);
        cv_plus(i) = v(i);
        prodVal(cv(i),mck(1)(i));
        prodVal(cv_plus(i),mck(1)(i+1));
    }
    cv(dof-1) = v(dof-1);
    prodVal(cv(dof-1),mck(1)(dof-1));

    // compute the acceleration
    // dof 0
    Array1D<double> temp_d(force(0));
    subtractVec(kbaru(0),temp_d);
    subtractVec(kbaru_plus(0),temp_d);
    addVec(kbaru(1),temp_d);
    subtractVec(cv(0),temp_d);
    subtractVec(cv_plus(0),temp_d);
    addVec(cv(1),temp_d);
    prodVal(temp_d,1/mck(0)(0));
    dv(0)=temp_d;
    // dof dof
    Array1D<double> temp_d1(force(dof-1));
    subtractVec(kbaru(dof-1),temp_d1);
    addVec(kbaru_plus(dof-2),temp_d1);
    subtractVec(cv(dof-1),temp_d1);
    addVec(cv_plus(dof-2),temp_d1);
    prodVal(temp_d1,1/mck(0)(dof-1));
    dv(dof-1)=temp_d1;
    // dof n
    for (int i=1;i<dof-1;i++){
        Array1D<double> temp_d2(force(i));
        subtractVec(kbaru(i),temp_d2);
        subtractVec(kbaru_plus(i),temp_d2);
        addVec(kbaru_plus(i-1),temp_d2);
        addVec(kbaru(i+1),temp_d2);
        subtractVec(cv(i),temp_d2);
        subtractVec(cv_plus(i),temp_d2);
        addVec(cv_plus(i-1),temp_d2);
        addVec(cv(i+1),temp_d2);
        prodVal(temp_d2,1/mck(0)(i));
        dv(i)=temp_d2;
    }
    
    return;
}

Array2D<double> postprocess_nGS(int dof, int nStep,  Array1D<Array2D<double> >& solution, PCSet& myPCSet, double dTym, int ord,Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array2D<double>& e2){
    int nPCTerms = myPCSet.GetNumberPCTerms();
    // Open files to write out solutions
    ostringstream s1;
    s1 << "nGS_dis_m_" << ord<<".dat";
    string dis_m(s1.str());
    ostringstream s2;
    s2 << "nGS_dis_s_" << ord<<".dat";
    string dis_s(s2.str());
    ostringstream s3;
    s3 << "nGS_vel_m_" << ord<<".dat";
    string vel_m(s3.str());
    ostringstream s4;
    s4 << "nGS_vel_s_" << ord<<".dat";
    string vel_s(s4.str());

    // Output solution (mean and std) 
    Array1D<double> temp_u(nPCTerms,0.e0);
    Array1D<double> temp_v(nPCTerms,0.e0);
    Array2D<double> StDv_u(nStep+1,dof,0.e0);
    Array2D<double> StDv_v(nStep+1,dof,0.e0);
    Array2D<double> mean_u(nStep+1,dof,0.e0);
    Array2D<double> mean_v(nStep+1,dof,0.e0);
    for (int j=0;j<dof;j++){
        for (int i=0;i<nStep+1;i++){
            for (int k=0;k<nPCTerms;k++){
                temp_v(k)=solution(j)(i,k);
                temp_u(k)=solution(j)(i,nPCTerms+k);
            }
            StDv_u(i,j) = myPCSet.StDv(temp_u);
            StDv_v(i,j) = myPCSet.StDv(temp_v);
        }
        Array1D<double> temp2(nStep+1,0.e0);
        getCol(solution(j),0,temp2);
        mean_v.replaceCol(temp2,j);
        getCol(solution(j),nPCTerms,temp2);
        mean_u.replaceCol(temp2,j);
    }
    
    write_datafile(mean_u,dis_m.c_str());
    write_datafile(StDv_u,dis_s.c_str());
    write_datafile(mean_v,vel_m.c_str());
    write_datafile(StDv_v,vel_s.c_str());

    Array1D<Array2D<double> > et(dof);
    Array2D<double> e =  nerror(ord, dof, nStep, et,mean_u,StDv_u,mean_v,StDv_v,mean_MCS,std_MCS, e2);

    return(e);
}
    
Array2D<double> nerror(int ord, int dof, int nStep, Array1D<Array2D<double> >& et,Array2D<double>& mean_u,Array2D<double>& StDv_u, Array2D<double>& mean_v, Array2D<double>& StDv_v, Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array2D<double>& e2){
    Array2D<double> e(dof,4,0.e0);
    Array2D<double> temp(dof,4,0.e0);

    for (int i=0;i<dof;i++){
        Array2D<double> error(nStep+1,4,0.e0);
        for (int ix=1;ix<nStep+1;ix++){
            error(ix,0)=abs(mean_v(ix,i)-mean_MCS(ix,i))/abs(mean_MCS(ix,i));
            error(ix,1)=abs(StDv_v(ix,i)-std_MCS(ix,i))/abs(std_MCS(ix,i));
            error(ix,2)=abs(mean_u(ix,i)-mean_MCS(ix,i+dof))/abs(mean_MCS(ix,i+dof));
            error(ix,3)=abs(StDv_u(ix,i)-std_MCS(ix,i+dof))/abs(std_MCS(ix,i+dof));
            e(i,0)+=error(ix,0);
            e(i,1)+=error(ix,1);
            e(i,2)+=error(ix,2);
            e(i,3)+=error(ix,3);
            e2(i,0)+=abs(mean_v(ix,i)-mean_MCS(ix,i));
            e2(i,1)+=abs(StDv_v(ix,i)-std_MCS(ix,i));
            e2(i,2)+=abs(mean_u(ix,i)-mean_MCS(ix,i+dof));
            e2(i,3)+=abs(StDv_u(ix,i)-std_MCS(ix,i+dof));
            temp(i,0)+=abs(mean_MCS(ix,i));
            temp(i,1)+=abs(std_MCS(ix,i));
            temp(i,2)+=abs(mean_MCS(ix,i+dof));
            temp(i,3)+=abs(std_MCS(ix,i+dof));
        }
        for (int id=0;id<4;id++){
            e(i,id)=e(i,id)/(nStep);
            e2(i,id)=e2(i,id)/temp(i,id);
        }
        et(i)=error;
        ostringstream name;
        name << "nerror_dof" << i << "_" << ord << ".dat";
        string name_str = name.str();
        write_datafile(error,name_str.c_str());
    }
    
    ostringstream name2;
    name2 << "e2_" << ord << ".dat";
    string name_str2 = name2.str();
    write_datafile(e2,name_str2.c_str());
    
    return(e);            
}
