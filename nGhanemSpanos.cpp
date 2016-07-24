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

void nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<Array2D<double> >& f_GS, Array1D<Array2D<double> >& solution){
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
            force_plus(i)=temp_force2;
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

void forward_duffing_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck, Array1D<Array1D<double> >& force_current, Array1D<Array1D<double> >& force_mid, Array1D<Array1D<double> >& force_plus,  double dTym, Array1D<Array1D<double> >& x){
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
   
Array1D<Array1D<double> > kbar(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck){
    // kbar 
    int nPC = myPCSet.GetNumberPCTerms();
    Array1D<Array1D<double> > kbar(dof);
    // kbar0
    Array1D<double> temp_u(nPC,0.e0);
    Array1D<double> temp_k(nPC,0.e0);
    myPCSet.IPow(u(0),temp_u,2);
    myPCSet.Prod(temp_u,epsilon(0),temp_k);
    temp_k(0)=temp_k(0)+1; 
    //prodVal(temp_k,mck(2)(0));
    Array1D<double> temp_k2(nPC,0.e0);
    myPCSet.Prod(temp_k,mck(2)(0),temp_k2);
    kbar(0) = temp_k2;
    // kbar rest 
    for (int i=1;i<dof;i++){
        Array1D<double> temp_k1(nPC,0.e0);
        Array1D<double> temp_u1(u(i));
        subtractVec(u(i-1),temp_u1); 
        myPCSet.IPow(temp_u1,temp_u,2);
        myPCSet.Prod(temp_u,epsilon(i),temp_k1);
        temp_k1(0)=temp_k1(0)+1; 
        Array1D<double> temp_k3(nPC,0.e0);
        //prodVal(temp_k1,mck(2)(i));
        myPCSet.Prod(temp_k1,mck(2)(i),temp_k3);
        kbar(i) = temp_k3;
    }

    return(kbar);
}

void dev_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& force, Array1D<Array1D<Array1D<double> > >& mck, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& v, Array1D<Array1D<double> >& du, Array1D<Array1D<double> >& dv, Array1D<Array1D<double> >& kbar){
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
        //cv(i) = v(i);
        //cv_plus(i) = v(i);
        //prodVal(cv(i),mck(1)(i));
        //prodVal(cv_plus(i),mck(1)(i+1));
        Array1D<double> cv_temp(nPCTerms,0.e0);
        cv(i) = cv_temp;
        myPCSet.Prod(v(i),mck(1)(i),cv(i));
        Array1D<double> cv_temp2(nPCTerms,0.e0);
        cv_plus(i) = cv_temp2;
        myPCSet.Prod(v(i),mck(1)(i+1),cv_plus(i));
    }
    //cv(dof-1) = v(dof-1);
    //prodVal(cv(dof-1),mck(1)(dof-1));
    Array1D<double> cv_temp3(nPCTerms,0.e0);
    cv(dof-1) = cv_temp3;
    myPCSet.Prod(v(dof-1),mck(1)(dof-1),cv(dof-1));

    // compute the acceleration
    // dof 0
    Array1D<double> temp_d(force(0));
    subtractVec(kbaru(0),temp_d);
    subtractVec(kbaru_plus(0),temp_d);
    addVec(kbaru(1),temp_d);
    subtractVec(cv(0),temp_d);
    subtractVec(cv_plus(0),temp_d);
    addVec(cv(1),temp_d);
    //prodVal(temp_d,1/mck(0)(0));
    Array1D<double> temp_d3(nPCTerms,0.e0); 
    myPCSet.Inv(mck(0)(0),temp_d3);
    Array1D<double> temp_d4(nPCTerms,0.e0); 
    myPCSet.Prod(temp_d,temp_d3,temp_d4);
    dv(0)=temp_d4;
    // dof dof
    Array1D<double> temp_d1(force(dof-1));
    subtractVec(kbaru(dof-1),temp_d1);
    addVec(kbaru_plus(dof-2),temp_d1);
    subtractVec(cv(dof-1),temp_d1);
    addVec(cv_plus(dof-2),temp_d1);
    //prodVal(temp_d1,1/mck(0)(dof-1));
    myPCSet.Inv(mck(0)(dof-1),temp_d3);
    myPCSet.Prod(temp_d1,temp_d3,temp_d4);
    dv(dof-1)=temp_d4;
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
        //prodVal(temp_d2,1/mck(0)(i));
        Array1D<double> temp_d5(nPCTerms,0.e0); 
        myPCSet.Inv(mck(0)(i),temp_d5);
        Array1D<double> temp_d6(nPCTerms,0.e0); 
        myPCSet.Prod(temp_d2,temp_d5,temp_d6);
        dv(i)=temp_d6;
    }
    
    return;
}

Array2D<double> postprocess_nGS(int dof, int nStep,  Array1D<Array2D<double> >& solution, PCSet& myPCSet, double dTym, int ord,Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array1D<Array2D<double> >& mstd_dis, Array1D<Array2D<double> >& mstd_vel, Array2D<double>& e2, Array1D<double>& e3){
    int nPCTerms = myPCSet.GetNumberPCTerms();
    // Open files to write out solutions
    ostringstream s1;
    s1 << "m_GS" << ord<<".dat";
    string m(s1.str());
    ostringstream s2;
    s2 << "s_GS" << ord<<".dat";
    string s(s2.str());
    ostringstream s3;

    // Output solution (mean and std) 
    Array1D<double> temp_u(nPCTerms,0.e0);
    Array1D<double> temp_v(nPCTerms,0.e0);
    Array2D<double> StDv_u(nStep+1,dof,0.e0);
    Array2D<double> StDv_v(nStep+1,dof,0.e0);
    Array2D<double> mean_u(nStep+1,dof,0.e0);
    Array2D<double> mean_v(nStep+1,dof,0.e0);
    Array2D<double> mean(nStep+1,2*dof,0.e0);
    Array2D<double> StDv(nStep+1,2*dof,0.e0);
    for (int j=0;j<dof;j++){
        for (int i=0;i<nStep+1;i++){
            for (int k=0;k<nPCTerms;k++){
                temp_v(k)=solution(j)(i,k);
                temp_u(k)=solution(j)(i,nPCTerms+k);
            }
            StDv_u(i,j) = myPCSet.StDv(temp_u);
            StDv_v(i,j) = myPCSet.StDv(temp_v);
            StDv(i,j) = StDv_v(i,j);
            StDv(i,dof+j) = StDv_u(i,j);
        }
        Array1D<double> temp2(nStep+1,0.e0);
        getCol(solution(j),0,temp2);
        mean_v.replaceCol(temp2,j);
        mean.replaceCol(temp2,j);
        mstd_vel(j).replaceRow(temp2,0);
        getCol(solution(j),nPCTerms,temp2);
        mean_u.replaceCol(temp2,j);
        mean.replaceCol(temp2,dof+j);
        mstd_dis(j).replaceRow(temp2,0);
        getCol(StDv_u,j,temp2);
        mstd_dis(j).replaceRow(temp2,1);
        getCol(StDv_v,j,temp2);
        mstd_vel(j).replaceRow(temp2,1);
    }
    
    write_datafile(mean,m.c_str());
    write_datafile(StDv,s.c_str());

    Array1D<Array2D<double> > et(dof);
    
    ostringstream info;
    info << "GS_"<<ord << "_d"<< dTym;
    string info_str = info.str();
    Array2D<double> e =  nerror(info_str, dof, nStep, et,mean,StDv,mean_MCS,std_MCS, e2, e3);

    return(e);
}
    
Array2D<double> nerror(string info, int dof, int nStep, Array1D<Array2D<double> >& et,Array2D<double>& mean,Array2D<double>& StDv, Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array2D<double>& e2, Array1D<double>& e3_norm){
    Array2D<double> e(dof,4,0.e0);
    Array2D<double> temp(dof,4,0.e0);

    for (int i=0;i<dof;i++){
        Array2D<double> error(nStep+1,4,0.e0);
        for (int ix=1;ix<nStep+1;ix++){
            error(ix,0)=abs(mean(ix,i)-mean_MCS(ix,i))/abs(mean_MCS(ix,i));
            error(ix,1)=abs(StDv(ix,i)-std_MCS(ix,i))/abs(std_MCS(ix,i));
            error(ix,2)=abs(mean(ix,i+dof)-mean_MCS(ix,i+dof))/abs(mean_MCS(ix,i+dof));
            error(ix,3)=abs(StDv(ix,i+dof)-std_MCS(ix,i+dof))/abs(std_MCS(ix,i+dof));
            e(i,0)+=error(ix,0);
            e(i,1)+=error(ix,1);
            e(i,2)+=error(ix,2);
            e(i,3)+=error(ix,3);
            e2(i,0)+=abs(mean(ix,i)-mean_MCS(ix,i));
            e2(i,1)+=abs(StDv(ix,i)-std_MCS(ix,i));
            e2(i,2)+=abs(mean(ix,i+dof)-mean_MCS(ix,i+dof));
            e2(i,3)+=abs(StDv(ix,i+dof)-std_MCS(ix,i+dof));
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
        name << "nerror_dof" << i << "_" << info << ".dat";
        string name_str = name.str();
        write_datafile(error,name_str.c_str());
    }

    Array2D<double> e3(4,nStep+1,0.e0);
    for (int ix=1;ix<nStep+1;ix++){
        for (int i=0;i<dof;i++){
            double temp1 = abs(mean(ix,i+dof)-mean_MCS(ix,i+dof));
            e3(0,ix)+=temp1*temp1;
            double temp2 = abs(StDv(ix,i+dof)-std_MCS(ix,i+dof));
            e3(1,ix)+=temp2*temp2;
            double temp3 = abs(mean_MCS(ix,i+dof));
            e3(2,ix)+=temp3*temp3;
            double temp4 = abs(std_MCS(ix,i+dof));
            e3(3,ix)+=temp4*temp4;
        }
        e3(0,ix)=sqrt(e3(0,ix));
        e3(1,ix)=sqrt(e3(1,ix));
        e3(2,ix)=sqrt(e3(2,ix));
        e3(3,ix)=sqrt(e3(3,ix));
    }
    Array1D<double> temp_e3(nStep+1,0.e0);
    getRow(e3,0,temp_e3); 
    Array1D<double> temp_e3_2(nStep+1,0.e0);
    getRow(e3,2,temp_e3_2); 
    e3_norm(0)=sum(temp_e3)/sum(temp_e3_2);
    getRow(e3,1,temp_e3); 
    getRow(e3,3,temp_e3_2); 
    e3_norm(1)=sum(temp_e3)/sum(temp_e3_2);
    
    ostringstream name2;
    name2 << "e2_" << info << ".dat";
    string name_str2 = name2.str();
    write_datafile(e2,name_str2.c_str());
    
    return(e);            
}
