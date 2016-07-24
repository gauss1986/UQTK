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
            //cout << "Coeff_D="<< coeff_D(0) << "," << coeff_D(1)<<endl;
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
        if (abs(inpParams(0)-4)<1e-10){ //VDP
            Array1D<double> epsilon_vector(nPCTerms,0.e0);
            // parse input parameters
            const double epsilon_m = inpParams(1);
            double epsilon_s = 0;
            if (inpParams.XSize()==3){
                epsilon_s = inpParams(2);
            }
            myPCSet.InitMeanStDv(epsilon_m,epsilon_s,1,epsilon_vector);
            // buff to store temperary results 
            Array1D<double> temp(nPCTerms,0.e0);                     
            Array1D<double> temp2(nPCTerms,0.e0);                     
            // nonlinear term
            myPCSet.Prod(x(1),x(1),temp2);
            myPCSet.Prod(x(0),temp2,temp);
            myPCSet.Prod(temp,epsilon_vector,temp2);
            // damping term
            myPCSet.Prod(epsilon_vector,x(0),temp);
            for (int ind=0;ind<nPCTerms;ind++){
                dev(1)(ind) = x(0)(ind);//velocity term
                dev(0)(ind) = force(ind)-temp2(ind)+temp(ind)-x(1)(ind);
            }
        }
}

Array1D<double> postprocess_GS(int noutput, int nPCTerms, int nStep,  Array2D<double>& solution, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS, Array2D<double>& stat, Array2D<double>& et){
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
    Array1D<double> e =  error(et,mean, StDv, mstd_MCS);

    stat.replaceRow(mean,0);
    stat.replaceRow(StDv,1);    
    return e;
}

Array2D<double> sampleGS(Array1D<int>& lout, int dim, int nStep, int nPCTerms, PCSet& myPCSet, Array2D<double>& solution, Array2D<double>& samPts, Array2D<double>& stat, Array1D<double>& e_GS_sample){
    int nspl = samPts.XSize();
    int noutput = lout.XSize();
    Array2D<double> GS_sampt(nspl,noutput,0.e0);
    int ind = 0;
    Array1D<double> temp(nPCTerms,0.e0);
    Array1D<double> normsq(nPCTerms,0.e0); 
    myPCSet.OutputNormSquare(normsq);
    for (int i=0;i<noutput;i++){
            getRow(solution,lout(i),temp);
            #pragma omp parallel default(none) shared(normsq,nspl,dim,myPCSet,samPts,temp,GS_sampt,ind) 
            {
            #pragma omp for 
            for (int j=0;j<nspl;j++){
                Array1D<double> samPt(dim,0.e0);
                getRow(samPts,j,samPt);
                if (!strcmp(myPCSet.GetPCType().c_str(),"LU")){
                    GS_sampt(j,ind)=myPCSet.EvalPC(temp, samPt);
                }
                else if(!strcmp(myPCSet.GetPCType().c_str(),"HG")){
                    //GS_sampt(j,ind)=myPCSet.EvalPC(temp, samPt);
                    GS_sampt(j,ind)=temp(0)+temp(1)*samPt(0)+temp(2)*samPt(1);
                    if(myPCSet.GetOrder()==2)
                        GS_sampt(j,ind)+=temp(3)*(samPt(0)*samPt(0)-normsq(1))+temp(4)*samPt(0)*samPt(1)+temp(5)*(samPt(1)*samPt(1)-normsq(2));
                }
            }            
            }
            ind++;
    }

    // Compute error of the sampling results against the direct assembly results
    Array1D<double> m_sample_GS(noutput,0.e0);
    Array1D<double> s_sample_GS(noutput,0.e0);
    Array2D<double> mstd_GS_noutput(2,noutput,0.e0);
   for (int i=0;i<noutput;i++){  
        Array1D<double> sample_GS_atT(nspl,0.e0);
        getCol(GS_sampt,i,sample_GS_atT);
        Array1D<double> sample_mstd_GS_temp = mStd(sample_GS_atT,nspl);
        m_sample_GS(i)=sample_mstd_GS_temp(0);
        s_sample_GS(i)=sample_mstd_GS_temp(1);
        mstd_GS_noutput(0,i)=stat(0,lout(i));
        mstd_GS_noutput(1,i)=stat(1,lout(i));
    }
    Array2D<double> et(nStep+1,2,0.e0);
    e_GS_sample =  error(et,m_sample_GS, s_sample_GS, mstd_GS_noutput);

    return GS_sampt;
}
