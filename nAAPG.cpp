#include <math.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "PCBasis.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "Utils.h"
#include "nAAPG.h"
#include "AAPG.h"
#include "MCS.h"
#include "nMCS.h"
#include "lapack.h"
#include "nGhanemSpanos.h"
#include "ticktock.h"

void nAAPG(int refine, int dof, int nkl, int dim, int nStep, int order, Array1D<int>& lout, double factor_OD, int AAPG_ord, Array1D<double>& fbar, Array1D<double>& fbar_fine,double dTym, Array1D<double>& epsilon_mean, string pcType, Array2D<double>& scaledKLmodes,Array2D<double>& scaledKLmodes_fine, Array2D<double>& stat_e,  Array2D<double>& stat_i, Array2D<double>& stat_m, Array2D<double>& stat_c, Array2D<double>& stat_k, Array1D<double>& normsq, Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array1D<Array2D<double> >& mck, bool PDF, Array2D<double>& samPts_norm, bool active_D, double p, Array2D<double>& performance, int ord_GS){

    // timing var
    Array1D<double> t(5,0.e0);
    
    // Compute zeroth-order solution
    printf("Zeroth-order term...\n");
    // Allocate zeroth order solution
    Array2D<double> uv_0(nStep+1,2*dof,0.e0);
    // Define the determinisitic force
    Array2D<double> f0(2*nStep+1,dof,0.e0);
    for (int i=0;i<dof;i++){
        for (int ix=0;ix<2*nStep+1;ix++){
            f0(ix,i)=fbar(ix)*mck(0)(i,0);    
        }
    }
    // Define deterministic initial condition
    Array1D<double> initial(2*dof,0.e0);
    getCol(stat_i,0,initial);
    // Define deterministic epsilon
    Array1D<double> epsilon(dof,0.e0);
    getCol(stat_e,0,epsilon);
    // Define deterministic mck
    Array1D<Array1D<double> > mck_det(3);
    Array1D<double> m_det(dof,0.e0);
    getCol(stat_m,0,m_det); 
    Array1D<double> c_det(dof,0.e0);
    getCol(stat_c,0,c_det); 
    Array1D<double> k_det(dof,0.e0);
    getCol(stat_k,0,k_det);
    mck_det(0)=m_det;
    mck_det(1)=c_det;
    mck_det(2)=k_det; 
    // Deterministic solution
    TickTock tt;
    tt.tick();
    Array1D<Array1D<double> > temp=ndet(dof,nStep,dTym,f0,epsilon,mck_det,initial); 
    tt.tock("Took");
    t(0) = tt.silent_tock();
    // Reform the result to uv_0
    for (int id=0;id<2*dof;id++){
        uv_0.replaceCol(temp(id),id); 
    }
    write_datafile(uv_0,"uv_0.dat");
    
    // Compute first order terms
    printf("First-order terms...");
    // generate PCSet and the number of terms in it
    PCSet PCSet_1("ISP",order,1,pcType,0.0,1.0); 
    int PCTerms_1 = PCSet_1.GetNumberPCTerms();
    Array1D<Array1D<Array2D<double> > > force_1(dim);
    Array1D<Array1D<Array1D<double> > > epsilon_1(dim);
    Array1D<Array1D<Array1D<double> > > init_1(dim);
    Array1D<Array1D<Array1D<Array1D<double> > > > mck_1(dim);
    // Generate the forcing, epsilon and initial conditions on each stochastic dim
    for (int idim=0;idim<dim;idim++){
        // forcing
        Array1D<Array2D<double> > force_temp(dof);
        Array2D<double> f1_temp(2*nStep*refine+1,PCTerms_1,0.e0);
        f1_temp.replaceCol(fbar_fine,0);
        if (idim<nkl){
            Array1D<double> tempKL(2*nStep*refine+1,0.e0);
            getCol(scaledKLmodes_fine,idim,tempKL);
            f1_temp.replaceCol(tempKL,1);    
        }
        for (int id=0;id<dof;id++){
            force_temp(id) = f1_temp;
            prodVal(force_temp(id),mck(0)(id,0));            
        }
        force_1(idim)=force_temp;
        //epsilon
        Array1D<Array1D<double> > epsilon_temp(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> e1_temp(PCTerms_1,0.e0);
            if (idim==nkl+i)
                PCSet_1.InitMeanStDv(stat_e(i,0),stat_e(i,1),1,e1_temp);
            else
                PCSet_1.InitMeanStDv(stat_e(i,0),0,1,e1_temp);
            epsilon_temp(i)=e1_temp;
        }
        epsilon_1(idim)=epsilon_temp;
        // initial conditions
        Array1D<Array1D<double> > init_temp(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> i1_temp(2*PCTerms_1,0.e0);
            Array1D<double> temp_init2(PCTerms_1,0.e0);
            Array1D<double> temp_init3(PCTerms_1,0.e0);
            PCSet_1.InitMeanStDv(stat_i(i,0),0,1,temp_init2);
            PCSet_1.InitMeanStDv(stat_i(i+dof,0),0,1,temp_init3);
            if (idim==nkl+dof+i){
                PCSet_1.InitMeanStDv(stat_i(i,0),stat_i(i,1),1,temp_init2);
            }
            else if (idim==nkl+2*dof+i){
                PCSet_1.InitMeanStDv(stat_i(i+dof,0),stat_i(i+dof,1),1,temp_init3);
            }
            merge(temp_init2,temp_init3,i1_temp); 
            init_temp(i)= i1_temp;
        }
        init_1(idim)=init_temp;
        // mck
        Array1D<Array1D<Array1D<double> > > mck_1_dof(3);
        Array1D<Array1D<double> > m_GS(dof);
        Array1D<Array1D<double> > c_GS(dof);
        Array1D<Array1D<double> > k_GS(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> temp_m(PCTerms_1,0.e0);
            Array1D<double> temp_c(PCTerms_1,0.e0);
            Array1D<double> temp_k(PCTerms_1,0.e0);
            PCSet_1.InitMeanStDv(stat_m(i,0),0,1,temp_m);
            PCSet_1.InitMeanStDv(stat_c(i,0),0,1,temp_c);
            PCSet_1.InitMeanStDv(stat_k(i,0),0,1,temp_k);
            if (idim==nkl+3*dof+i){
                PCSet_1.InitMeanStDv(stat_m(i,0),stat_m(i,1),1,temp_m);
            }
            else if (idim==nkl+4*dof+i){
                PCSet_1.InitMeanStDv(stat_c(i,0),stat_c(i,1),1,temp_c);
            }
            else if (idim==nkl+5*dof+i){
                PCSet_1.InitMeanStDv(stat_k(i,0),stat_k(i,1),1,temp_k);
            }
            m_GS(i)=temp_m;
            c_GS(i)=temp_c;
            k_GS(i)=temp_k;
        }
        mck_1_dof(0)=m_GS;
        mck_1_dof(1)=c_GS;
        mck_1_dof(2)=k_GS;
        mck_1(idim)=mck_1_dof;
    }

    // allocate uv_1
    Array1D<Array1D<Array2D<double> > > uv_1(2*dof);
    for (int id=0;id<2*dof;id++){
        Array1D<Array2D<double> > temp(dim);
        for (int i=0;i<dim;i++){
            Array2D<double> temp2(nStep+1,PCTerms_1,0.e0);
            temp(i) = temp2;
        }
        uv_1(id)=temp;
    }
    // compute the std in each sub-dim and select the 'active' dim based on P5 of X. Yang et al. Adaptive ANOVA decomp. of stocha. imcompre
    Array1D<double> var(dim,0.e0); 
    tt.tick();
    #pragma omp parallel default(none) shared(refine,dof,PCSet_1,mck_1,nStep,dTym,PCTerms_1,epsilon_1,init_1,force_1,dim,uv_1,var)
    {
    #pragma omp for
    for (int i=0;i<dim;i++){
        Array1D<Array2D<double> >  uv_solution(dof);
        // MCS is assumed deterministic for now
        nGS(dof, PCSet_1, epsilon_1(i), mck_1(i), nStep*refine, init_1(i), dTym/refine, force_1(i), uv_solution);
        for (int id=0;id<dof;id++){
            for (int iPC=0;iPC<PCTerms_1;iPC++){
                for (int ix=0;ix<nStep+1;ix++){
                    uv_1(id)(i)(ix,iPC)=uv_solution(id)(ix*refine,iPC);
                    uv_1(dof+id)(i)(ix,iPC)=uv_solution(id)(ix*refine,PCTerms_1+iPC);
                }
            }
        }
        // adaptive AAPG routines
        Array1D<double> sol_temp(2*PCTerms_1,0.e0);
        Array1D<double> dis_temp(PCTerms_1,0.e0);
        for (int it=0;it<nStep+1;it++){
            for (int id=0;id<dof;id++){
                getRow(uv_solution(id),it,sol_temp);
                    for (int iPC=0;iPC<PCTerms_1;iPC++) dis_temp(iPC) = sol_temp(PCTerms_1+iPC);
                double temp =PCSet_1.StDv(dis_temp);  
                var(i)=var(i)+temp*temp*dTym;
            }  
        }
    }
    }
    tt.tock("Took");
    t(1)=tt.silent_tock();
    write_datafile_1d(var,"AAPG1var_ndof.dat");

    // sort ind and find the active dims based on threshold p
    Array1D<int> ind(dim,0);
    for (int i=0;i<dim;i++) ind(i)=i;
    if (active_D){
        // sort the variance in ascending order
        shell_sort_ind(var,ind);
        cout << "Var is sorted now" << endl;
        double var_sum = sum(var);
        double var_temp = 0.e0;
        int i_temp = 0;
        while ((var_temp+=var(i_temp))<((1-p)*var_sum)&&i_temp<dim){
            ind.erase(0);
            i_temp++;
        }
        shell_sort(ind);
    }
    int N_adof = ind.XSize();
    cout << "There are " << N_adof << "/" << dim << " active Dims:" << endl;
    for (int i=0;i<N_adof;i++) cout << ind(i) << endl;
    cout << endl;
    write_datafile_1d(ind,"ind2.dat");

    // Second order term
    Array1D<int> indi_2(N_adof*(N_adof-1)/2,0);
    Array1D<int> indj_2(N_adof*(N_adof-1)/2,0);
    printf("Second-order terms...");
    // generate PCSet and the number of terms in it
    PCSet PCSet_2("ISP",order,2,pcType,0.0,1.0); 
    int PCTerms_2 = PCSet_2.GetNumberPCTerms();
    int N2 = N_adof*(N_adof-1)/2;
    Array1D<Array1D<Array2D<double> > > force_2(N2);
    Array1D<Array1D<Array1D<double> > > epsilon_2(N2);
    Array1D<Array1D<Array1D<double> > > init_2(N2);
    Array1D<Array1D<Array1D<Array1D<double> > > > mck_2(N2);
    // Generate the forcing, epsilon and initial conditions on each stochastic dim
    int k=0;
    for (int i=0;i<N_adof-1;i++){
        int idim = ind(i);
        for (int i2=i+1;i2<N_adof;i2++){
        int idim2 = ind(i2);
        // forcing
        Array1D<Array2D<double> > force_temp(dof);
        Array2D<double> f2_temp(2*nStep+1,PCTerms_2,0.e0);
        f2_temp.replaceCol(fbar,0);
        if (idim<nkl){
            Array1D<double> tempKL(2*nStep+1,0.e0);
            getCol(scaledKLmodes,idim,tempKL);
            f2_temp.replaceCol(tempKL,1);    
        }
        if (idim2<nkl){
            Array1D<double> tempKL(2*nStep+1,0.e0);
            getCol(scaledKLmodes,idim2,tempKL);
            f2_temp.replaceCol(tempKL,2);    
        }
        for (int id=0;id<dof;id++){
            force_temp(id) = f2_temp;
            prodVal(force_temp(id),mck(0)(id,0));            
        }
        force_2(k)=force_temp;        
        //epsilon
        Array1D<Array1D<double> > epsilon_temp(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> e2_temp(PCTerms_2,0.e0);
            e2_temp(0)=stat_e(i,0);
            //PCSet_2.InitMeanStDv(stat_e(i,0),0,1,e2_temp);
            if (idim==nkl+i)
                PCSet_2.InitMeanStDv(stat_e(i,0),stat_e(i,1),1,e2_temp);
            if(idim2==nkl+i)
                PCSet_2.InitMeanStDv(stat_e(i,0),stat_e(i,1),2,e2_temp);
            epsilon_temp(i)=e2_temp;
        }
        epsilon_2(k)=epsilon_temp;
        // initial conditions
        Array1D<Array1D<double> > init_temp(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> i2_temp(2*PCTerms_2,0.e0);
            Array1D<double> temp_init2(PCTerms_2,0.e0);
            Array1D<double> temp_init3(PCTerms_2,0.e0);
            PCSet_2.InitMeanStDv(stat_i(i,0),0,1,temp_init2);
            PCSet_2.InitMeanStDv(stat_i(i+dof,0),0,1,temp_init3);
            if (idim==nkl+dof+i){
                PCSet_2.InitMeanStDv(stat_i(i,0),stat_i(i,1),1,temp_init2);
            }
            if (idim==nkl+2*dof+i){
                PCSet_2.InitMeanStDv(stat_i(i+dof,0),stat_i(i+dof,1),1,temp_init3);
            }
            if (idim2==nkl+dof+i){
                PCSet_2.InitMeanStDv(stat_i(i,0),stat_i(i,1),2,temp_init2);
            }
            if (idim2==nkl+2*dof+i){
                PCSet_2.InitMeanStDv(stat_i(i+dof,0),stat_i(i+dof,1),2,temp_init3);
            }
            merge(temp_init2,temp_init3,i2_temp); 
            init_temp(i)= i2_temp;
        }
        init_2(k) = init_temp;
	    indi_2(k) = idim;
        indj_2(k) = idim2;
        // mck
        Array1D<Array1D<Array1D<double> > > mck_2_dof(3);
        Array1D<Array1D<double> > m_GS(dof);
        Array1D<Array1D<double> > c_GS(dof);
        Array1D<Array1D<double> > k_GS(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> temp_m(PCTerms_2,0.e0);
            Array1D<double> temp_c(PCTerms_2,0.e0);
            Array1D<double> temp_k(PCTerms_2,0.e0);
            PCSet_2.InitMeanStDv(stat_m(i,0),0,1,temp_m);
            PCSet_2.InitMeanStDv(stat_c(i,0),0,1,temp_c);
            PCSet_2.InitMeanStDv(stat_k(i,0),0,1,temp_k);
            if (idim==nkl+3*dof+i){
                PCSet_2.InitMeanStDv(stat_m(i,0),stat_m(i,1),1,temp_m);
            }
            if (idim==nkl+4*dof+i){
                PCSet_2.InitMeanStDv(stat_c(i,0),stat_c(i,1),1,temp_c);
            }
            if (idim==nkl+5*dof+i){
                PCSet_2.InitMeanStDv(stat_k(i,0),stat_k(i,1),1,temp_k);
            }
            if (idim2==nkl+3*dof+i){
                PCSet_2.InitMeanStDv(stat_m(i,0),stat_m(i,1),2,temp_m);
            }
            if (idim2==nkl+4*dof+i){
                PCSet_2.InitMeanStDv(stat_c(i,0),stat_c(i,1),2,temp_c);
            }
            if (idim2==nkl+5*dof+i){
                PCSet_2.InitMeanStDv(stat_k(i,0),stat_k(i,1),2,temp_k);
            }
            m_GS(i)=temp_m;
            c_GS(i)=temp_c;
            k_GS(i)=temp_k;
        }
        mck_2_dof(0)=m_GS;
        mck_2_dof(1)=c_GS;
        mck_2_dof(2)=k_GS;
        mck_2(k)=mck_2_dof;
        k++;
        }
    }
    write_datafile_1d(indi_2,"indi2.dat");
    write_datafile_1d(indj_2,"indj2.dat");
    cout << "Finished generating the forcing, epsilon and initial conditions on each stochastic dim." << endl;

    // allocate uv_2
    Array1D<Array2D<Array2D<double> > > uv_2(2*dof);
    //Array2D<Array2D<double> > uv_2(dim,dim); 
    for (int id=0;id<2*dof;id++){
        Array2D<Array2D<double> > temp(dim,dim);
        for (int i=0;i<dim;i++){
            for (int j=0;j<dim;j++){
                Array2D<double> temp2(nStep+1,PCTerms_2,0.e0);
                temp(i,j) = temp2;
            }
        }
        uv_2(id)=temp;
    }
    tt.tick();
    #pragma omp parallel default(none) shared(N2,dof,PCSet_2,mck_2,nStep,dTym,PCTerms_2,epsilon_2,init_2,force_2,dim,uv_2,indi_2,indj_2)
    {
    #pragma omp for
    for (int i=0;i<N2;i++){
        Array1D<Array2D<double> >  uv_solution(dof);
        // MCS is assumed deterministic for now
        nGS(dof, PCSet_2, epsilon_2(i), mck_2(i), nStep, init_2(i), dTym, force_2(i), uv_solution);
        //printf("Finished nGS on AAPG2 dof %d\n",i);
        for (int id=0;id<dof;id++){
            for (int iPC=0;iPC<PCTerms_2;iPC++){
                for (int ix=0;ix<nStep+1;ix++){
                    uv_2(id)(indi_2(i),indj_2(i))(ix,iPC)=uv_solution(id)(ix,iPC);
                    uv_2(dof+id)(indi_2(i),indj_2(i))(ix,iPC)=uv_solution(id)(ix,PCTerms_2+iPC);
                }
            }
        }
    }
    }
    tt.tock("Took");
    t(2)=tt.silent_tock();

    // Third order term
    Array3D<Array2D<double> > uv_3(dim,dim,dim); 
    int PCTerms_3 = 0;
    Array1D<int> indi_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indj_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indk_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
 
    // initialize the first/second order mean/std of the AAPG solutions
    Array2D<double> m1(nStep+1,2*dof,0.e0);
    Array2D<double> m2(nStep+1,2*dof,0.e0);
    Array2D<double> m3(nStep+1,2*dof,0.e0);
    Array2D<double> std1(nStep+1,2*dof,0.e0);
    Array2D<double> std2(nStep+1,2*dof,0.e0);
    Array2D<double> std3(nStep+1,2*dof,0.e0);
 
    printf("Assemble the solutions...\n");
    tt.tick();
     #pragma omp parallel default(none) shared(indi_2,indj_2,indi_3,indj_3,indk_3,AAPG_ord,uv_1,uv_2,uv_3,dim,nStep,PCTerms_1,PCTerms_2,PCTerms_3,order,dTym,factor_OD,samPts_norm,lout,PDF,pcType,stat_i,dof,uv_0,mean_MCS,std_MCS,m1,m2,m3,std1,std2,std3)
    {
    #pragma omp for
    for (int i=0;i<dof*2;i++){
        //cout << "DOF " << i << endl;
        printf("DOF %d\n",i);
        ostringstream name1;
        name1<< "disvel_dof" << i;
        string name = name1.str();
        Array1D<double> uv0_temp(nStep+1,0.e0);
        getCol(uv_0,i,uv0_temp);
        Array1D<double> m1_temp(uv0_temp);
        Array1D<double> m2_temp(uv0_temp);
        Array1D<double> m3_temp(uv0_temp);
        Array1D<double> s1_temp(nStep+1,0.e0);
        Array1D<double> s2_temp(nStep+1,0.e0);
        Array1D<double> s3_temp(nStep+1,0.e0);
        Array2D<double> mstd_MCS(nStep+1,2,0.e0);
        Array1D<double> MCS_temp(nStep+1,0.e0);
        Array1D<Array1D<double> > e_sample(AAPG_ord);
        getCol(mean_MCS,i,MCS_temp);
        mstd_MCS.replaceCol(MCS_temp,0);
        getCol(std_MCS,i,MCS_temp);
        mstd_MCS.replaceCol(MCS_temp,1);
        PostProcess(indi_2,indj_2, indi_3, indj_3, indk_3, AAPG_ord, uv0_temp, uv_1(i), uv_2(i), uv_3, m1_temp, m2_temp, m3_temp, s1_temp, s2_temp, s3_temp, dim, nStep, PCTerms_1, PCTerms_2, PCTerms_3, order, dTym, factor_OD, mstd_MCS, samPts_norm, name, lout, e_sample, PDF, pcType,stat_i);
        for (int j=0;j<nStep+1;j++){
            m1(j,i)=m1_temp(j);
            m2(j,i)=m2_temp(j);
            m3(j,i)=m3_temp(j);
            std1(j,i)=s1_temp(j);
            std2(j,i)=s2_temp(j);
            std3(j,i)=s3_temp(j);
        }
    }
    }
    tt.tock("Took");
    t(4)=tt.silent_tock();
    
    // Compute the error compare to MCS
    Array2D<double> e2_AAPG1(dof,4,0.e0);
    Array1D<double> e3_AAPG1(2,0.e0);
    string info = "AAPG1";
    Array1D<Array2D<double> > et_AAPG1(dof);
    Array2D<double> e1_AAPG1 =  nerror(info,dof,nStep,et_AAPG1,m1,std1,mean_MCS,std_MCS,e2_AAPG1,e3_AAPG1);
    string info2 = "AAPG2";
    Array2D<double> e2_AAPG2(dof,4,0.e0);
    Array1D<double> e3_AAPG2(2,0.e0);
    Array1D<Array2D<double> > et_AAPG2(dof);
    Array2D<double> e1_AAPG2 =  nerror(info2,dof,nStep,et_AAPG2,m2,std2,mean_MCS,std_MCS,e2_AAPG2,e3_AAPG2);

    // print out the error
    cout << "AAPG1 Error kind 1 is" << endl;
    for (int i=0;i<dof;i++){
        cout << "Dof " << i << ", m_v:" << e1_AAPG1(i,0) << ",s_v:" << e1_AAPG1(i,1) << "," <<",m_u:" << e1_AAPG1(i,2) << ",s_u:"<<e1_AAPG1(i,3) << "." << endl;
    }
    cout << "Error kind 2 is" << endl;
    for (int i=0;i<dof;i++){
        cout << "Dof " << i << ", m_v:" << e2_AAPG1(i,0) << ",s_v:" << e2_AAPG1(i,1) << "," <<",m_u:" << e2_AAPG1(i,2) << ",s_u:"<<e2_AAPG1(i,3) << "." << endl;
    }
    cout << "AAPG2 Error kind 1 is" << endl;
    for (int i=0;i<dof;i++){
        cout << "Dof " << i << ", m_v:" << e1_AAPG2(i,0) << ",s_v:" << e1_AAPG2(i,1) << "," <<",m_u:" << e1_AAPG2(i,2) << ",s_u:"<<e1_AAPG2(i,3) << "." << endl;
    }
    cout << "Error kind 2 is" << endl;
    for (int i=0;i<dof;i++){
        cout << "Dof " << i << ", m_v:" << e2_AAPG2(i,0) << ",s_v:" << e2_AAPG2(i,1) << "," <<",m_u:" << e2_AAPG2(i,2) << ",s_u:"<<e2_AAPG2(i,3) << "." << endl;
    }
    cout << "Error kind new is" << endl;
    cout << "E1_mean=" << e3_AAPG1(0) <<", E1_std=" << e3_AAPG1(1) << endl;
    cout << "E2_mean=" << e3_AAPG2(0) <<", E2_std=" << e3_AAPG2(1) << endl;
    write_datafile_1d(e3_AAPG1,"errornew_AAPG1.dat");
    write_datafile_1d(e3_AAPG2,"errornew_AAPG2.dat");
    performance(1+ord_GS,0) = e3_AAPG1(0);
    performance(1+ord_GS,1) = e3_AAPG1(1);
    performance(1+ord_GS,2) = t(0)+t(1);
    performance(2+ord_GS,0) = e3_AAPG2(0);
    performance(2+ord_GS,1) = e3_AAPG2(1);
    performance(2+ord_GS,2) = t(0)+t(1)+t(2);

    write_datafile(m1,"m_AAPG1.dat");
    write_datafile(std1,"s_AAPG1.dat");
    write_datafile(m2,"m_AAPG2.dat");
    write_datafile(std2,"s_AAPG2.dat");
    ostringstream name1;
    name1<< "e1_AAPG1_" << refine <<".dat";
    string name1_str = name1.str();
    write_datafile(e1_AAPG1,name1_str.c_str());
    ostringstream name2;
    name2 << "e2_AAPG1_" << refine <<".dat";
    string name2_str = name2.str();
    write_datafile(e2_AAPG1,name2_str.c_str());
    write_datafile(et_AAPG1(0),"etat0_AAPG1.dat");
    ostringstream name3;
    name3 << "e1_AAPG2_" << refine <<".dat";
    string name3_str = name3.str();
    write_datafile(e1_AAPG2,name3_str.c_str());
    ostringstream name4;
    name4 << "e2_AAPG2_" << refine << ".dat";
    string name4_str = name4.str();
    write_datafile(e2_AAPG2,name4_str.c_str());
    write_datafile(et_AAPG2(0),"etat0_AAPG2.dat");
 
    write_datafile_1d(t,"t_AAPG.dat");
    return;
}

