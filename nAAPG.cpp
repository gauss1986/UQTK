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

Array1D<double> nAAPG(int dof, int nkl, int dim, int nStep, int order, int noutput, double factor_OD, int AAPG_ord, bool act_D, Array1D<Array1D<double> >& mck, Array1D<double>& fbar, double dTym, Array1D<double>& epsilon_mean, string pcType, Array2D<double>& scaledKLmodes,  Array2D<double>& stat_e,  Array2D<double>& stat_i, Array1D<double>& normsq, Array2D<double>& mean_MCS, Array2D<double>& std_MCS){
    Array2D<double> samPts_norm(2,2,0.e0);
    bool PDF = false;
    Array1D<Array1D<double> > e_sample;

    // timing var
    Array1D<double> t(AAPG_ord+1,0.e0);
    
    // Compute zeroth-order solution
    printf("Zeroth-order term...\n");
    // Allocate zeroth order solution
    Array2D<double> uv_0(nStep+1,2*dof,0.e0);
    // Define the determinisitic force
    Array2D<double> f0(2*nStep+1,dof,0.e0);
    for (int i=0;i<dof;i++){
        for (int ix=0;ix<2*nStep+1;ix++){
            f0(ix,i)=fbar(ix)*mck(0)(i);    
        }
    }
    // Define deterministic initial condition
    Array1D<double> initial(2*dof,0.e0);
    getCol(stat_i,0,initial);
    // Deterministic solution
    TickTock tt;
    tt.tick();
    Array1D<Array1D<double> > temp=ndet(dof,nStep,dTym,f0,epsilon_mean,mck,initial); 
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
    // Generate the forcing, epsilon and initial conditions on each stochastic dim
    for (int idim=0;idim<dim;idim++){
        // forcing
        Array1D<Array2D<double> > force_temp(dof);
        Array2D<double> f1_temp(2*nStep+1,PCTerms_1,0.e0);
        f1_temp.replaceCol(fbar,0);
        if (idim<nkl){
            Array1D<double> tempKL(2*nStep+1,0.e0);
            getCol(scaledKLmodes,idim,tempKL);
            f1_temp.replaceCol(tempKL,1);    
        }
        for (int id=0;id<dof;id++){
            force_temp(id) = f1_temp;
            prodVal(force_temp(id),mck(0)(id));            
        }
        force_1(idim)=force_temp;
        //epsilon
        Array1D<Array1D<double> > epsilon_temp(dof);
        for (int i=0;i<dof;i++){
            Array1D<double> e1_temp(PCTerms_1,0.e0);
            e1_temp(0)=stat_e(i,0);
            if (idim==nkl+i)
                PCSet_1.InitMeanStDv(stat_e(i,0),stat_e(i,1),1,e1_temp);
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
    tt.tick();
    #pragma omp parallel default(none) shared(dof,PCSet_1,mck,nStep,dTym,PCTerms_1,epsilon_1,init_1,force_1,dim,uv_1)
    {
    #pragma omp for
    for (int i=0;i<dim;i++){
        Array1D<Array2D<double> >  uv_solution(dof);
        // MCS is assumed deterministic for now
        nGS(dof, PCSet_1, epsilon_1(i), mck, nStep, init_1(i), dTym, force_1(i), uv_solution);
        for (int id=0;id<dof;id++){
            for (int iPC=0;iPC<PCTerms_1;iPC++){
                for (int ix=0;ix<nStep+1;ix++){
                    uv_1(id)(i)(ix,iPC)=uv_solution(id)(ix,iPC);
                    uv_1(dof+id)(i)(ix,iPC)=uv_solution(id)(ix,PCTerms_1+iPC);
                }
            }
        }
    }
    }
    tt.tock("Took");
    t(1)=tt.silent_tock();

    Array1D<int> ind(dim,0);
    if (!act_D){
        for (int i=0;i<dim;i++)
            ind(i)=i;
    }
    int N_adof = ind.XSize();

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
    // Generate the forcing, epsilon and initial conditions on each stochastic dim
    int k=0;
    for (int idim=0;idim<N_adof-1;idim++){
        for (int idim2=idim+1;idim2<N_adof;idim2++){
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
            prodVal(force_temp(id),mck(0)(id));            
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
	    indi_2(k) = ind(idim);
        indj_2(k) = ind(idim2);
        k++;
        }
    }
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
    #pragma omp parallel default(none) shared(N2,dof,PCSet_2,mck,nStep,dTym,PCTerms_2,epsilon_2,init_2,force_2,dim,uv_2,indi_2,indj_2)
    {
    #pragma omp for
    for (int i=0;i<N2;i++){
        Array1D<Array2D<double> >  uv_solution(dof);
        // MCS is assumed deterministic for now
        nGS(dof, PCSet_2, epsilon_2(i), mck, nStep, init_2(i), dTym, force_2(i), uv_solution);
        //ostringstream name_AAPG2;
        //name_AAPG2 << "uv2_dim" << i << ".dat";
        //string name_AAPG2str = name_AAPG2.str();
        //write_datafile(uv_solution(0),name_AAPG2str.c_str());
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
    string name = "disvel";
    tt.tick();
    for (int i=0;i<dof*2;i++){
        cout << "DOF " << i << endl;
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
        getCol(mean_MCS,i,MCS_temp);
        mstd_MCS.replaceCol(MCS_temp,0);
        getCol(std_MCS,i,MCS_temp);
        mstd_MCS.replaceCol(MCS_temp,1);
        PostProcess(indi_2,indj_2, indi_3, indj_3, indk_3, AAPG_ord, uv0_temp, uv_1(i), uv_2(i), uv_3, m1_temp, m2_temp, m3_temp, s1_temp, s2_temp, s3_temp,  normsq, dim, nStep, PCTerms_1, PCTerms_2, PCTerms_3, order, dTym, factor_OD, mstd_MCS, samPts_norm, name, noutput, e_sample, PDF);
        m1.replaceCol(m1_temp,i);
        m2.replaceCol(m2_temp,i);
        m3.replaceCol(m3_temp,i);
        std1.replaceCol(s1_temp,i);
        std2.replaceCol(s2_temp,i);
        std3.replaceCol(s3_temp,i);
    }
    tt.tock("Took");
    t(AAPG_ord+1)=tt.silent_tock();
    
    // Compute the error compare to MCS
    Array2D<double> e2_AAPG1(dof,4,0.e0);
    ostringstream info;
    info << "AAPG1_"<< "GS_"<<order << "d"<< dTym;
    string info_str = info.str();
    Array1D<Array2D<double> > et_AAPG1(dof);
    Array2D<double> e1_AAPG1 =  nerror(info_str,dof,nStep,et_AAPG1,m1,std1,mean_MCS,std_MCS,e2_AAPG1);
    ostringstream info2;
    info2 << "AAPG2_"<< "GS_"<<order << "_d"<< dTym;
    string info2_str = info2.str();
    Array2D<double> e2_AAPG2(dof,4,0.e0);
    Array1D<Array2D<double> > et_AAPG2(dof);
    Array2D<double> e1_AAPG2 =  nerror(info2_str,dof,nStep,et_AAPG2,m2,std2,mean_MCS,std_MCS,e2_AAPG2);

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

    write_datafile(m1,"m1.dat");
    write_datafile(std1,"s1.dat");
    write_datafile(m2,"m2.dat");
    write_datafile(std2,"s2.dat");
    write_datafile(e1_AAPG1,"e1_AAPG1.dat");
    write_datafile(e2_AAPG1,"e2_AAPG1.dat");
    write_datafile(et_AAPG1(0),"etat0_AAPG1.dat");
    write_datafile(e1_AAPG2,"e1_AAPG2.dat");
    write_datafile(e2_AAPG2,"e2_AAPG2.dat");
    write_datafile(et_AAPG2(0),"etat0_AAPG2.dat");
 
    write_datafile_1d(t,"t_nAAPG.dat");
    return(t);
}

