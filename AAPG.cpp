#include <math.h>
#include <cmath>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "PCBasis.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "Utils.h"
#include "AAPG.h"
#include "MCS.h"
#include "lapack.h"
#include "GhanemSpanos.h"
#include "ticktock.h"

Array1D<double> AAPG(int dof, Array1D<double> inpParams, Array1D<double>& fbar, double dTym, int order, string pcType, int noutput, int dim, int nStep, Array2D<double>& scaledKLmodes, Array1D<double>& normsq, double factor_OD, int AAPG_ord, bool act_D, double p, Array2D<double>& mstd_MCS, FILE* err_dump, Array2D<double>& sample_mstd_2D, Array2D<double>& samPts_norm, Array2D<double>& e_AAPG, Array1D<Array1D<double> >& e_sample_dis, Array1D<Array1D<double> >& e_sample_vel, Array1D<int>& init_D, Array1D<int>& coeff_D, bool PDF){
    // timing var
    Array1D<double> t(6,0.e0);
    
    // Compute zeroth-order term with DET initial condition and force
    printf("Zeroth-order term...\n");
    // Variable to store solution at every step
    Array1D<double> Temp(2,0.e0);
    getCol(sample_mstd_2D,0,Temp); 
    // store zeroth order solution
    Array2D<double> x_0(nStep+1,2,0.e0);
    // Initialize x_0
    x_0.replaceRow(Temp,0);
    // Time marching steps
    TickTock tt;
    tt.tick();
    Array1D<double> tempf(3,0.e0);
    for (int ix=0;ix<nStep;ix++){
        tempf(0)=fbar(2*ix);
        tempf(1)=fbar(2*ix+1);
        tempf(2)=fbar(2*ix+2);
        forward_duffing_dt(inpParams,tempf,dTym,Temp);
        x_0.replaceRow(Temp,ix+1);
    }
    tt.tock("Took");
    t(0) = tt.silent_tock();
    // abstract the displacement terms
    Array1D<double> dis_0(nStep+1,0.e0);
    getCol(x_0,1,dis_0);
    Array1D<double> vel_0(nStep+1,0.e0);
    getCol(x_0,0,vel_0);
    
    // Compute first order terms
    printf("First-order terms...");
    // generate PCSet and the number of terms in it
    PCSet PCSet_1("ISP",order,1,pcType,0.0,1.0); 
    int PCTerms_1 = PCSet_1.GetNumberPCTerms();
    // initialize the displacement and force terms
    Array1D<Array2D<double> > dis_1(dim);
    Array1D<Array2D<double> > vel_1(dim);
    Array1D<Array2D<double> > force_1(dim);
    // Generate the forcing on each dim
    for (int i=0;i<dim;i++){
        Array2D<double> f_1(2*nStep+1,PCTerms_1,0.e0);
        for (int it=0;it<2*nStep+1;it++){
            f_1(it,0) = fbar(it);
            f_1(it,1) = scaledKLmodes(it,i);
        }
        force_1(i) = f_1;
    }
  
    tt.tick();
    #pragma omp parallel default(none) shared(sample_mstd_2D,dof, dim,PCSet_1,order,PCTerms_1,nStep,dTym,inpParams,force_1,vel_1,dis_1,init_D, coeff_D)
    {
    #pragma omp for
    for (int i=0;i<dim;i++){
        Array1D<Array2D<double> > temp(dof);
        Array1D<double> inpParams_1(5,0.e0);
        // Initial condition
        Array1D<double> temp_init(PCTerms_1,0.e0);
        Array1D<Array1D<double> > initial_GS1(dof);
        for (int j=0;j<dof;j++){
            initial_GS1(j) = temp_init;
            initial_GS1(j)(0) = sample_mstd_2D(j,0);
        }
        for (int j=0;j<dof;j++){
            if (i==init_D(j)){
                PCSet_1.InitMeanStDv(sample_mstd_2D(j,0),sample_mstd_2D(j,1),1,initial_GS1(j));        
            }
        }
        // coeffs on this ANOVA dim
        inpParams_1(0)=inpParams(0);
        inpParams_1(1)=inpParams(1);
        inpParams_1(2)=inpParams(2);
        for (int j=0;j<2;j++){
            if (i==coeff_D(j)){
               inpParams_1(j+3)=inpParams(j+3); 
            }
        }
        Array1D<int> active_1D(2,0);
        GS(dof, PCSet_1, active_1D, PCTerms_1, nStep, initial_GS1, dTym, inpParams_1, force_1(i),temp);
        dis_1(i) = temp(1);
        vel_1(i) = temp(0);
    }
    }
    tt.tock("Took");
    t(1)=tt.silent_tock();

    Array1D<int> ind(dim,0);
    if (!act_D){
        for (int i=0;i<dim;i++)
            ind(i)=i;
    }
    else{
    // compute the std in each sub-dim and select the 'active' dim based on P5 of X. Yang et al. Adaptive ANOVA decomp. of stocha. imcompre
    Array1D<double> var(dim,0.e0); 
    for (int i=0;i<dim;i++){
        Array1D<double> dis_temp(PCTerms_1,0.e0);
        ind(i) = i;
        for (int it=0;it<nStep+1;it++){
            getRow(dis_1(i),it,dis_temp);
            var(i)=var(i)+PCSet_1.StDv(dis_temp);  
        }
        cout << "Sum of var on dim No." <<i <<" is "<< var(i) << endl;
    }
    // sort the variance in ascending order
    shell_sort_ind(var,ind);
    cout << "Var is sorted now" << endl;
    for (int i=0;i<dim;i++){
        cout << "Sum of var on dim NO. " <<ind(i) <<" is "<< var(i) << endl;
    }
    // select the active dims
    double var_sum = sum(var);
    double temp = 0.e0;
    int i_temp = 0;
    while ((temp+=var(i_temp))<((1-p)*var_sum)&&i_temp<dim){
        cout << "Dim No." << ind(0) << " is non-active." << " Percentage:"<<temp/var_sum*100<<"%"<< endl;
        ind.erase(0);
        i_temp++;
    }
    shell_sort(ind);
    }
    int N_adof = ind.XSize();

    // Second order term
    Array2D<Array2D<double> > vel_2(dim,dim); 
    Array2D<Array2D<double> > dis_2(dim,dim); 
    int PCTerms_2 = 0;
    Array1D<int> indi_2(N_adof*(N_adof-1)/2,0);
    Array1D<int> indj_2(N_adof*(N_adof-1)/2,0);
    if (AAPG_ord >= 2){
        printf("Second-order terms...");
        PCSet PCSet_2("ISP",order,2,pcType,0.0,1.0); 
        PCTerms_2 = PCSet_2.GetNumberPCTerms();
        // Initial condition
        Array1D<Array1D<Array1D<double> > > initial_GS2(N_adof*(N_adof-1)/2);
        Array1D<Array1D<double> > temp2(dof);
        Array1D<double> temp_init2(PCTerms_2,0.e0);
        Array1D<Array2D<double> > force_2(N_adof*(N_adof-1)/2);
        Array1D<Array1D<double> > inpParams_2(N_adof*(N_adof-1)/2);
        Array1D<double> temp_coeff2(5,0.e0);
        int k = 0;
        Array1D<Array1D<int> > active_2D(N_adof*(N_adof-1)/2);
        Array1D<int> temp_active(2,0);
        for (int i=0;i<N_adof-1;i++){
            for (int j=i+1;j<N_adof;j++){
                initial_GS2(k)=temp2;
                // Initial condition
                for (int L=0;L<dof;L++){
                    initial_GS2(k)(L) = temp_init2;
                    initial_GS2(k)(L)(0) = sample_mstd_2D(L,0);
                }
                for (int l=0;l<dof;l++){
                    if (i==init_D(l))
                        PCSet_2.InitMeanStDv(sample_mstd_2D(l,0),sample_mstd_2D(l,1),1,initial_GS2(k)(l));
                    if (j==init_D(l))
                        PCSet_2.InitMeanStDv(sample_mstd_2D(l,0),sample_mstd_2D(l,1),2,initial_GS2(k)(l));
                }
                // force
                Array2D<double> f_2(2*nStep+1,PCTerms_2,0.e0);
                for (int it=0;it<2*nStep+1;it++){
                    f_2(it,0) = fbar(it);
                    f_2(it,1) = scaledKLmodes(it,ind(i));
                    f_2(it,2) = scaledKLmodes(it,ind(j));
                }
                // coeffs on this ANOVA dim
                inpParams_2(k)=temp_coeff2;
                inpParams_2(k)(0)=inpParams(0);
                inpParams_2(k)(1)=inpParams(1);
                inpParams_2(k)(2)=inpParams(2);
                active_2D(k)= temp_active;
                for (int l=0;l<2;l++){
                    if (i==coeff_D(l)){
                        inpParams_2(k)(l+3)=inpParams(l+3);
                        active_2D(k)(l)=0; 
                    }
                    if (j==coeff_D(l)){
                        inpParams_2(k)(l+3)=inpParams(l+3); 
                        active_2D(k)(l)=1; 
                    }
                }
	            force_2(k)=f_2;
	            indi_2(k) = ind(i);
    	        indj_2(k) = ind(j);
	            k++;
	    }
    }
    tt.tick();
    #pragma omp parallel  default(none) shared(dof,k,dis_2,vel_2,indi_2,indj_2,PCSet_2, PCTerms_2,nStep,initial_GS2,dTym,inpParams_2,force_2, active_2D)
    {
    #pragma omp for
    for (int i=0;i<k;i++){
        Array1D<Array2D<double> > temp(dof);
        GS(dof, PCSet_2, active_2D(i), PCTerms_2, nStep, initial_GS2(i), dTym, inpParams_2(i), force_2(i), temp); 
        vel_2(indi_2(i),indj_2(i)) = temp(0);
        dis_2(indi_2(i),indj_2(i)) = temp(1);
    }
    }
    tt.tock("Took");   
    t(2)=tt.silent_tock();
    cout << "Saving on computational cost using only active dims: " << (1-(N_adof*(N_adof-1)*1.0/dim/(dim-1)))*100 << "%" << endl;
    }
 
    // Third order term
    Array3D<Array2D<double> > dis_3(dim,dim,dim); 
    Array3D<Array2D<double> > vel_3(dim,dim,dim); 
    int PCTerms_3 = 0;
    Array1D<int> indi_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indj_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indk_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    if (AAPG_ord >= 3){
    printf("Third-order terms...");
    PCSet PCSet_3("ISP",order,3,pcType,0.0,1.0); 
    PCTerms_3 = PCSet_3.GetNumberPCTerms();
    // Initial condition
    Array1D<Array1D<Array1D<double> > > initial_GS3(N_adof*(N_adof-1)*(N_adof-2)/6);
    Array1D<Array1D<double> > temp3(dof);
    Array1D<double> temp_init3(PCTerms_3,0.e0);
    Array1D<Array2D<double> > force_3(N_adof*(N_adof-1)*(N_adof-2)/6);
    Array1D<double> temp_coeff3(5,0.e0);
    Array1D<Array1D<double> > inpParams_3(N_adof*(N_adof-1)*(N_adof-2)/6);
    Array1D<Array1D<int> > active_3D(N_adof*(N_adof-1)*(N_adof-2)/6);
    Array1D<int> temp_active(2,0);
    int l=0;
    for (int i=0;i<N_adof-2;i++){
        for (int j=i+1;j<N_adof-1;j++){
            for (int k=j+1;k<N_adof;k++){
                initial_GS3(l)=temp3;
                // Initial condition
                for (int L=0;L<dof;L++){
                    initial_GS3(l)(L) = temp_init3;
                    initial_GS3(l)(L)(0) = sample_mstd_2D(L,0);
                }
                for (int L=0;L<dof;L++){
                    if (i==init_D(L))
                        PCSet_3.InitMeanStDv(sample_mstd_2D(L,0),sample_mstd_2D(L,1),1,initial_GS3(l)(L));
                    if (j==init_D(L))
                        PCSet_3.InitMeanStDv(sample_mstd_2D(L,0),sample_mstd_2D(L,1),2,initial_GS3(l)(L));
                    if (k==init_D(L))
                        PCSet_3.InitMeanStDv(sample_mstd_2D(L,0),sample_mstd_2D(L,1),3,initial_GS3(l)(L));
                }
                // force
                Array2D<double> f_3(2*nStep+1,PCTerms_3,0.e0);
                for (int it=0;it<2*nStep+1;it++){
                    f_3(it,0) = fbar(it);
                    f_3(it,1) = scaledKLmodes(it,ind(i));
                    f_3(it,2) = scaledKLmodes(it,ind(j));
                    f_3(it,3) = scaledKLmodes(it,ind(k));
                }
		        force_3(l)=f_3;
    		    indi_3(l)=ind(i);
		        indj_3(l)=ind(j);
		        indk_3(l)=ind(k);
                // parameters
                active_3D(l)= temp_active;
                inpParams_3(l)=temp_coeff3;
                inpParams_3(l)(0)=inpParams(0);
                inpParams_3(l)(1)=inpParams(1);
                inpParams_3(l)(2)=inpParams(2);
                for (int L=0;L<2;L++){
                    if (i==coeff_D(L)){
                        inpParams_3(l)(L+3)=inpParams(L+3);
                        active_3D(l)(L)=0; 
                    }
                    if (j==coeff_D(L)){
                        inpParams_3(l)(L+3)=inpParams(L+3); 
                        active_3D(l)(L)=1; 
                    }
                    if (k==coeff_D(L)){
                        inpParams_3(l)(L+3)=inpParams(L+3);
                        active_3D(l)(L)=2; 
                    }
                }
		        l++;
                }
            }
	    }
    //start = clock();
    tt.tick();
    cout << "Finished generating initial conditions and excition force. Now starting parallel computing of sub-problems..."<< endl<<flush;
    #pragma omp parallel default(none) shared(dof, l,indi_3,indj_3,indk_3,PCSet_3, PCTerms_3,nStep,initial_GS3,dTym,inpParams_3,vel_3,dis_3,force_3,active_3D)
    {
    #pragma omp for
    for (int i=0;i<l;i++){
        Array1D<Array2D<double> > temp(dof);
        //Array1D<int> active_3D(2,0);
        //active_3D(0)=0;
        //active_3D(1)=1;
        GS(dof, PCSet_3, active_3D(i), PCTerms_3, nStep, initial_GS3(i), dTym, inpParams_3(i), force_3(i),temp);         
	    dis_3(indi_3(i),indj_3(i),indk_3(i)) = temp(1);
	    vel_3(indi_3(i),indj_3(i),indk_3(i)) = temp(0);
    }
    }
    tt.tock("Took");
    t(3)=tt.silent_tock();
    cout << "Saving on computational cost using only active dims: " << (1-(N_adof*(N_adof-1)*(N_adof-2)*1.0/dim/(dim-1)/(dim-2)))*100 << "%" << endl;
    }
 
    // initialize the first/second order mean/std of the AAPG solutions
    Array1D<double> dis_1_mean = dis_0;
    Array1D<double> dis_2_mean = dis_0;
    Array1D<double> dis_3_mean = dis_0;
    Array1D<double> std1(nStep+1,0.e0);
    Array1D<double> std2(nStep+1,0.e0);
    Array1D<double> std3(nStep+1,0.e0);
   
    printf("\n");
    printf("Assemble the solutions...\n");
    string name = "dis";
    printf("Dis...\n");
    tt.tick();
    PostProcess(indi_2,indj_2, indi_3, indj_3, indk_3, AAPG_ord, dis_0, dis_1, dis_2, dis_3, dis_1_mean, dis_2_mean, dis_3_mean, std1, std2, std3,  normsq, dim, nStep, PCTerms_1, PCTerms_2, PCTerms_3, order, dTym, factor_OD, mstd_MCS, samPts_norm, name, noutput, e_sample_dis, PDF, pcType, sample_mstd_2D);
    tt.tock("Took");
    t(4)=tt.silent_tock();
    printf("Vel...\n");
    string name2 = "vel";
    Array1D<double> vel_1_mean = vel_0;
    Array1D<double> vel_2_mean = vel_0;
    Array1D<double> vel_3_mean = vel_0;
    Array1D<double> std_vel_1(nStep+1,0.e0);
    Array1D<double> std_vel_2(nStep+1,0.e0);
    Array1D<double> std_vel_3(nStep+1,0.e0);
    tt.tick();
    PostProcess(indi_2,indj_2, indi_3, indj_3, indk_3, AAPG_ord, vel_0, vel_1, vel_2, vel_3, vel_1_mean, vel_2_mean, vel_3_mean, std_vel_1, std_vel_2, std_vel_3,  normsq, dim, nStep, PCTerms_1, PCTerms_2, PCTerms_3, order, dTym, factor_OD, mstd_MCS, samPts_norm, name2, noutput, e_sample_vel, PDF, pcType, sample_mstd_2D);
    tt.tock("Took");
    t(5)=tt.silent_tock();
   
    // print out mean/std valus at specific points for comparison
    printf("\n");
    printf("First-order AAPG results:\n");
    if (AAPG_ord >= 1){
    	for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/noutput) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_1_mean(ix),std1(ix));
            }
    	}
        Array2D<double> et1(nStep+1,2,0.e0);
        Array1D<double> e1 =  error(et1,dis_1_mean, std1, mstd_MCS);    
    	write_datafile_1d(dis_1_mean,"dis_1_mean.dat");
    	write_datafile_1d(std1,"dis_1_std.dat");
    	write_datafile_1d(vel_1_mean,"vel_1_mean.dat");
    	write_datafile_1d(std_vel_1,"vel_1_std.dat");
    	write_datafile(et1,"e_dis_AAPG1.dat");
        fprintf(err_dump,"%lg %lg\n",e1(0),e1(1)); 
        e_AAPG.replaceRow(e1,0);
    }
    if(AAPG_ord >= 2){
	printf("Second-order AAPG results:\n");
    	for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/noutput) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_2_mean(ix),std2(ix));
            }
    	}
        Array2D<double> et2(nStep+1,2,0.e0);
        Array1D<double> e2 =  error(et2,dis_2_mean, std2, mstd_MCS);    
    	write_datafile_1d(dis_2_mean,"dis_2_mean.dat");
    	write_datafile_1d(std2,"dis_2_std.dat");
    	write_datafile_1d(vel_2_mean,"vel_2_mean.dat");
    	write_datafile_1d(std_vel_2,"vel_2_std.dat");
    	write_datafile(et2,"e_dis_AAPG2.dat");
        fprintf(err_dump,"%lg %lg\n",e2(0),e2(1)); 
        e_AAPG.replaceRow(e2,1);
    }
    if(AAPG_ord >= 3){
        printf("Third-order AAPG results:\n");
        for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/noutput) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_3_mean(ix),std3(ix));
        	}
    	}
        Array2D<double> et3(nStep+1,2,0.e0);
        Array1D<double> e3 =  error(et3,dis_3_mean, std3, mstd_MCS);    
    	write_datafile_1d(dis_3_mean,"dis_3_mean.dat");
    	write_datafile_1d(std3,"dis_3_std.dat");
    	write_datafile_1d(vel_3_mean,"vel_3_mean.dat");
    	write_datafile_1d(std_vel_3,"vel_3_std.dat");
    	write_datafile(et3,"e_dis_AAPG3.dat");
        fprintf(err_dump,"%lg %lg\n",e3(0),e3(1)); 
        e_AAPG.replaceRow(e3,2);
    }
    return(t);
}

void PostProcess(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, int AAPG_ord, Array1D<double>& sol_0, Array1D<Array2D<double> >& sol_1, Array2D<Array2D<double> >& sol_2, Array3D<Array2D<double> >& sol_3, Array1D<double>& sol_1_mean, Array1D<double>& sol_2_mean, Array1D<double>& sol_3_mean, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3, Array1D<double>& normsq, int dim, int nStep, int PCTerms_1, int PCTerms_2, int PCTerms_3, int order, double dTym, double factor_OD, Array2D<double>& mstd_MCS, Array2D<double>& samPts_norm, string name,int noutput, Array1D<Array1D<double> >& e_sample, bool PDF, string pcType, Array2D<double>& sample_mstd_2D){
    TickTock tt;
    tt.tick();
    // Post-process the AAPG solutions
    // initialization
    int nAAPGTerms = 1+dim+dim*(dim-1)/2+dim*(dim-1)*(dim-2)/6;
    int  nPCTerms = computeNPCTerms(dim, order);
    //int nPCTerms = myPCSet.GetNumberPCTerms();
    Array2D<double> sol_1_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> sol_2_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> sol_3_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> coeffAAPG1(1+dim,nStep+1,1.0);
    Array2D<double> coeffAAPG2(1+dim+dim*(dim-1)/2,nStep+1,1.0);
    Array2D<double> coeffAAPG3(nAAPGTerms,nStep+1,1.0);

    tt.tock("Initialization Took");
 
    // assemble and save the mean values
    printf("Assembling the mean...\n");
    tt.tick();
    assemblemean(indi_2, indj_2, indi_3, indj_3, indk_3, sol_0, sol_1, sol_2, sol_3, dim, nStep, AAPG_ord, sol_1_mean, sol_2_mean, sol_3_mean, coeffAAPG1, coeffAAPG2, coeffAAPG3);
    tt.tock("Assemble mean took");   
 
    // add the mean to the assembled solution
    sol_1_assembled.replaceCol(sol_1_mean,0);
    sol_2_assembled.replaceCol(sol_2_mean,0); 
    sol_3_assembled.replaceCol(sol_3_mean,0); 

    // assemble the rest PC terms
    printf("Assembling the rest...\n");
    tt.tick();
    assemblerest(indi_2, indj_2, indi_3, indj_3, indk_3, sol_1, sol_2, sol_3, sol_1_assembled, sol_2_assembled, sol_3_assembled, PCTerms_1, PCTerms_2, PCTerms_3, dim, order, nStep, AAPG_ord, coeffAAPG1, coeffAAPG2, coeffAAPG3);
    tt.tock("Assemble rest took");    

    // compute and print the std values
    printf("Computing the std...\n");
    tt.tick();
    computeStd(nStep, nPCTerms, sol_1_assembled,sol_2_assembled, sol_3_assembled, normsq, std1, std2, std3);
    tt.tock("Compute std took");

    // output sol_1_assembled for debug
    //ostringstream s0;
    //s0 << name << "_1_assembled"<<".dat";
    //write_datafile(sol_1_assembled,s0.str().c_str());
    // output sol_2_assembled for debug
    //ostringstream s4;
    //s4 << name << "_2_assembled"<<".dat";
    //write_datafile(sol_2_assembled,s4.str().c_str());

    // sample result
    Array2D<double> stat1(2,nStep+1,0.e0);
    Array2D<double> stat2(2,nStep+1,0.e0);
    Array2D<double> stat3(2,nStep+1,0.e0);
    Array1D<double> temp_stat(nStep+1,0.e0);
    stat1.replaceRow(sol_1_mean,0);
    getCol(sol_2_assembled,0,temp_stat);
    stat2.replaceRow(temp_stat,0);
    getCol(sol_3_assembled,0,temp_stat);
    stat3.replaceRow(temp_stat,0);
    stat1.replaceRow(std1,1);
    stat2.replaceRow(std2,1);
    stat3.replaceRow(std3,1);

    ostringstream s6;
    s6 << name << "stat1"<<".dat";
    write_datafile(stat1,s6.str().c_str());
    ostringstream s5;
    s5 << name << "stat2"<<".dat";
    write_datafile(stat2,s5.str().c_str());

    if (PDF){
        PCSet myPCSet("ISP",order,dim,pcType,0.0,1.0); 
        const int nPCTerms = myPCSet.GetNumberPCTerms();
        Array1D<Array1D<double> > initial_GS(2);
        Array1D<double> temp(nPCTerms,0.e0);
        initial_GS(0)=temp;
        initial_GS(1) = temp;
        myPCSet.InitMeanStDv(sample_mstd_2D(0,0),sample_mstd_2D(0,1),1,initial_GS(0));
        myPCSet.InitMeanStDv(sample_mstd_2D(1,0),sample_mstd_2D(1,1),2,initial_GS(1));
        Array2D<double> AAPG_sol_sample_1=sampleGS(noutput,dim, nStep, nPCTerms, myPCSet, sol_1_assembled, samPts_norm, stat1, e_sample(0));
        ostringstream s;
        s << "AAPG" << name << "sample_1"<<".dat";
        write_datafile(AAPG_sol_sample_1,s.str().c_str());
        if (AAPG_ord >= 2){
            Array2D<double> AAPG_sol_sample_2=sampleGS(noutput,dim, nStep, nPCTerms, myPCSet, sol_2_assembled, samPts_norm, stat2, e_sample(1));
            ostringstream s2;
            s2 << "AAPG" << name << "sample_2"<<".dat";
            write_datafile(AAPG_sol_sample_2,s2.str().c_str());
        }
        if (AAPG_ord >= 3){
            Array2D<double> AAPG_sol_sample_3=sampleGS(noutput,dim, nStep, nPCTerms, myPCSet, sol_3_assembled, samPts_norm, stat3, e_sample(2));
            ostringstream s3;
            s3 << "AAPG" << name << "sample_3"<<".dat";
            write_datafile(AAPG_sol_sample_3,s3.str().c_str());
        }
    }

    return;     
}

Array1D<int> identicalrow(int PCTerms, int dim, Array2D<int>& Pbtot, Array2D<int>& Pbmore){
    // find out identical row
    Array1D<int> ind(PCTerms,0);
    int ii = 0;
    int jj = 0;

    while (ii < PCTerms){
        int kk=0;
        while (kk<dim&&Pbtot(jj,kk)==Pbmore(ii,kk))
            kk++;
        if (kk>dim-2){
            ind(ii) = jj;
            ii++;
        }
        jj++;
    }
    return (ind);
}

void index(int dim, int order, int PCTerms_1, int PCTerms_2, int PCTerms_3, Array2D<int>& ind1, Array2D<Array1D<int> >& ind2, Array3D<Array1D<int> >& ind3){
    // Get corresponding row index in the global matrix for sub-problem solutions
    int Px = computeNPCTerms(dim, order);
    Array2D<int> Pbtot(Px,dim);
    computeMultiIndex(dim,order,Pbtot);
    write_datafile(Pbtot,"indexglobal.dat");
    // 1D psibasis
    Array2D<int> Pb1(PCTerms_1,1);
    computeMultiIndex(1,order,Pb1);
    // 2D psibasis
    Array2D<int> Pb2(PCTerms_2,2);
    computeMultiIndex(2,order,Pb2);
    // 3D psibasis
    Array2D<int> Pb3(PCTerms_3,3);
    computeMultiIndex(3,order,Pb3);
   
    for (int i=0;i<dim;i++){
        Array2D<int> Pb1_more(PCTerms_1,dim,0);
        for (int ind=0;ind<PCTerms_1;ind++){
            Pb1_more(ind,i)=Pb1(ind,0);
        }
        Array1D<int> ind=identicalrow(PCTerms_1,dim,Pbtot,Pb1_more);
        ind1.replaceCol(ind,i);
    }
    
    for (int i=0;i<dim-1;i++){
        for (int j=i+1;j<dim;j++){
            Array2D<int> Pb2_more(PCTerms_2,dim,0);
            for (int ind=0;ind<PCTerms_2;ind++){
                Pb2_more(ind,i)=Pb2(ind,0);
                Pb2_more(ind,j)=Pb2(ind,1);
            }
            ind2(i,j)=identicalrow(PCTerms_2,dim,Pbtot,Pb2_more);
        }
    }

    for (int i=0;i<dim-2;i++){
        for (int j=i+1;j<dim-1;j++){
            for (int k=j+1;k<dim;k++){
                Array2D<int> Pb3_more(PCTerms_3,dim,0);
                for (int ind=0;ind<PCTerms_3;ind++){
                    Pb3_more(ind,i)=Pb3(ind,0);    
                    Pb3_more(ind,j)=Pb3(ind,1);    
                    Pb3_more(ind,k)=Pb3(ind,2);    
                }
                ind3(i,j,k)=identicalrow(PCTerms_3,dim,Pbtot,Pb3_more);
            }
        }
    }
}

void assemblemean(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<double>& sol_0, Array1D<Array2D<double> >& sol_1, Array2D<Array2D<double> >& sol_2, Array3D<Array2D<double> >& sol_3, int dim, int nStep, int AAPG_ord, Array1D<double>& sol_1_mean, Array1D<double>& sol_2_mean, Array1D<double>& sol_3_mean, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3){
    // save first-order AAPG component functions for the assembly of second-order AAPG mean value
    //Array1D<Array1D<double> > sol_1_mean_ind(dim);
    Array2D<double> sol_1_mean_ind(dim, nStep+1, 0.e0);
    if (AAPG_ord >= 1){
        // First-order terms
        for (int i=0;i<dim;i++){
            // Retrieve mean
            Array1D<double> Temp1(nStep+1,0.e0);
            getCol(sol_1(i),0,Temp1);
            // Compute correction terms on the mean
            subtractVec(sol_0,Temp1);
            // link Temp1 to sol_1_mean_ind
            //sol_1_mean_ind(i) = Temp1;        
            sol_1_mean_ind.replaceRow(Temp1,i);
            // Add correction terms to the mean
            Array1D<double> Temp2(nStep+1,0.e0);
            getRow(coeff1,i+1,Temp2);
            dotprodVec(Temp2,Temp1);
            addVec(Temp1,sol_1_mean);
        }
    }
    // second-order terms
    int n = dim+1;
    sol_2_mean = sol_1_mean;
    Array2D<Array1D<double> > sol_2_mean_ind(dim, dim); 
    if (AAPG_ord >= 2){
        for (unsigned int i=0;i<indi_2.XSize();i++){
                Array1D<double> Temp(nStep+1,0.e0);
                getCol(sol_2(indi_2(i),indj_2(i)),0,Temp);
                // subtract the zero-th order term
                subtractVec(sol_0,Temp); 
                // subtract the first order term
                Array1D<double> Temp1(nStep+1,0.e0);
                getRow(sol_1_mean_ind,indi_2(i),Temp1);
                subtractVec(Temp1,Temp);
                getRow(sol_1_mean_ind,indj_2(i),Temp1);
                subtractVec(Temp1,Temp);
                // save the mean value for third-order use
                sol_2_mean_ind(indi_2(i),indj_2(i))=Temp;
                // Add correction terms to the mean
                Array1D<double> Temp2(nStep+1,0.e0);
                getRow(coeff2,n,Temp2);
                dotprodVec(Temp2,Temp);
                addVec(Temp,sol_2_mean);
                n++;
        }
    }
    // third-order terms
    sol_3_mean = sol_2_mean;
    if (AAPG_ord >= 3){
        for (unsigned ii=0;ii<indi_3.XSize();ii++){
                    Array1D<double> Temp(nStep+1,0.e0);
                    getCol(sol_3(indi_3(ii),indj_3(ii),indk_3(ii)),0,Temp);
                    // subtract the zero-th order terms
                    subtractVec(sol_0,Temp); 
                    // subtract the first order terms
                    Array1D<double> Temp1(nStep+1,0.e0);
                    getRow(sol_1_mean_ind,indi_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    getRow(sol_1_mean_ind,indj_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    getRow(sol_1_mean_ind,indk_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    // subtract the second order terms
                    subtractVec(sol_2_mean_ind(indi_3(ii),indj_3(ii)),Temp);
                    subtractVec(sol_2_mean_ind(indj_3(ii),indk_3(ii)),Temp);
                    subtractVec(sol_2_mean_ind(indi_3(ii),indk_3(ii)),Temp);
                    // Add correction terms to the mean
                    Array1D<double> Temp2(nStep+1,0.e0);
                    getRow(coeff3,n,Temp2);
                    dotprodVec(Temp2,Temp);
                    addVec(Temp,sol_3_mean);
                    n++;
        }
    }
    return;
}

void assemblerest(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<Array2D<double> >& sol_1, Array2D<Array2D<double> >& sol_2, Array3D<Array2D<double> >& sol_3, Array2D<double>& sol_1_assembled, Array2D<double>& sol_2_assembled, Array2D<double>& sol_3_assembled, int PCTerms_1, int PCTerms_2, int PCTerms_3, int dim, int order, int nStep, int AAPG_ord, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3){

    printf("Computing the index...\n");
    Array2D<int> ind1(PCTerms_1,dim,0);
    Array2D<Array1D<int> > ind2(dim,dim);
    Array3D<Array1D<int> > ind3(dim,dim,dim);
    index(dim, order, PCTerms_1, PCTerms_2, PCTerms_3, ind1, ind2, ind3);

    // Assemble first order AAPG solutions
    printf("Assembling the first order terms...\n");
    for (int i=0;i<dim;i++){
        for (int j=1;j<PCTerms_1;j++){
            // retrieve the col-index on the global PC matrix
            int ind = ind1(j,i);
            for (int it=0;it<nStep+1;it++){
                sol_1_assembled(it,ind)=sol_1_assembled(it,ind)+coeff1(i+1,it)*sol_1(i)(it,j);
            }
        }        
    }
        
    // Assemble second order AAPG solutions
    printf("Assembling the second order terms...\n");
    sol_2_assembled = sol_1_assembled;
    int n = dim+1;
    Array1D<double> coeffn(nStep+1,0.e0);
    if (AAPG_ord >= 2){
        for (unsigned int i=0;i<indi_2.XSize();i++){
                getRow(coeff2,n,coeffn);
                // Add second-order AAPG terms
                for (int p=1;p<PCTerms_2;p++){
                    // retrieve the index in the global matrix
                    int ind = ind2(indi_2(i),indj_2(i))(p);
                    for (int it=0;it<nStep+1;it++){
                        sol_2_assembled(it,ind)=sol_2_assembled(it,ind)+coeffn(it)*sol_2(indi_2(i),indj_2(i))(it,p);
                     }
                }

                // Subtract first-order AAPG terms
                for (int p=1;p<PCTerms_1;p++){
                    subtractfirst(p,indi_2(i),ind1,sol_1,nStep,coeffn,sol_2_assembled);
                    subtractfirst(p,indj_2(i),ind1,sol_1,nStep,coeffn,sol_2_assembled);
                }
                n++;    
        }
    }    
    
    // Assemble third order AAPG solutions
    sol_3_assembled = sol_2_assembled;
    if (AAPG_ord >= 3){
        printf("Assembling the third order terms...\n");
        for (unsigned int ii=0;ii<indi_3.XSize();ii++){
                    // retrieve the index in the global matrix
                    Array1D<int> ind = ind3(indi_3(ii),indj_3(ii),indk_3(ii));
                    Array1D<double> coeffn(nStep+1,0.e0);
                    getRow(coeff3,n,coeffn);
                    
                    // Add third-order AAPG terms
                    for (int p=1;p<PCTerms_3;p++){
                        for (int it=0;it<nStep+1;it++)
                            sol_3_assembled(it,ind(p))=sol_3_assembled(it,ind(p))+coeffn(it)*sol_3(indi_3(ii),indj_3(ii),indk_3(ii))(it,p);
                    }

                    // Subtract second-order AAPG terms
                    for (int p=1;p<PCTerms_2;p++){
                        subtractsecond(p,indi_3(ii),indj_3(ii),ind2,sol_2,nStep,coeffn,sol_3_assembled); 
                        subtractsecond(p,indj_3(ii),indk_3(ii),ind2,sol_2,nStep,coeffn,sol_3_assembled); 
                        subtractsecond(p,indi_3(ii),indk_3(ii),ind2,sol_2,nStep,coeffn,sol_3_assembled); 
                    }

                    // Subtract first-order AAPG terms
                    coeffn.Resize(nStep+1,-1); 
                    for (int p=1;p<PCTerms_1;p++){
                        subtractfirst(p,indi_3(ii),ind1,sol_1,nStep,coeffn,sol_3_assembled);
                        subtractfirst(p,indj_3(ii),ind1,sol_1,nStep,coeffn,sol_3_assembled);
                        subtractfirst(p,indk_3(ii),ind1,sol_1,nStep,coeffn,sol_3_assembled);
                    }
                    n++;
        }
    }
    return;
}

void computeStd(int nStep, int nPCTerms, Array2D<double>& sol_1_assembled, Array2D<double>& sol_2_assembled, Array2D<double>& sol_3_assembled, Array1D<double>& normsq, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3){
    // compute std and save
    for (int i=0;i<nStep+1;i++){
        double temp1 = 0.e0;
        double temp2 = 0.e0;
        double temp3 = 0.e0;
        for (int j=1;j<nPCTerms;j++){
            temp1 += sol_1_assembled(i,j)*sol_1_assembled(i,j)*normsq(j);
            temp2 += sol_2_assembled(i,j)*sol_2_assembled(i,j)*normsq(j);
            temp3 += sol_3_assembled(i,j)*sol_3_assembled(i,j)*normsq(j);
        }
        std1(i) = sqrt(temp1);
        std2(i) = sqrt(temp2);
        std3(i) = sqrt(temp3);
    }
    //write_datafile_1d(std1,"AAPG1_std.dat");
    //write_datafile_1d(std2,"AAPG2_std.dat");
    //write_datafile_1d(std3,"AAPG3_std.dat");
    return;
}

void subtractfirst(int p,int ii, Array2D<int>& ind1,Array1D<Array2D<double> >& sol_1,int nStep,Array1D<double>& coeffn,Array2D<double>& sol_assembled){
    int ind = ind1(p,ii);
    for (int it=0;it<nStep+1;it++){
        sol_assembled(it,ind) = sol_assembled(it,ind)-coeffn(it)*sol_1(ii)(it,p);
    }

    return;
}

void subtractsecond(int p,int ii,int jj, Array2D<Array1D<int> >& ind2, Array2D<Array2D<double> >& sol_2, int nStep, Array1D<double>& coeffn, Array2D<double>& sol_assembled){
    int ind = ind2(ii,jj)(p);
    for (int it=0;it<nStep+1;it++){
        sol_assembled(it,ind) = sol_assembled(it,ind)-coeffn(it)*sol_2(ii,jj)(it,p);
    }
    
    return;
}
