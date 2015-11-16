#include <math.h>
#include "Array2D.h"
#include "Array1D.h"
#include "PCSet.h"
#include "arraytools.h"
#include "uqtktools.h"
#include "UtilsDuffing.h"
#include "AAPG.h"
#include "MCS.h"
#include "lapack.h"
#include "GhanemSpanos.h"
#include "ticktock.h"

Array1D<double> AAPG(Array1D<double> inpParams, double fbar, double dTym, int order, string pcType, int dim, int nStep, Array2D<double>& scaledKLmodes, double dis0, double vel0, PCSet& myPCSet, double factor_OD, int AAPG_ord, bool act_D){
    // timing var
    Array1D<double> t(5,0.e0);
    
    // Compute zeroth-order term
    printf("Zeroth-order term...\n");
    // variable to store solution at every step. Zero initial condition.
    Array1D<double> Temp(2,0.e0);
    // deterministic force
    Array1D<double> f_0(2,0.e0);
    f_0(0)=fbar;
    // store zeroth order solution
    Array2D<double> x_0(nStep+1,2,0.e0);
    // Time marching steps
    TickTock tt;
    tt.tick();
    for (int ix=0;ix<nStep;ix++){
        forward_duffing_dt(inpParams,f_0,dTym,Temp);
        x_0.replaceRow(Temp,ix+1);
    }
    tt.tock("Took");
    t(0) = tt.silent_tock();
    // abstract the displacement terms
    Array1D<double> dis_0(nStep+1,0.e0);
    getCol(x_0,1,dis_0);
    
    // Compute first order terms
    printf("First-order terms...");
    // initialize the displacement and force terms
    Array1D<Array2D<double> > dis_1(dim);
    Array1D<Array2D<double> > force_1(dim);
    // generate PCSet and the number of terms in it
    PCSet PCSet_1("ISP",order,1,pcType,0.0,1.0); 
    int PCTerms_1 = PCSet_1.GetNumberPCTerms();
    // Generate the forcing on each dim
    for (int i=0;i<dim;i++){
        Array2D<double> f_1(nStep+1,PCTerms_1,0.e0);
        for (int it=0;it<nStep+1;it++)
            f_1(it,0) = fbar;
        // abstract the corresponding forcing term and append it to f_1
        Array1D<double> tempf(nStep+1,0.e0);
        getCol(scaledKLmodes,i,tempf);
        f_1.replaceCol(tempf,1);
        force_1(i) = f_1;
    }
    
    tt.tick();
    #pragma omp parallel default(none) shared(dim,PCSet_1,order,PCTerms_1,pcType,nStep,dis0,vel0,dTym,inpParams,force_1,dis_1)
    {
    #pragma omp for
    for (int i=0;i<dim;i++){
        dis_1(i) = GS(PCSet_1, order, 1, PCTerms_1, pcType, nStep, dis0, vel0, dTym, inpParams, force_1(i));
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
    double p = 0.99; //ratio of the sum(var) on active dims 
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
    Array2D<Array2D<double> > dis_2(dim,dim); 
    int PCTerms_2 = 0;
    Array1D<int> indi_2(N_adof*(N_adof-1)/2,0);
    Array1D<int> indj_2(N_adof*(N_adof-1)/2,0);
    if (AAPG_ord >= 2){
    printf("Second-order terms...");
    PCSet PCSet_2("ISP",order,2,pcType,0.0,1.0); 
    PCTerms_2 = PCSet_2.GetNumberPCTerms();
    Array2D<double> f_2(nStep+1,PCTerms_2,0.e0);
    for (int it=0;it<nStep+1;it++)
        f_2(it,0) = fbar;
    Array1D<Array2D<double> > force_2(N_adof*(N_adof-1)/2);
    int k = 0;
    for (int i=0;i<N_adof-1;i++){
        for (int j=i+1;j<N_adof;j++){
            Array1D<double> tempf(nStep+1,0.e0);
            getCol(scaledKLmodes,ind(i),tempf);
            f_2.replaceCol(tempf,1);
            getCol(scaledKLmodes,ind(j),tempf);
            f_2.replaceCol(tempf,2);
	        force_2(k)=f_2;
	        indi_2(k) = ind(i);
    	    indj_2(k) = ind(j);
	        k++;
	}
    }
    tt.tick();
    #pragma omp parallel  default(none) shared(k,dis_2,indi_2,indj_2,PCSet_2,order,PCTerms_2,pcType,nStep,dis0,vel0,dTym,inpParams,force_2)
    {
    #pragma omp for
    for (int i=0;i<k;i++){ 
        dis_2(indi_2(i),indj_2(i)) = GS(PCSet_2, order, 2, PCTerms_2, pcType, nStep, dis0, vel0, dTym, inpParams, force_2(i));
    }
    }
    tt.tock("Took");   
    t(2)=tt.silent_tock();
    cout << "Saving on computational cost using only active dims: " << (1-(N_adof*(N_adof-1)*1.0/dim/(dim-1)))*100 << "%" << endl;
    }
 
    // Third order term
    Array3D<Array2D<double> > dis_3(dim,dim,dim); 
    int PCTerms_3 = 0;
    Array1D<int> indi_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indj_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    Array1D<int> indk_3(N_adof*(N_adof-1)*(N_adof-2)/6,0);
    if (AAPG_ord >= 3){
    printf("Third-order terms...");
    PCSet PCSet_3("ISP",order,3,pcType,0.0,1.0); 
    PCTerms_3 = PCSet_3.GetNumberPCTerms();
    Array2D<double> f_3(nStep+1,PCTerms_3,0.e0);
    Array1D<Array2D<double> > force_3(N_adof*(N_adof-1)*(N_adof-2)/6);
    for (int ix=0;ix<nStep+1;ix++)
        f_3(ix,0) = fbar;
    int l=0;
    for (int i=0;i<N_adof-2;i++){
        for (int j=i+1;j<N_adof-1;j++){
            for (int k=j+1;k<N_adof;k++){
                Array1D<double> tempf(nStep+1,0.e0);
                getCol(scaledKLmodes,ind(i),tempf);
                f_3.replaceCol(tempf,1);
                getCol(scaledKLmodes,ind(j),tempf);
                f_3.replaceCol(tempf,2);
                getCol(scaledKLmodes,ind(k),tempf);
                f_3.replaceCol(tempf,3);
		        force_3(l)=f_3;
    		    indi_3(l)=ind(i);
		        indj_3(l)=ind(j);
		        indk_3(l)=ind(k);
		        l++;
                }
            }
	    }
    //start = clock();
    tt.tick();
    #pragma omp parallel default(none) shared(l,indi_3,indj_3,indk_3,PCSet_3,order,PCTerms_3,pcType,nStep,dis0,vel0,dTym,inpParams,dis_3,force_3)
    {
    #pragma omp for
    for (int i=0;i<l;i++)        
	dis_3(indi_3(i),indj_3(i),indk_3(i)) = GS(PCSet_3, order, 3, PCTerms_3, pcType, nStep, dis0, vel0, dTym, inpParams, force_3(i));
    }
    tt.tock("Took");
    t(3)=tt.silent_tock();
    cout << "Saving on computational cost using only active dims: " << (1-(N_adof*(N_adof-1)*(N_adof-2)*1.0/dim/(dim-1)/(dim-2)))*100 << "%" << endl;
    }
 
    printf("\nAssemble the solutions...\n");
    tt.tick();
    PostProcess(indi_2,indj_2, indi_3, indj_3, indk_3, AAPG_ord, dis_0, dis_1, dis_2, dis_3, myPCSet, fbar, dim, nStep, PCTerms_1, PCTerms_2, PCTerms_3, order, dTym, pcType, inpParams, scaledKLmodes, factor_OD);
    tt.tock("Took");
    t(4)=tt.silent_tock();
   
    return(t);
}

void PostProcess(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, int AAPG_ord, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, PCSet& myPCSet, double fbar, int dim, int nStep, int PCTerms_1, int PCTerms_2, int PCTerms_3, int order, double dTym, string pcType, Array1D<double>& inpParams, Array2D<double>& scaledKLmodes, double factor_OD){
    TickTock tt;
    tt.tick();
    // Post-process the AAPG solutions
    // initialization
    int nAAPGTerms = 1+dim+dim*(dim-1)/2+dim*(dim-1)*(dim-2)/6;
    int nPCTerms = myPCSet.GetNumberPCTerms();
    Array2D<double> dis_1_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> dis_2_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> dis_3_assembled(nStep+1,nPCTerms,0.e0);
    Array2D<double> coeffAAPG1(1+dim,nStep+1,1.0);
    Array2D<double> coeffAAPG2(1+dim+dim*(dim-1)/2,nStep+1,1.0);
    Array2D<double> coeffAAPG3(nAAPGTerms,nStep+1,1.0);

    // initialize the first/second order mean/std of the AAPG solutions
    Array1D<double> dis_1_mean = dis_0;
    Array1D<double> dis_2_mean = dis_0;
    Array1D<double> dis_3_mean = dis_0;
   
    tt.tock("Initialization Took");
 
    // assemble and save the mean values
    printf("Assembling the mean...\n");
    tt.tick();
    assemblemean(indi_2, indj_2, indi_3, indj_3, indk_3, dis_0, dis_1, dis_2, dis_3, dim, nStep, AAPG_ord, dis_1_mean, dis_2_mean, dis_3_mean, coeffAAPG1, coeffAAPG2, coeffAAPG3);
    tt.tock("Assemble mean took");   
 
    // add the mean to the assembled solution
    dis_1_assembled.replaceCol(dis_1_mean,0);
    dis_2_assembled.replaceCol(dis_2_mean,0); 
    dis_3_assembled.replaceCol(dis_3_mean,0); 

    // assemble the rest PC terms
    printf("Assembling the rest...\n");
    tt.tick();
    assemblerest(indi_2, indj_2, indi_3, indj_3, indk_3, dis_1, dis_2, dis_3, dis_1_assembled, dis_2_assembled, dis_3_assembled, PCTerms_1, PCTerms_2, PCTerms_3, dim, order, nStep, AAPG_ord, coeffAAPG1, coeffAAPG2, coeffAAPG3);
    tt.tock("Assemble rest took");    

    // compute and print the std values
    Array1D<double> std1(nStep+1,0.e0);
    Array1D<double> std2(nStep+1,0.e0);
    Array1D<double> std3(nStep+1,0.e0);
    printf("Computing the std...\n");
    tt.tick();
    computeStd(nStep, nPCTerms, dis_1_assembled,dis_2_assembled, dis_3_assembled, myPCSet, std1, std2, std3);
    tt.tock("Compute std took");
    
    // print out mean/std valus at specific points for comparison
    printf("First-order AAPG results:\n");
    if (AAPG_ord >= 1){
    	for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/10) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_1_mean(ix),std1(ix));
            }
    	}
    	write_datafile_1d(dis_1_mean,"dis_1_mean.dat");
    	write_datafile_1d(std1,"dis_1_std.dat");
    }
    if(AAPG_ord >= 2){
	printf("Second-order AAPG results:\n");
    	for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/10) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_2_mean(ix),std2(ix));
            }
    	}
    	write_datafile_1d(dis_2_mean,"dis_2_mean.dat");
    	write_datafile_1d(std2,"dis_2_std.dat");
    }
    if(AAPG_ord >= 3){
        printf("Third-order AAPG results:\n");
        for (int ix=0;ix<nStep+1;ix++){
            if (ix % ((int) nStep/10) == 0){
            	WriteMeanStdDevToStdOut(ix,ix*dTym,dis_3_mean(ix),std3(ix));
        	}
    	}
    	write_datafile_1d(dis_3_mean,"dis_3_mean.dat");
    	write_datafile_1d(std3,"dis_3_std.dat");
    }
    
    return;     
}

Array1D<int> index1(int dim, int i, int order){
    // full psibasis
    
    int Px = computeNPCTerms(dim, order);
    Array2D<int> Pbtot(Px,dim);
    computeMultiIndex(dim,order,Pbtot);

    // 1D psibasis
    int Px1 = computeNPCTerms(1, order);
    Array2D<int> Pb1(Px1,1);
    //int temp2 = computeMultiIndex(1,order,Pb1);
    computeMultiIndex(1,order,Pb1);
    Array1D<int> Pb1_1D(Px1,1);
    array2Dto1D(Pb1,Pb1_1D);
    // expand to full psibasis size
    Array2D<int> Pb1_more(Px1,dim,0);
    Pb1_more.replaceCol(Pb1_1D,i);
    
    // find out identical row
    Array1D<int> index(Px1,0);
    int ii = 1;
    int jj = 1;

    while (ii <= Px1){
        int test = 0;
        for (int kk=1;kk<dim+1;kk++){
            test += abs(Pbtot(jj-1,kk-1)-Pb1_more(ii-1,kk-1));    
        }    
        if (test == 0){
            index(ii-1) = jj-1;
            ii++;    
        }
        jj++;
    }

    return (index);    
}

Array1D<int> index2(int dim, int i, int j, int order){
    int Px = computeNPCTerms(dim, order);
    Array2D<int> Pbtot(Px,dim);
    computeMultiIndex(dim,order,Pbtot);

    // 2D psibasis
    int Px2 = computeNPCTerms(2, order);
    Array2D<int> Pb2(Px2,2);
    computeMultiIndex(2,order,Pb2);
    // expand to full psibasis size
    Array2D<int> Pb2_more(Px2,dim,0);
    Array1D<int> temp3(Px2,0.e0);
    getCol(Pb2,0,temp3);
    Pb2_more.replaceCol(temp3,i);
    getCol(Pb2,1,temp3);
    Pb2_more.replaceCol(temp3,j);
    
    // find out identical row
    Array1D<int> index(Px2,0);
    int ii = 1;
    int jj = 1;

    while (ii <= Px2){
        int test = 0;
        for (int kk=1;kk<dim+1;kk++){
            test += abs(Pbtot(jj-1,kk-1)-Pb2_more(ii-1,kk-1));    
        }    
        if (test == 0){
            index(ii-1) = jj-1;
            ii++;    
        }
        jj++;
    }

    return (index);    
}

Array1D<int> index3(int dim, int i, int j, int k, int order){
    int Px = computeNPCTerms(dim, order);
    Array2D<int> Pbtot(Px,dim);
    computeMultiIndex(dim,order,Pbtot);

    // 3D psibasis
    int Px3 = computeNPCTerms(3, order);
    Array2D<int> Pb3(Px3,3);
    computeMultiIndex(3,order,Pb3);
    // expand to full psibasis size
    Array2D<int> Pb3_more(Px3,dim,0);
    Array1D<int> temp3(Px3,0.e0);
    getCol(Pb3,0,temp3);
    Pb3_more.replaceCol(temp3,i);
    getCol(Pb3,1,temp3);
    Pb3_more.replaceCol(temp3,j);
    getCol(Pb3,2,temp3);
    Pb3_more.replaceCol(temp3,k);
    
    // find out identical row
    Array1D<int> index(Px3,0);
    int ii = 1;
    int jj = 1;

    while (ii <= Px3){
        int test = 0;
        for (int tt=1;tt<dim+1;tt++){
            test += abs(Pbtot(jj-1,tt-1)-Pb3_more(ii-1,tt-1));    
        }    
        if (test == 0){
            index(ii-1) = jj-1;
            ii++;    
        }
        jj++;
    }

    return (index);    
}

void assemblemean(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, int dim, int nStep, int AAPG_ord, Array1D<double>& dis_1_mean, Array1D<double>& dis_2_mean, Array1D<double>& dis_3_mean, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3){
    // save first-order AAPG component functions for the assembly of second-order AAPG mean value
    //Array1D<Array1D<double> > dis_1_mean_ind(dim);
    Array2D<double> dis_1_mean_ind(dim, nStep+1, 0.e0);
    if (AAPG_ord >= 1){
        // First-order terms
        for (int i=0;i<dim;i++){
            // Retrieve mean
            Array1D<double> Temp1(nStep+1,0.e0);
            getCol(dis_1(i),0,Temp1);
            // Compute correction terms on the mean
            subtractVec(dis_0,Temp1);
            // link Temp1 to dis_1_mean_ind
            //dis_1_mean_ind(i) = Temp1;        
            dis_1_mean_ind.replaceRow(Temp1,i);
            // Add correction terms to the mean
            Array1D<double> Temp2(nStep+1,0.e0);
            getRow(coeff1,i+1,Temp2);
            dotprodVec(Temp2,Temp1);
            addVec(Temp1,dis_1_mean);
        }
    }
    // second-order terms
    int n = dim+1;
    dis_2_mean = dis_1_mean;
    Array2D<Array1D<double> > dis_2_mean_ind(dim, dim); 
    if (AAPG_ord >= 2){
        for (unsigned int i=0;i<indi_2.XSize();i++){
                Array1D<double> Temp(nStep+1,0.e0);
                getCol(dis_2(indi_2(i),indj_2(i)),0,Temp);
                // subtract the zero-th order term
                subtractVec(dis_0,Temp); 
                // subtract the first order term
                Array1D<double> Temp1(nStep+1,0.e0);
                getRow(dis_1_mean_ind,indi_2(i),Temp1);
                subtractVec(Temp1,Temp);
                getRow(dis_1_mean_ind,indj_2(i),Temp1);
                subtractVec(Temp1,Temp);
                // save the mean value for third-order use
                dis_2_mean_ind(indi_2(i),indj_2(i))=Temp;
                // Add correction terms to the mean
                Array1D<double> Temp2(nStep+1,0.e0);
                getRow(coeff2,n,Temp2);
                dotprodVec(Temp2,Temp);
                addVec(Temp,dis_2_mean);
                n++;
        }
    }
    // third-order terms
    dis_3_mean = dis_2_mean;
    if (AAPG_ord >= 3){
        for (unsigned ii=0;ii<indi_3.XSize();ii++){
                    Array1D<double> Temp(nStep+1,0.e0);
                    getCol(dis_3(indi_3(ii),indj_3(ii),indk_3(ii)),0,Temp);
                    // subtract the zero-th order terms
                    subtractVec(dis_0,Temp); 
                    // subtract the first order terms
                    Array1D<double> Temp1(nStep+1,0.e0);
                    getRow(dis_1_mean_ind,indi_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    getRow(dis_1_mean_ind,indj_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    getRow(dis_1_mean_ind,indk_3(ii),Temp1);
                    subtractVec(Temp1,Temp);
                    // subtract the second order terms
                    subtractVec(dis_2_mean_ind(indi_3(ii),indj_3(ii)),Temp);
                    subtractVec(dis_2_mean_ind(indj_3(ii),indk_3(ii)),Temp);
                    subtractVec(dis_2_mean_ind(indi_3(ii),indk_3(ii)),Temp);
                    // Add correction terms to the mean
                    Array1D<double> Temp2(nStep+1,0.e0);
                    getRow(coeff3,n,Temp2);
                    dotprodVec(Temp2,Temp);
                    addVec(Temp,dis_3_mean);
                    n++;
        }
    }
    return;
}

void assemblerest(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, int PCTerms_1, int PCTerms_2, int PCTerms_3, int dim, int order, int nStep, int AAPG_ord, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3){
    // Get corresponding row index in the global matrix for sub-problem solutions
    printf("Computing the index...\n");
    Array2D<int> ind1(PCTerms_1,dim,0);
    for (int i=0;i<dim;i++){
        Array1D<int> ind = index1(dim, i, order);
        ind1.replaceCol(ind,i);
    }
    Array3D<int> ind2(dim,dim,PCTerms_2,0);
    for (int i=0;i<dim-1;i++){
        for (int j=i+1;j<dim;j++){
            Array1D<int> ind = index2(dim, i, j, order);
            for (int k=0;k<PCTerms_2;k++)
                ind2(i,j,k) = ind(k);
        }
    }
    Array3D<Array1D<int> > ind3(dim,dim,dim);
    for (int i=0;i<dim-2;i++){
        for (int j=i+1;j<dim-1;j++){
            for (int k=j+1;k<dim;k++){
                ind3(i,j,k) = index3(dim, i, j, k, order);
            }
        }
    }

    // Assemble first order AAPG solutions
    printf("Assembling the first order terms...\n");
    Array1D<double> coeffn(nStep+1,0.e0);
    for (int i=0;i<dim;i++){
        // retrieve the solution on dim i
        //Array2D<double> Temp(nStep+1,PCTerms_1,0.e0);
        //get2D(dis_1,i,Temp);
        getRow(coeff1,i,coeffn);
        for (int j=1;j<PCTerms_1;j++){
            // retrieve the col-index on the global PC matrix
            int ind = ind1(j,i);
            // retrieve the j-th term of the solution
            Array1D<double> Temp1(nStep+1,0.e0);
            getCol(dis_1(i),j,Temp1);
            // add the retrieved j-th term to the corresponding col
            Array1D<double> Temp2(nStep+1,0.e0);
            getCol(dis_1_assembled,ind,Temp2);
            dotprodVec(coeffn,Temp1);
            addVec(Temp1,Temp2);
            dis_1_assembled.replaceCol(Temp2,ind); 
            //for (int it=1;it<nStep+1;it++)
           //     dis_1_assembled(it,ind)=Temp2(it);
        }        
    }
        
    // Assemble second order AAPG solutions
    printf("Assembling the second order terms...\n");
    dis_2_assembled = dis_1_assembled;
    int n = dim+1;
    if (AAPG_ord >= 2){
        for (unsigned int i=0;i<indi_2.XSize();i++){
                getRow(coeff2,n,coeffn);
                // Add second-order AAPG terms
                for (int p=1;p<PCTerms_2;p++){
                    // retrieve the index in the global matrix
                    int ind = ind2(indi_2(i),indj_2(i),p);
                    // retrieve the solution on dim (ii,jj)
                    Array1D<double> Temp1(nStep+1,0.e0);
                    getCol(dis_2(indi_2(i),indj_2(i)),p,Temp1);
                    // add PC terms to the corresponding place in dis_2_assembled
                    Array1D<double> Temp2(nStep+1,0.e0);
                    getCol(dis_2_assembled,ind,Temp2);
                    dotprodVec(coeffn,Temp1);
                    addVec(Temp1,Temp2);
                    dis_2_assembled.replaceCol(Temp2,ind);
                    //for (int it=1;it<nStep+1;it++)
                    //    dis_2_assembled(it,ind)=Temp2(it);
                }

                // Subtract first-order AAPG terms
                for (int p=1;p<PCTerms_1;p++){
                    subtractfirst(p,indi_2(i),ind1,dis_1,nStep,coeffn,dis_2_assembled);
                    subtractfirst(p,indj_2(i),ind1,dis_1,nStep,coeffn,dis_2_assembled);
                }
                n++;    
        }
    }    
    
    // Assemble third order AAPG solutions
    dis_3_assembled = dis_2_assembled;
    if (AAPG_ord >= 3){
        printf("Assembling the third order terms...\n");
        for (unsigned int ii=0;ii<indi_3.XSize();ii++){
                    // retrieve the index in the global matrix
                    Array1D<int> ind = ind3(indi_3(ii),indj_3(ii),indk_3(ii));
                    Array1D<double> coeffn(nStep+1,0.e0);
                    getRow(coeff3,n,coeffn);
                    
                    // Add third-order AAPG terms
                    for (int p=1;p<PCTerms_3;p++){
                        // retrieve the solution on dim (ii,jj,kk)
                        Array1D<double> Temp1(nStep+1,0.e0);
                        getCol(dis_3(indi_3(ii),indj_3(ii),indk_3(ii)),p,Temp1);
                        // add PC terms to the corresponding place in dis_3_assembled
                        Array1D<double> Temp2(nStep+1,0.e0);
                        getCol(dis_3_assembled,ind(p),Temp2);
                        dotprodVec(coeffn,Temp1);
                        addVec(Temp1,Temp2);
                        //dis_3_assembled.replaceCol(Temp2,ind(p));
                        for (int it=1;it<nStep+1;it++)
                            dis_3_assembled(it,ind(p))=Temp2(it);
                    }

                    // Subtract second-order AAPG terms
                    for (int p=1;p<PCTerms_2;p++){
                        subtractsecond(p,indi_3(ii),indj_3(ii),ind2,dis_2,nStep,coeffn,dis_3_assembled); 
                        subtractsecond(p,indj_3(ii),indk_3(ii),ind2,dis_2,nStep,coeffn,dis_3_assembled); 
                        subtractsecond(p,indi_3(ii),indk_3(ii),ind2,dis_2,nStep,coeffn,dis_3_assembled); 
                    }

                    // Subtract first-order AAPG terms
                    coeffn.Resize(nStep+1,-1); 
                    for (int p=1;p<PCTerms_1;p++){
                        subtractfirst(p,indi_3(ii),ind1,dis_1,nStep,coeffn,dis_3_assembled);
                        subtractfirst(p,indj_3(ii),ind1,dis_1,nStep,coeffn,dis_3_assembled);
                        subtractfirst(p,indk_3(ii),ind1,dis_1,nStep,coeffn,dis_3_assembled);
                    }
                    n++;
        }
    }
    return;
}

void computeStd(int nStep, int nPCTerms, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, PCSet& myPCSet, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3){
    // compute std and save
    // First-order 
    for (int i=0;i<nStep+1;i++){
        Array1D<double> dis_temp(nPCTerms,0.e0);
        getRow(dis_1_assembled,i,dis_temp);
        std1(i) = myPCSet.StDv(dis_temp);
    }
    write_datafile_1d(std1,"AAPG1_std.dat");
    // Second-order
    for (int i=0;i<nStep+1;i++){
        Array1D<double> dis_temp(nPCTerms,0.e0);
        getRow(dis_2_assembled,i,dis_temp);
        std2(i) = myPCSet.StDv(dis_temp);
    }
    write_datafile_1d(std2,"AAPG2_std.dat");
    // Third-order
    for (int i=0;i<nStep+1;i++){
        Array1D<double> dis_temp(nPCTerms,0.e0);
        getRow(dis_3_assembled,i,dis_temp);
        std3(i) = myPCSet.StDv(dis_temp);
    }
    write_datafile_1d(std3,"AAPG3_std.dat");
    return;
}

void subtractfirst(int p,int ii, Array2D<int>& ind1,Array1D<Array2D<double> >& dis_1,int nStep,Array1D<double>& coeffn,Array2D<double>& dis_assembled){
    int indii = ind1(p,ii);
    // retrieve PC solution on ii dim
    Array1D<double> Temp1(nStep+1,0.e0);
    getCol(dis_1(ii),p,Temp1);
    // incorporate coeffn
    dotprodVec(coeffn,Temp1);
    // subtract the first order terms and put back into dis_assembled
    Array1D<double> Temp2(nStep+1,0.e0);
    getCol(dis_assembled,indii,Temp2);
    subtractVec(Temp1,Temp2);
    dis_assembled.replaceCol(Temp2,indii);
    //for (int it=0;it<nStep+1;it++)
    //    dis_assembled(it,indii)=Temp2(it);

    return;
}

void subtractsecond(int p,int ii,int jj, Array3D<int>& ind2, Array2D<Array2D<double> >& dis_2, int nStep, Array1D<double>& coeffn, Array2D<double>& dis_assembled){
    int ind = ind2(ii,jj,p);
    // retrieve PC solution on ii/jj dim
    Array1D<double> Temp1(nStep+1,0.e0);
    getCol(dis_2(ii,jj),p,Temp1);
    // incorporate coeffn
    dotprodVec(coeffn,Temp1);
    // subtract the second order terms and put back into dis_assembled
    Array1D<double> Temp2(nStep+1,0.e0);
    getCol(dis_assembled,ind,Temp2);
    subtractVec(Temp1,Temp2);
    dis_assembled.replaceCol(Temp2,ind);
    //for (int it=0;it<nStep+1;it++)
    //    dis_assembled(it,ind)=Temp2(it);
    
    return;
}
