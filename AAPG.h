Array1D<double> AAPG(int dof, Array1D<double> inpParams, Array1D<double>& fbar, double dTym, int order, string pcType, int noutput, int dim, int nStep, Array2D<double>& scaledKLmodes, double factor_OD, int AAPG_ord, bool act_D, double p, Array2D<double>& mstd_MCS, FILE* err_dump, Array2D<double>& sample_mstd_2D, Array2D<double>& samPts_norm, Array2D<double>& e_AAPG, Array1D<Array1D<double> >& e_sample_dis, Array1D<Array1D<double> >& e_sample_vel, Array1D<int>& init_D, Array1D<int>& coeff_D, bool PDF);
void PostProcess(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, int AAPG_ord, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, Array1D<double>& dis_1_mean, Array1D<double>& dis_2_mean, Array1D<double>& dis_3_mean, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3, int dim, int nStep, int PCTerms_1, int PCTerms_2, int PCTerms_3, int order, double dTym, double factor_OD, Array2D<double>& mstd_MCS, Array2D<double>& samPts_norm, string name, int noutput, Array1D<Array1D<double> >& e_sample, bool PDF, string pcType, Array2D<double>& sample_mstd_2D);
Array1D<int> identicalrow(int PCTerms, int dim, Array2D<int>& Pbtot, Array2D<int>& Pbmore);
Array2D<int> index1(int dim, int PCTerms_1, int order, Array2D<int>& Pbtot, Array2D<int>& Pb1);
Array2D<Array1D<int> > index2(int dim, int PCTerms_2, int order, Array2D<int>& Pbtot, Array2D<int>& Pb2);
Array3D<Array1D<int> > index3(int dim, int PCTerms_3, int order, Array2D<int>& Pbtot, Array2D<int>& Pb3);
void index(int dim, int order, int PCTerms_1, int PCTerms_2, int PCTerms_3, Array2D<int>& ind1, Array2D<Array1D<int> >& ind2, Array3D<Array1D<int> >& ind3);
Array2D<int> computeBasis(int dim,int order);
void assemblemean(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, int dim, int nStep, int AAPG_ord, Array1D<double>& dis_1_mean, Array1D<double>& dis_2_mean, Array1D<double>& dis_3_mean, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3);
void assemblerest(Array1D<int>& indi_2, Array1D<int>& indj_2, Array1D<int>& indi_3, Array1D<int>& indj_3, Array1D<int>& indk_3, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, int PCTerms_1, int PCTerms_2, int PCTerms_3, int dim, int order, int nStep, int AAPG_ord, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3);
void computeStd(int nStep, int nPCTerms, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, Array1D<double>& normsq, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3);
void subtractfirst(int p,int ii, Array2D<int>& ind1,Array1D<Array2D<double> >& dis_1,int nStep,Array1D<double>& coeffn,Array2D<double>& dis_assembled);
void subtractsecond(int p,int ii,int jj, Array2D<Array1D<int> >& ind2, Array2D<Array2D<double> >& dis_2, int nStep, Array1D<double>& coeffn, Array2D<double>& dis_assembled);
