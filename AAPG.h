void AAPG(Array1D<double> inpParams, double fbar, double dTym, int order, string pcType, int dim, int nStep, Array2D<double>& scaledKLmodes, double dis0, double vel0, PCSet& myPCSet, double factor_OD, int AAPG_ord);
void PostProcess(int AAPG_ord, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, PCSet& myPCSet, double fbar, int dim, int nStep, int PCTerms_1, int PCTerms_2, int PCTerms_3, int order, double dTym, string pcType, Array1D<double>& inpParams, Array2D<double>& scaledKLmodes, double factor_OD);
Array1D<int> index1(int dim, int i, int order);
Array1D<int> index2(int dim, int i, int j, int order);
Array1D<int> index3(int dim, int i, int j, int k, int order);
Array2D<int> computeBasis(int dim,int order);
void sampledet(double fbar, int PA, int dim, int nStep, double dTym, Array1D<double>& inpParams, Array2D<double>& samPts, Array2D<double>& scaledKLmodes, Array2D<double>& dis_det);
void sampleAAPG(int AAPG_ord, Array2D<double>& samPts, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, int PCTerms_1, int PCTerms_2, int PA, int dim, int order, string pcType, int nStep, Array3D<double>& samAAPG);
void coeff(int AAPG_ord, Array2D<double>& coeffAAPG1, Array2D<double>& coeffAAPG2, double fbar, int PA, int dim, int order, int nStep, double dTym, string pcType, Array1D<double>& inpParams, Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, int PCTerms_1, int PCTerms_2, int nAAPGTerms, Array2D<double>& scaledKLmodes);
void assemblemean(Array1D<double>& dis_0, Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, int dim, int nStep, int AAPG_ord, Array1D<double>& dis_1_mean, Array1D<double>& dis_2_mean, Array1D<double>& dis_3_mean, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3);
void assemblerest(Array1D<Array2D<double> >& dis_1, Array2D<Array2D<double> >& dis_2, Array3D<Array2D<double> >& dis_3, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, int PCTerms_1, int PCTerms_2, int PCTerms_3, int dim, int order, int nStep, int AAPG_ord, Array2D<double>& coeff1, Array2D<double>& coeff2, Array2D<double>& coeff3);
void computeStd(int nStep, int nPCTerms, Array2D<double>& dis_1_assembled, Array2D<double>& dis_2_assembled, Array2D<double>& dis_3_assembled, PCSet& myPCSet, Array1D<double>& std1, Array1D<double>& std2, Array1D<double>& std3);
void subtractfirst(int p,int ii, Array2D<int>& ind1,Array1D<Array2D<double> >& dis_1,int nStep,Array1D<double>& coeffn,Array2D<double>& dis_assembled);
void subtractsecond(int p,int ii,int jj, Array3D<int>& ind2, Array2D<Array2D<double> >& dis_2, int nStep, Array1D<double>& coeffn, Array2D<double>& dis_assembled);
