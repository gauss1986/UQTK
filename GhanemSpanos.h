void GS(int dof, PCSet& myPCSet, int order, int dim, int nPCTerms, string pcType, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS, Array1D<Array2D<double> >& solution);
void forward_duffing_GS(int dof, PCSet& myPCSet, Array1D<double>& inpParams, Array1D<double>& force_current, Array1D<double>& force_mid, Array1D<double>& force_plus,  double dTym, Array1D<Array1D<double> >& x);
void RHS_GS(int dof, PCSet& myPCSet, Array1D<double>& force, Array1D<Array1D<double> >& x, Array1D<double>& inpParams, Array1D<Array1D<double> >& dev);
Array1D<double> postprocess_GS(int noutput, int nPCTerms, int nStep,  Array2D<double>& solution, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS);
Array2D<double> sampleGS(int noutput, int dim, int nStep, int nPCTerms, PCSet& myPCSet, Array2D<double>& solution, Array2D<double>& samPts);
