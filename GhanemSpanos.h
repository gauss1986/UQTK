Array2D<double> GS(PCSet& myPCSet, int order, int dim, int nPCTerms, string pcType, int nStep, double dis_0, double vel_0, double dTym, Array1D<double>& inpParams, Array2D<double>& f_GS);
void forward_duffing_GS(PCSet& myPCSet, Array1D<double>& inpParams, Array1D<double>& force, double dTym, Array1D<double>& dis, Array1D<double>& vel);
void RHS_GS(PCSet& myPCSet, Array1D<double>& force, Array1D<double>& dis, Array1D<double>& vel, Array1D<double>& dudt, Array1D<double>& dvdt, double epsilon, double zeta);
Array1D<double> postprocess_GS(int nPCTerms, int nStep, double dis0, Array2D<double>& dis_GS, PCSet& myPCSet, double dTym, FILE* GS_dump, FILE* GSstat_dump, Array2D<double>& mstd_MCS);
