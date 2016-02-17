void nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<Array2D<double> >& f_GS, Array1D<Array2D<double> >& solution);
void forward_duffing_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck, Array1D<Array1D<double> >& force_current, Array1D<Array1D<double> >& force_mid, Array1D<Array1D<double> >& force_plus,  double dTym, Array1D<Array1D<double> >& x);
Array1D<Array1D<double> > kbar(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<double> >& mck);
void dev_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& force, Array1D<Array1D<double> >& mck, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& v, Array1D<Array1D<double> >& du, Array1D<Array1D<double> >& dv, Array1D<Array1D<double> >& kbar);
Array2D<double> postprocess_nGS(int dof, int noutput, int nPCTerms, int nStep,  Array1D<Array2D<double> >& solution, PCSet& myPCSet, double dTym, int ord, Array1D<Array2D<double> >& mstd_MCS_u, Array1D<Array2D<double> >& mstd_MCS_v);
