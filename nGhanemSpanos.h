void nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck, int nStep, Array1D<Array1D<double> >& initial, double dTym, Array1D<Array2D<double> >& f_GS, Array1D<Array2D<double> >& solution);
void forward_duffing_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck, Array1D<Array1D<double> >& force_current, Array1D<Array1D<double> >& force_mid, Array1D<Array1D<double> >& force_plus,  double dTym, Array1D<Array1D<double> >& x);
Array1D<Array1D<double> > kbar(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& epsilon, Array1D<Array1D<Array1D<double> > >& mck);
void dev_nGS(int dof, PCSet& myPCSet, Array1D<Array1D<double> >& force, Array1D<Array1D<Array1D<double> > >& mck, Array1D<Array1D<double> >& u, Array1D<Array1D<double> >& v, Array1D<Array1D<double> >& du, Array1D<Array1D<double> >& dv, Array1D<Array1D<double> >& kbar);
Array2D<double> postprocess_nGS(int dof, int nStep,  Array1D<Array2D<double> >& solution, PCSet& myPCSet, double dTym, int ord,Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array1D<Array2D<double> >& mstd_dis, Array1D<Array2D<double> >& mstd_vel, Array2D<double>& e2, Array1D<double>& e3);
Array2D<double> nerror(string info, int dof,int nStep,Array1D<Array2D<double> >& et,Array2D<double>& mean,Array2D<double>& StDv, Array2D<double>& mean_MCS, Array2D<double>& std_MCS, Array2D<double>& e2, Array1D<double>& e3);
