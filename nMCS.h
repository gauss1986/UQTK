Array1D<Array1D<double> > ndet(int dof, int nStep, double dTym, Array2D<double>& sampleforce, Array1D<double>& epsilon, Array1D<Array1D<double> >& mck, Array1D<double>& initial);
void nforward_duffing_dt(Array1D<double>& epsilon, int dof, Array2D<double>& force, Array1D<Array1D<double> >& mck, double dTym, Array1D<double>& u, Array1D<double>& v);
void nRHS(Array1D<double>& acc, int dof, Array1D<double>& epsilon, Array1D<Array1D<double> >& mck, Array1D<double>& force, Array1D<double>& u, Array1D<double>& v);
Array2D<double>  nsample_force(int dof, Array2D<double>& samPts,int iq, Array1D<double>& fbar,int nkl,Array2D<double>& scaledKLmodes, Array1D<double>& m);    
void mstd_MCS(Array1D<double>& result, double& mean, double& std);
