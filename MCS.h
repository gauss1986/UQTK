Array2D<double> det(int dof, int nspl,  int nStep, int nkl, double dTym, Array1D<double>& totalforce, Array1D<double>& inpParams, Array1D<double>& initial);
void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force, double dTym, Array1D<double>& x);
Array1D<double> RHS(double force, Array1D<double>& x, Array1D<double>& inpParams);
Array1D<double> mStd(Array1D<double>& x,int nspl);
Array1D<double> error(Array1D<double>& dis, Array1D<double>& StDv, Array2D<double>& mstd_MCS);
Array1D<double>  sample_force(Array2D<double>& samPts,int iq,int nStep,double fbar,int nkl,Array2D<double>& scaledKLmodes, Array1D<double>& inpParams); 
