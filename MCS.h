double MCS(int nspl, int dim, int nStep, int nkl, double dTym, double fbar, Array2D<double>& scaledKLmodes, Array1D<double>& inpParams, Array2D<double>& samPts, Array2D<double>& dis_MC);
void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force, double dTym, Array1D<double>& x);
void RHS(Array1D<double>& force, Array1D<double>& x, Array1D<double>& dxdt, double epsilon, double zeta);
Array1D<double> mStd(Array1D<double>& x,int nspl);
