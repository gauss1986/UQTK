double MCS(int nspl, int dim, int nStep, int nkl, double dTym, double fbar, Array2D<double>& scaledKLmodes, Array1D<double>& inpParams, Array2D<double>& samPts, Array2D<double>& dis_MC);
void forward_duffing_dt(Array1D<double>& inpParams, Array1D<double>& force, double dTym, Array1D<double>& x);
Array1D<double> RHS(Array1D<double>& force, Array1D<double>& x, Array1D<double>& inpParams);
Array1D<double> mStd(Array1D<double>& x,int nspl);
Array1D<double> error(Array1D<double>& dis, Array1D<double>& StDv, Array2D<double>& mstd_MCS);
void sample_duffing(Array2D<double>& samPts,int iq,int nStep,double fbar,int nkl,Array2D<double>& scaledKLmodes,Array2D<double>& totalforce);
void sample_lorenz(Array2D<double>& samPts,int iq,int nStep,double fbar,double Amp,double w,int nkl,Array2D<double>& scaledKLmodes,Array2D<double>& totalforce); 
