/* Structure for ARMAX/CARMA model parameters */


struct probCDT{

  double *TRAJ;
  int trajLength;
  int arOrderMax;
  int maOrderMax;
  int *maBIN;
  int *topoBIN;
  int nNodes;
  double wnVariance;
  int numParams;

};
