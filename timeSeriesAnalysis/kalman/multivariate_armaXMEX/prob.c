/* Structure for ARMAX/CARMA model parameters */


struct probCDT{

  double *TRAJ;
  int trajLength;
  int arOrderMax;
  int maOrderMax;
  int *maBIN;
  int *topoBIN;
  int nNodes;
  int nInputs;
  double wnVariance;
  int numParams;
  int numMissing;

};
