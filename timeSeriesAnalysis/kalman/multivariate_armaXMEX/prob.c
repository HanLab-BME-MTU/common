/* Structure for ARMAX/CARMA model parameters */



struct probCDT{

  double *TRAJ;
  int trajLength;
  int nNodes;
  int numMissing;

  int numParams;
  int arOrderMax;
  int maOrderMax;
  int *topoBIN;
  int *maBIN;

  int nInputs;

  double wnVariance;
};

/*

struct movie {
  double *TRAJ;
  int trajLength;
  int nNodes;
  int numMissing;
};

*/

/*
 * struct probCDT
 * --------------
 * The first field will be modified as:
 *
 * struct movie *data
 *
 * Data contains an array of movies.
 */

/*

struct probCDT{

  struct movie data;
  int nMovies;

  int numParams;
  int arOrderMax;
  int maOrderMax;
  int *topoBIN;
  int *maBIN;

  int nInputs;

  double wnVariance;
};

*/


