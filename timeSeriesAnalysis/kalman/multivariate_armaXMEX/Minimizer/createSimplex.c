#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

void createSimplex(double *paramV,int numParams,double radius,double **p)
{

  
  /*double offSet = .3; */

  /* Initialize the random number generator with time and process ID */
  srand(time(NULL) + getpid());

  double  randMax = RAND_MAX;
  int j,k;
  double cRand;

  for ( j = 1; j <= numParams + 1; j ++){
    for (k = 1; k <= numParams; k++){

	/* Get a random number between -1 and 1 */
	cRand = ( (double) rand()  /   randMax) *2 - 1;
	p[j][k] = *(paramV + k) + cRand * radius;
	/* p[j][k] = *(paramV + k) + offSet; */

    }
  }
}
