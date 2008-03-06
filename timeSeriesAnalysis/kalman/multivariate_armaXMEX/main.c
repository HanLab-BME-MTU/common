/*
 * File: main.c
 * ------------
 * This main module is to demonstrate how covKalmanInit should be used
 * in the context of large program, as well as to compare the results
 * with those from MATLAB implementation.
 * 
 * Note:
 *
 * (1) The dimension of the initial covariance matrix is m x m, where
 *     m = max(p, q + 1), p = autoregressive order, q = moving average order.
 *     (Jones, equation 2.3)
 *
 * (2) Memory allocation should be done before calling covKalmanInit.
 *
 * Reference: Jones, R.H. (1980) Technometrics 22(3): 389-395. 
 *
 * Pei-hsin Hsu, March 2008
 */


#include <stdio.h>
#include <stdlib.h>
#include "covKalmanInit.h"


/* Private function prototypes */

static double **NewMat(int row, int col);
static void FreeMat(double **mat, int row);
static void PrintMat(double **mat, int row, int col);


/*
 * Function: main
 * --------------
 * (1) Local variable cov is the initial state covariance matrix
 *     whose dimension is m x m, where m = max(p, q + 1). Wrong
 *     dimension will cause segmentation fault.
 * (2) The result is passed by reference.
 */

main(void)
{
  double ar[2] = {1.5, 1.2};
  double ma[2] = {0.8, 0.4};
  double G[3] = {0.9, 0.6, 0.3};
  double **cov;

  cov = NewMat(3, 3); /* memory allocation */

  covKalmanInit(ar, ma, G, 2, 2, 3, cov);

  PrintMat(cov, 3, 3);

  FreeMat(cov, 3);
}


/* 
 * Function: NewMat
 * ----------------
 * Allocate memory for a new matrix
 */

static double **NewMat(int row, int col)
{
  int i;
  double **mat;

  mat = (double **) malloc(sizeof(double *) * row);
  for (i = 0; i < row; i++) mat[i]  = (double *) malloc(sizeof(double) * col);
  return mat;
}


/* 
 * Function: FreeMat
 * -----------------
 * Free the memory associated with each row ptr,
 * as well as the matrix ptr
 */

static void FreeMat(double **mat, int row)
{
  int i;

  for (i = 0; i < row; i++) free(mat[i]);
  free(mat);
}


static void PrintMat(double **mat, int row, int col)
{
  int i, j;

  printf("\n");
  for (i = 0; i < row; i++){
    for (j = 0; j < col; j++) printf("%-12.6g", mat[i][j]);
    printf("\n");
  }
  printf("\n");
}
