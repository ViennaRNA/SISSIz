/*

  SISSIz -- Dinucleotide controlled random alignments and RNA gene prediction

  Copyright (c) Tanja Gesell & Stefan Washietl, 2008

  See file COPYING for licence.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "expfit.h"
#include "utils.h"
#include "lm.h"

void expFunction (double *parVector, double *xVector, int m, int n, void *adat);
void expDer (double *parVector, double *jacobian, int m, int n, void *adat);

struct data {
  int n;
  double * pDat;
  double * dDat;
};

int fitExpFunction(double* parA, double* parB, const int N, double* pData, double* dData){

  double* parVector;
  double* xVector;
  double* jacobian;

  double parInitGuess[2] = {-1.33,0.75};

  int iter, i, status;

  struct data dat = {N, pData, dData};

  parVector=(double*)space(sizeof(double)*2);
  xVector=(double*)space(sizeof(double)*N);
  jacobian=(double*)space(sizeof(double)*2*N);

  parVector[0]=parInitGuess[0];
  parVector[1]=parInitGuess[1];

  status=dlevmar_der(&expFunction,&expDer,parVector,jacobian,2,N,50,NULL,NULL,NULL,NULL,&dat);
  
  *parA=parVector[0];
  *parB=parVector[1];


  free(parVector);
  free(xVector);
  free(jacobian);

  
  return 1;
}


void expFunction (double *parVector, double *xVector, int m, int n, void *adat){

  double *pDat = ((struct data *)adat)->pDat;
  double *dDat = ((struct data *)adat)->dDat;

  double A = parVector[0];
  double B = parVector[1];

  size_t i;

  for (i = 0; i < n; i++){
    double pEst = B * (1-exp(A * dDat[i]));
    xVector[i]=(pEst - pDat[i]);
  }
}

void expDer (double *parVector, double *jacobian, int m, int n, void *adat){

  double *pDat = ((struct data *)adat)->pDat;
  double *dDat = ((struct data *)adat)->dDat;

  double A = parVector[0];
  double B = parVector[1];

  size_t i,j;

  for (i=j=0; i < n; i+=1){

    double e = exp(A*dDat[i]);

    jacobian[j++]=(-1)*B*e*dDat[i];
    jacobian[j++]=(1-e); 

  }
}


