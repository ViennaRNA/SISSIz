#ifndef _MUTATE_H_
#define _MUTATE_H_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include "dkiss.h"
#include "global.h"
#include "treefile.h"
#include "models.h"
#include "evolve.h"



  #ifdef SPRNG
  #   if SPRNG != 0
  #      undef  SPRNG
  #      define SPRNG 1
  #   endif
  #else
  #   define SPRNG 0
  #endif


#if SPRNG
#  include "sprng/sprng.h"
#endif  /* SPRNG */

#define MUTATION 0
#define sissiDI01 1
 #define NEARZERO 1.0e-10

extern char  filename[256];
extern int   printfile, fastqx;

extern int catnei, numseq;
double timeoverall(double rsum) ;   /* TIME FOR THE NEXT MUATATION                        */
/*char   *copyseq(char *rseq);*/
char   *changeseq(char *seqi, int secpos, int mutation);
void   AllocateMemory();
void   RandomSequence_independent(char *seq);
void   RandomSequence(char *seq);
void   RandomSequenceMarkov1(char *seq);
void   CreateSequences(TTree *tree);
void   SetSequence(char *seq, char *source);
extern void RandomSequence_Markov(char *seq);
extern void RandomSequence_independent(char *seq);
extern void MutateSequenceInd(char *seq, double len, double *Qij);
extern void MutateSequenceIndqx(char *seq, double len, double *Qij);

extern double sqx;
extern int scalesequencerate;

/* prototypes */
/*void        MutateSequence(char *seq, double len, double *Qij, double **QSHij, double **QTij,double **QT2ij,double **QT3ij);*/
void MutateSequenceZ(char *seq, double len, double **QT2ij);
void MutateSequenceZqx(char *seq, double len, double **QT2ij);

 /* Compute rates */
double *rate(char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij);
double *rateDI(char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij); /*sissiDI01 1 */
double *rateDIz(char *seq, double **QT2ij, double *rat);
double  rate_sum(double rat[]);
double *rate_independent(char *seq, double *Qij, double *ra);
extern int fastrates;
extern double *scalingfactor;
extern double *scalingfactor_sum;
extern int sitefactor;
extern int verboseScalingFactor;
extern int verboseunnormScalingFactor;
/* PROBABILITY TO SELECT THE POSITION  */  
/* double *selsite(double rat[], double rsum, double *proselsite); */ 
inline   int selsite(double rat[], double rsum, double *proselsite);           
int    rand_choose(double proselsite[]);                 /* RANDOM CHOOSE: Select one positon for the next mutation  */
/* calculation of substitutions:  */
double *probmut( int secpos, double *Qij, char *seq);    /* INDEPENDENT case*/
double *probmutqx( int secpos, double *Qij, char *seq);    /* INDEPENDENT case*/
int     choosenuc(double *probmutij);
double *probmutSH(int secpos, double **D, char *seq);    /* nk=1 */
int     choosedoublet(double *probmutij_doub);
double *probmutTriplet(int secpos, double **D, char *seq); /* nk=2 */
int choosetriplet(double *probmutij_trip);
int tripletnucleotidmutation(int tripletmutation, int secpos, char *seq);
/*     MIX       */
double *probmutSHmix(int secpos, double **D, char *seq);
int     choosenucmix(double *probmutij, double *probmutij_dmix, int secpos, char *seq );
int     tripletnucleotidmutation(int tripletmutation, int secpos, char *seq);
/* Triplett  with categorie  */
/*double *cprobmutTriplet(int pos, double **QTij,double **QT2ij,double **QT3ij, char *seq); */
double *cprobmutTripletz(int pos,double **QT2ij, char *seq);     
int cchoosetriplet(double *probmutij_trip);
int ctripletnucleotidmutation(int tripletmutation, int secpos, char *seq); 
 
double** computePij(double **QT2ij);

#endif /* _MUTATE_H_ */
