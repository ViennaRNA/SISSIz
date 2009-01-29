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
#include "mutate.h"

#define sissiDI01 1
#define NEARZERO 1.0e-10

FILE  *dinucright_fv;
FILE  *dinucleft_fv;
/* FILE  *scalingf_fv;*/
FILE  *dh_fv;
FILE  *qx_fv;
char  filename[256];
int   printfile, fastqx;

int catnei;
int numseq;                 /* number of recent sequence, resp.number of taxons */

double *scalingfactor;
double *scalingfactor_sum;
int    sitefactor;
int    verboseScalingFactor;
int    verboseunnormScalingFactor;
double sqx;
int    scalesequencerate;

 /* Compute rates */
 double *rate(char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij);
 double *rateDI(char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij);
   void RandomSequence_Markov(char *seq);
 double *rateDIz(char *seq, double **QT2ij,double *rat);

 double  rate_sum(double rat[]);
 double  rate_sum_sf(double rat[]);
double  *rate_independent(char *seq, double *Qij,double *rat);
 int     fastrates;
 /* PROBABILITY TO SELECT THE POSITION  */
/* double *selsite(double rat[], double rsum, double *proselsite);     */   
 inline int selsite(double rat[], double rsum, double *proselsite);        
 int    rand_choose(double proselsite[]);                 /* RANDOM CHOOSE: Select one positon for the next mutation  */
 /* calculation of substitutions:  */
 double *probmut( int secpos, double *Qij, char *seq);    /* INDEPENDENT case*/
  double *probmutqx( int secpos, double *Qij, char *seq);    /* INDEPENDENT case*/
 int     choosenuc(double *probmutij);
 double *probmutSH(int secpos, double **D, char *seq);    /* nk=1 */
 int doubletnucleotidmutation(int doubletmutation, int secpos, char *seq);
 int     choosedoublet(double *probmutij_doub);
 double *probmutTriplet(int secpos, double **D, char *seq); /* nk=2 */
 int choosetriplet(double *probmutij_trip);
 int tripletnucleotidmutation(int tripletmutation, int secpos, char *seq);
 /*     MIX       */
 double *probmutSHmix(int secpos, double **D, char *seq);
 int     choosenucmix(double *probmutij, double *probmutij_dmix, int secpos, char *seq );
 int     tripletnucleotidmutation(int tripletmutation, int secpos, char *seq);
 /* Triplett  with categorie  */
/* double *cprobmutTriplet(int pos, double **QTij,double **QT2ij,double **QT3ij, char *seq); */
 double *cprobmutTripletz(int pos,double **QT2ij, char *seq);      
 double *cprobmutTripletzqx(int pos,double **QT2ij, char *seq); 
 int cchoosetriplet(double *probmutij_trip);
 int ctripletnucleotidmutation(int tripletmutation, int secpos, char *seq); 
 /*************************************************/
 double timeoverall(double rsum) ;   /* TIME FOR THE NEXT MUATATION                        */
/* char *copyseq(char *rseq);*/
 char *changeseq(char *seqi, int secpos, int mutation);
 void AllocateMemory();
 void RandomSequence_independent(char *seq);
 void RandomSequence(char *seq);
  void RandomSequenceMarkov1(char *seq);
 void CreateSequences(TTree *tree);
 void SetSequence(char *seq, char *source);
 double getrandomdouble();
void RandomSequence_independent(char *seq);
  /*************************************************/
double *newcurrentrate(int changepos, double *currentrat, char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij); 
double *changerate(int currentpos, double *currentrat, char *seq, double *Qij, double **QSHij, double **QTij, double **QT2ij, double **QT3ij);
double *changerateDIz(int secpos, double *rat, char *seq, double **QT2ij);
double minusrsumDIz(int secpos, double *rat, double rsum);
double plusrsumDIz(int secpos, double *rat, double rsum);
int    selsite_uniform();
double *changerateInd(int secpos, double *rat, char *seq, double *Qij);
double minusrsumInd(int secpos, double *rat, double rsum);
double plusrsumInd(int secpos, double *rat, double rsum);
void   MutateSequenceInd(char *seq, double len, double *Qij);
void   MutateSequenceIndqx(char *seq, double len, double *Qij);

double** computePij(double **QT2ij);


 double getrandomdouble()
 {
#  if SPRNG
     if(randomgeneratordkiss == 1)
#  endif    
	   return(dkiss());               /*random with dkiss()*/
#  if SPRNG	   
     else	
       return(sprng(randstream));  
#  endif     
 /*#       else
	   return(sprng(randstream));*/
 /*#       endif*/
 }

 /****************************************************************************************************************/
 /*                                 TIME                                                                         */
 /****************************************************************************************************************/
 /******************************************************/
 /* TIME FOR THE NEXT MUATATION                        */
 /* IS EXPONENTIAL DISTRIBUTE OVER THE SUM OF THE RATES*/
 /******************************************************/
 double timeoverall(double rsum)
 {
         double t_expsum;
         double a;
		 
           a=1-getrandomdouble();
		   /*1-dkiss();
		   if(a==0 || a==1) a=0.5;
          /* t_expsum=(log(a)*2)/rsum;*/
		   t_expsum=(log(a))/rsum;
		  /* if(verbosewaiting==1) fprintf(stderr,"t_expsum: %f\n", t_expsum);*/
	
		return(t_expsum);
  }
 /****************************/
 /* COPY SEQI TO SEQ J       */
 /****************************/
/* char *copyseq(char *rseq){
            int pos;
            char *seqi;

            seqi=(char*)malloc(seqlen*sizeof(char));
            for(pos=0;pos<seqlen;pos++)
              seqi[pos]=rseq[pos];
            
			return(seqi);
  }*/
 /**************************************************************/
 /*CHANGE NUCLEOTID POSITION I TO J  INDEPENDENT=DEPENDENT CASE*/
 /**************************************************************/
 char *changeseq(char *seqi, int secpos, int mutation)
 {
   seqi[secpos]=nucleotide[mutation];
   return(seqi);
 }
 /**************************************************************/
 /*                    SEQUENCE                                */
 /**************************************************************/		
 /*
  * Random Sequence indpendent
  */		
/*  void RandomSequence_independent(char *seq)
  {
	char *P, *rseq;
	double p[4], x;
    int i, nr, ink, pos;

	rseq=(char*)calloc(seqlen,sizeof(char));
	P=seq;
	
	   p[0]=0;
       p[1]=mpi[0];
       p[2]=mpi[0]+mpi[1];
       p[3]=mpi[0]+mpi[1]+mpi[2];

	for(pos=0; pos<seqlen; pos++)
	{
			do{
			   nr=0;
			   x=getrandomdouble();
			   for(ink=1; ink<4;ink++)
				  if(x>p[ink]) nr=ink;
			}while(nr<0 || nr>3);

			rseq[pos]=nucleotide[nr];
	}
           	
	for (i=0; i<seqlen; i++) {
		*P=rseq[i];
		P++;
	}
	if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", rseq[pos]);
	
	free(rseq);
 }*/
  /**************************************************************/		
  /*
  * Random Sequence Markov 1 for Dinucleotides
  */
  void RandomSequenceMarkov1(char *seq){
	double p[4],markovd[4][4],mp[4][4], x, summarkov;
    int i,j,k, nr,nr2, ink, pos;
	/*double **D;
	double sumdinuc;*/
	
	p[0]=0;
	p[1]=mpi[0];
	p[2]=mpi[0]+mpi[1];
	p[3]=mpi[0]+mpi[1]+mpi[2];

	do{
	    nr=0;
		x=getrandomdouble();
		for(ink=1; ink<4;ink++)
				  if(x>p[ink]) nr=ink;
				 /* fprintf(stderr, "\n %d ", nr);*/
	}while(nr<0 || nr>3); 

	seq[0]=nucleotide[nr];
	
	
	
	
/*	fprintf(stderr, "\n");
	for(i=0;i<16;i++) fprintf(stderr, "%lf ", dipi[i]);
	fprintf(stderr, "\n");*/
	/*for(i=0;i<4;i++) fprintf(stderr, "%lf ", mpi[i]);
	fprintf(stderr, "\n");*/

	   
	   
	
	   
	   
	 k=0;
	 for(i=0;i<4;i++){
	    summarkov=0.00;
	    for(j=0;j<4;j++){
		   markovd[i][j]=(dipi[j+(4*k)]/mpi[i]);
		 /*  fprintf(stderr, "markovd %lf = dipi %lf/ mpi %lf sum %lf \n",    markovd[i][j], dipi[j+(4*k)], mpi[i], summarkov); */
		   summarkov=summarkov+ markovd[i][j];
	    }
		if((summarkov < 1.00000 - NEARZERO) || (summarkov > 1.00000 + NEARZERO)){	
		  if(!quiet) fprintf(stderr, "Warning: %lf rowsum of RANDOM MARKOV", summarkov);
		  for(j=0;j<4;j++){
		      markovd[i][j]=markovd[i][j]/summarkov;
	       }
        }
		k++;
      }  
	  summarkov=0.00;
	  if(verboseMarkov==1){
		   for(j=0;j<4;j++){
		      fprintf(stderr, "%lf ", markovd[i][j]);
		       summarkov=summarkov+ markovd[i][j];
		   }
		   fprintf(stderr, "=%lf\n",	summarkov);   
	 }	   
	    		   

	   
	 /*  fprintf(stderr, "\nMarkovd hallo d\n");
	   for(i=0;i<4;i++){ 
	     for(j=0;j<4;j++){
	     fprintf(stderr, "\n%lf ", markovd[i][j]);
		 }
	   }*/	 
	
	  for(i=0;i<4;i++){
	    mp[i][0]=0;
	    mp[i][1]=markovd[i][0];
	    mp[i][2]=markovd[i][0]+markovd[i][1];
	    mp[i][3]=markovd[i][0]+markovd[i][1]+markovd[i][2];
	 }
	 /*   summarkov=0.00;
	   for(i=0;i<4;i++){
		 summarkov=0.00;
	     for(j=0;j<4;j++){
		    fprintf(stderr, "\n%lf ", mp[i][j]);
		 }
		 summarkov= 1.00-mp[i][3];
		 fprintf(stderr, "rest %lf= %lf\n", (summarkov), markovd[i][3]);	
	  }		*/
	for(pos=1; pos<seqlen; pos++)
	{
			do{
			   nr2=0;
			   x=getrandomdouble();
			   for(ink=1; ink<4;ink++)
				  if(x>mp[nr][ink]) nr2=ink;
			}while(nr2<0 || nr2>3);

			/*rseq[pos]=nucleotide[nr2];*/
			seq[pos]=nucleotide[nr2];
			nr=nr2;
	}
  /*  for (i=0; i<seqlen; i++) {
		*P=rseq[i];
		P++;
	}*/
	
  /* if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", rseq[pos]);
   if(verboserandom==1) fprintf(stderr, "\n");
   if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", rseq[pos]);*/
   
    if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", seq[pos]);
   if(verboserandom==1) fprintf(stderr, "\n");
   if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", seq[pos]);
   
   
  /*   D=(double**)malloc(4*sizeof(double*));
	  for(i=0;i<4;i++)
        D[i]=(double *)malloc(4*sizeof(double));
	  for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[i][j]=0.00;
	  }*/
	 /*LEFT*/	
/*	 for(i=1;i<numBases;i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i-1])]++;
	 for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=(numBases-1);
	 }
	 sumdinuc=0.00;
	 for(i=0;i<4;i++)
	      for(j=0;j<4;j++) sumdinuc += D[j][i];
	 if(sumdinuc !=1.00){
	   	  for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
	 }	  
	 
	 for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
		    fprintf(stderr,"%f ", D[j][i]);
		  }	
	  }	
	  fprintf(stderr, "\n");*/


   
 /* free(rseq);*/
   
 }
/******************************************************************************/
/******************************************************************************/
void CreateSequences(TTree *tree)
{
	TNode *P;
	
	P=tree->nodeList;
	while (P!=NULL) {
	    if(P->sequence != NULL) free(P->sequence);
		P->sequence=CAllocMem(numBases+1, "sequences", "CreateSequences", 0);
		P=P->next;
	}
}
/******************************************************************************/
void SetSequence(char *seq, char *source)
{
	int i;
	char *P, *Q;
	
	P=seq;
	Q=source;
	
	for(i=0;i<seqlen;i++){
	*P=source[i];
	P++;
	}
}
	
 /**************************************************************/
/*             MUTATE  SEQUENCE                               */
/**************************************************************/	
void MutateSequenceZ(char *seq, double len, double **QT2ij){
  /*char    *P;*/
  char    *ancestorseq;
  double  *rat;                                  
  double  rsum;
 /* double  rsum_sf;*/
  double  t=0;
  int secpos,i,j;
  /*int *secneighpositions;    */    
  int     substitutions, osubstitutions;
 /* double  *pro_selsite;*/
  double  *probmuttriple;
  int     mutation;
  int     tripletmutation; 
  double  **D;
  char dinucFilename[256];
  char dinucRightFilename[256];
 /* char scalingfFilename[256];*/
  char dhFilename[256];
  double sumdinuc;
  double *proselsite;  /*stefan camb*/

     

  proselsite=(double *)malloc(seqlen*sizeof(double));  /*stefan camb*/
  ancestorseq =(char*)malloc(numBases*sizeof(char));
  for(i=0; i<numBases; i++)
		ancestorseq[i]=seq[i];
		
  rat=(double *)calloc(seqlen,sizeof(double));
  
	/*P=seq;*/
	/*if(simindependent==1) rat = rate_independent(seq,Qij)
    else */
	rat  = rateDIz(seq,QT2ij,rat); 
    rsum =0.0;	
	rsum = rate_sum(rat);
	   
  /*if(verbosewaiting==1) fprintf(stderr,"\n rsum1= %f\n",rsum);*/
	/*************Time*****************************************/
	t    = timeoverall(rsum);
	substitutions=0;
	
	/*for(i=0;i<seqlen;i++){
    if(sitefactor==0) proselsite[i] = 1/(double)seqlen ;
    else  proselsite[i]= scalingfactor[i] ;
		}*/
  
  /**********************************              
	        EVOLUTION PROCESS                        
          ********************************* */
  if(t>len && verboseWarning==1){
		fprintf(stderr, "No substitution! (waiting time:%f > branch length:%f)\n", t, len);
  }
	
  while(t<=len)
    {
      substitutions++; 
      /* pro_selsite=NULL; */ 
      /*pro_selsite=selsite(rat,rsum);    */    /*probability for position selection*/
      secpos=-1;
      /*selsite(rat,rsum,proselsite);  
        secpos  = rand_choose(proselsite); */  /*actuell position*/
      secpos=selsite(rat,rsum,proselsite);
      /*printf("%d\n",secpos);*/
      /*********************************/
	    /*which subsitution*/
      mutation=-1;
      /*probmuttriple        = cprobmutTriplet(secpos, QTij,QT2ij, QT3ij, seq);*/
      probmuttriple        = cprobmutTripletz(secpos, QT2ij, seq);
      
      tripletmutation      = cchoosetriplet(probmuttriple);
      /*fprintf(stderr, "triplet %d : \n", tripletmutation);*/
      mutation             = ctripletnucleotidmutation(tripletmutation, secpos, seq);
      /*		fprintf(stderr, "->single mutation %d\n",  mutation);*/
      
      /*********************************/
      seq=changeseq(seq,secpos,mutation);
      /*P=seq;*/
      /*******************************************************/
      /* free(secneighpositions);*/
      /*  free(rat);*/
      /* rat  = rateDIz(seq,QT2ij); */
      /*double *changerateDIz(int secpos, double *rat, char *seq, double **QT2ij)*/
      if(fastrates==1){
        /* fprintf(stderr, "rsum %lf ",rsum);*/
        rsum=minusrsumDIz(secpos,rat, rsum);
        /*	fprintf(stderr, "  %lf",rsum);*/
        rat=changerateDIz(secpos, rat, seq,QT2ij);
        rsum=plusrsumDIz(secpos,rat, rsum);
        /*	fprintf(stderr, "  %lf\n",rsum);*/
      }	
      else{ 
        rat  = rateDIz(seq,QT2ij,rat);
        rsum= -1;
        rsum = rate_sum(rat);
      }	 
      /* if(verbosewaiting==1) fprintf(stderr, "\n rsum= %f\n",rsum);*/
      /* fprintf(stderr, "rsum= %f\n",rsum);*/
      t += timeoverall(rsum);
      /********************************************************/
    }	
	free(rat);
	/*free(probmuttriple);*/
	osubstitutions=0;
	for(i=0;i<numBases;i++){
		if(ancestorseq[i]!=seq[i]) osubstitutions=osubstitutions+1;
		/*if(P[i]!=seq[i]) osubstitutions=osubstitutions+1;*/
	}
	if(verboseEvolutiontime==1){
	  if(len==0.00) fprintf(stderr, "branch length: %f - no substituion\n", len);
	  else{
      if(printfile==1) {	
        sprintf(dhFilename, "%s_dh.sissi", filename);
        if((dh_fv=fopen(dhFilename,"a"))==NULL)
          fprintf(dh_fv, "dinuc.txt could not open\n");	
        fprintf(dh_fv ,  "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
      }else fprintf(stderr , "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
    }
	}
	/*print DiNucContent*/
	if((verboseDiNucContent==1 || (verboseDiNucContent==1 &&printfile==1)) && len>0){
	  /*allocate 4*4 matrix*/
	  D=(double**)malloc(4*sizeof(double*));
	  for(i=0;i<4;i++)
      D[i]=(double *)malloc(4*sizeof(double));
	  for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[i][j]=0.00;
	  }
    /*LEFT*/	
    for(i=1;i<numBases;i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i-1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/	
    if(printfile ==1 ){
      /* for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(stderr, "%c%c ", nucleotide[i],nucleotide[j]);*/
      sprintf(dinucFilename, "%s_dinuc_left.sissi", filename);
      if((dinucleft_fv=fopen(dinucFilename,"a"))==NULL)
	      fprintf(stderr, "dinuc.txt could not open\n");
      /* fprintf(dinucleft_fv,"%f \n", sumdinuc);	  
         for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(dinucleft_fv, "%c%c ", nucleotide[i],nucleotide[j]);
         fprintf(dinucleft_fv," \n");	*/
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucleft_fv,"%f ", D[j][i]);
        }	
      }	
      fprintf(dinucleft_fv, "\n");
	  } 
	  /*print*/
	  else{    
      for(i=0;i<4;i++){
		    for(j=0;j<4;j++) {
          fprintf(stderr ,"%f ", D[j][i]);
          fprintf(stderr, "%c%c ", nucleotide[j],nucleotide[i]);
        }   
      }
      fprintf(stderr, "\n"); 
	  }   	 
    /*Right*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]=0.00;
    }
    for(i=0;i<(numBases-1);i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i+1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/
    if(printfile==1){	
      sprintf(dinucRightFilename, "%s_dinuc_right.sissi", filename);
      if((dinucright_fv=fopen(dinucRightFilename,"a"))==NULL)
        fprintf(stderr, "dinuc.txt could not open\n");	  
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          fprintf(dinucright_fv,"%f ", D[j][i]);
        }	
      }	 
      fprintf(dinucright_fv, "\n");	
    }/*print */
    else{
      for(i=0;i<4;i++){
        for(j=0;j<4;j++) 
          fprintf(stderr,"%f ", D[j][i]);
      }	 
      fprintf(stderr, "\n\n"); 
    }	 
    /*free memory and close files*/
    /* free(rat);*/
      
    if(printfile==1){
      fclose(dh_fv);
      fclose(dinucleft_fv);
      fclose(dinucright_fv);
    }
    
    free(D);  
 }
  /*stefan mod.*/
	free(ancestorseq);			
  free(proselsite);
 

}


void MutateSequenceZqx(char *seq, double len, double **QT2ij){
 
   char    *ancestorseq;
   double  *rat;                                  
   double  rsum;
   double  t=0;
   int trip;
   int     secpos,i,j,k; 
   int     substitutions, osubstitutions;
   double  *probmuttriple;
   int     mutation;
   int     tripletmutation; 
   double  **D;
   char    dinucFilename[256];
   char    dinucRightFilename[256];
   char    dhFilename[256];
   double  sumdinuc;
   double  *proselsite;  /*stefan camb*/

     

    proselsite=(double *)malloc(seqlen*sizeof(double));  /*stefan camb*/
	ancestorseq =(char*)malloc(numBases*sizeof(char));
    for(i=0; i<numBases; i++)
         ancestorseq[i]=seq[i];
	rat=(double *)calloc(seqlen,sizeof(double));
    /*************Time*****************************************/
    rsum =0.0;	
	rsum = -seqlen;
	if(verbosewaiting==1) fprintf(stderr,"\n rsum1= %f\n",rsum);
	/*************Time*****************************************/
	t    = timeoverall(rsum);
	substitutions=0;
   /**********************************              
	        EVOLUTION PROCESS                        
    ********************************* */
    if(t>len && verboseWarning==1){
		fprintf(stderr, "No substitution! (waiting time:%f > branch length:%f)\n", t, len);
     }
	 /*hier also possible fixed D or poisson D (number of substitutions): Nick and Arndt are agree, that is not so proper*/
	 /*anyway do option*/
     while(t<=len)
     {
     secpos=-1;
     secpos=selsite_uniform();
      /*********************************/
	  /*subsitution and which*/
      mutation=-1;

      //      probmuttriple        = cprobmutTripletzqx(secpos, QT2ij, seq);

      //probmuttriple=(double *)malloc(64*sizeof(double));

      //for (k=0;k<64;k++) {
        // if (probmuttriple[k]<0.0001) continue;
        //      printf("%.3f ",probmuttriple[k]);
           //}
      
      //printf("\n");

      if (secpos==0){
        trip=tripletTO64(seq[seqlen-1],seq[secpos],seq[secpos+1]);
      } else {
        if(secpos==(seqlen-1)) {
          trip=tripletTO64(seq[secpos-1],seq[secpos],seq[0]);
        }
        else {
          trip=tripletTO64(seq[secpos-1],seq[secpos],seq[secpos+1]);
        }
      }

      probmuttriple=Pij[trip];


      //      for (k=0;k<64;k++) {
      //  if (Pij[trip][k]<0.0001) continue;
      //  printf("%.3f ",Pij[trip][k]);
      //}
      //      printf("\n\n");


	  tripletmutation      = cchoosetriplet(probmuttriple);
      mutation             = ctripletnucleotidmutation(tripletmutation, secpos, seq);
	  if(nucleotide[mutation] != seq[secpos]){
	     substitutions++;
	     seq=changeseq(seq,secpos,mutation);
         t += timeoverall(rsum);
	 } 
	 /********************************************************/
    }	
	free(rat);
	osubstitutions=0;
	for(i=0;i<numBases;i++){
		if(ancestorseq[i]!=seq[i]) osubstitutions=osubstitutions+1;
	}
	if(verboseEvolutiontime==1){
	  if(len==0.00) fprintf(stderr, "branch length: %f - no substituion\n", len);
	  else{
      if(printfile==1) {	
        sprintf(dhFilename, "%s_dh.sissi", filename);
        if((dh_fv=fopen(dhFilename,"a"))==NULL)
          fprintf(dh_fv, "dinuc.txt could not open\n");	
        fprintf(dh_fv ,  "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
      }else fprintf(stderr , "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
    }
	}
	/*print DiNucContent*/
	if((verboseDiNucContent==1 || (verboseDiNucContent==1 &&printfile==1)) && len>0){
	  /*allocate 4*4 matrix*/
	  D=(double**)malloc(4*sizeof(double*));
	  for(i=0;i<4;i++)
      D[i]=(double *)malloc(4*sizeof(double));
	  for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[i][j]=0.00;
	  }
    /*LEFT*/	
    for(i=1;i<numBases;i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i-1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/	
    if(printfile ==1 ){
      /* for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(stderr, "%c%c ", nucleotide[i],nucleotide[j]);*/
      sprintf(dinucFilename, "%s_dinuc_left.sissi", filename);
      if((dinucleft_fv=fopen(dinucFilename,"a"))==NULL)
	      fprintf(stderr, "dinuc.txt could not open\n");
      /* fprintf(dinucleft_fv,"%f \n", sumdinuc);	  
         for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(dinucleft_fv, "%c%c ", nucleotide[i],nucleotide[j]);
         fprintf(dinucleft_fv," \n");	*/
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucleft_fv,"%f ", D[j][i]);
        }	
      }	
      fprintf(dinucleft_fv, "\n");
	  } 
	  /*print*/
	  else{    
      for(i=0;i<4;i++){
		    for(j=0;j<4;j++) {
          fprintf(stderr ,"%f ", D[j][i]);
          fprintf(stderr, "%c%c ", nucleotide[j],nucleotide[i]);
        }   
      }
      fprintf(stderr, "\n"); 
	  }   	 
    /*Right*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]=0.00;
    }
    for(i=0;i<(numBases-1);i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i+1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/
    if(printfile==1){	
      sprintf(dinucRightFilename, "%s_dinuc_right.sissi", filename);
      if((dinucright_fv=fopen(dinucRightFilename,"a"))==NULL)
        fprintf(stderr, "dinuc.txt could not open\n");	  
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          fprintf(dinucright_fv,"%f ", D[j][i]);
        }	
      }	 
      fprintf(dinucright_fv, "\n");	
    }/*print */
    else{
      for(i=0;i<4;i++){
        for(j=0;j<4;j++) 
          fprintf(stderr,"%f ", D[j][i]);
      }	 
      fprintf(stderr, "\n\n"); 
    }	 
    /*free memory and close files*/
    /* free(rat);*/
      
    if(printfile==1){
      fclose(dh_fv);
      fclose(dinucleft_fv);
      fclose(dinucright_fv);
    }
    
    free(D);  
 }
  /*stefan mod.*/
	free(ancestorseq);			
  free(proselsite);
 

}



 
void MutateSequenceInd(char *seq, double len, double *Qij){
  char    *ancestorseq;
  double  *rat;                                  
  double  rsum;
  double  t=0;
  int     secpos,i,j;    
  int     substitutions, osubstitutions;
  double  *probmutij;
  int     mutation;
  double  **D;
  char dinucFilename[256];
  char dinucRightFilename[256];
  char dhFilename[256];
  double sumdinuc;
  double *proselsite;  /*stefan camb*/
     
  proselsite=(double *)malloc(seqlen*sizeof(double));  /*stefan camb*/
  ancestorseq =(char*)malloc(numBases*sizeof(char));
  for(i=0; i<numBases; i++) ancestorseq[i]=seq[i];
  rat=(double *)calloc(seqlen,sizeof(double));
  rate_independent(seq,Qij,rat);
  rsum = 0.0;	   
  rsum = rate_sum(rat);
  /*************Time*****************************************/
  t    = timeoverall(rsum);
  substitutions=0;
  /**********************************              
	        EVOLUTION PROCESS                        
          ********************************* */
  if(t>len && verboseWarning==1){
		fprintf(stderr, "No substitution! (waiting time:%f > branch length:%f)\n", t, len);
  }
  while(t<=len)
    {
      substitutions++; 
      secpos=-1;
	    secpos=selsite(rat,rsum,proselsite);
      /*********************************/
      mutation=-1;
      /*probmuttriple        = cprobmutTriplet(secpos, QTij,QT2ij, QT3ij, seq);*/
      probmutij= probmut(secpos,Qij,seq);
	    mutation=choosenuc(probmutij);
      free(probmutij);
      /*********************************/
	    seq=changeseq(seq,secpos,mutation);
      if(fastrates==1){
        /* fprintf(stderr, "rsum %lf ",rsum);*/
        rsum=minusrsumInd(secpos,rat, rsum);
        /*	fprintf(stderr, "  %lf",rsum);*/
        rat=changerateInd(secpos, rat, seq,Qij);
        rsum=plusrsumInd(secpos,rat, rsum);
        /*	fprintf(stderr, "  %lf\n",rsum);*/
      }	
      else{ 
        rat  = rate_independent(seq,Qij,rat);
        rsum= -1;
        rsum = rate_sum(rat);
      }	 
      /* if(verbosewaiting==1) fprintf(stderr, "\n rsum= %f\n",rsum);*/
      /* fprintf(stderr, "rsum= %f\n",rsum);*/
      t += timeoverall(rsum);
      /********************************************************/
    }	
	free(rat);
	/*free(probmuttriple);*/
	osubstitutions=0;
	for(i=0;i<numBases;i++){
		if(ancestorseq[i]!=seq[i]) osubstitutions=osubstitutions+1;
		/*if(P[i]!=seq[i]) osubstitutions=osubstitutions+1;*/
	}
	if(verboseEvolutiontime==1){
	  if(len==0.00) fprintf(stderr, "branch length: %f - no substituion\n", len);
	  else{
      if(printfile==1) {	
        sprintf(dhFilename, "%s_dh.sissi", filename);
        if((dh_fv=fopen(dhFilename,"a"))==NULL)
          fprintf(dh_fv, "dinuc.txt could not open\n");	
        fprintf(dh_fv ,  "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
      }else fprintf(stderr , "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
    }
	}
	/*print DiNucContent*/
	if((verboseDiNucContent==1 || (verboseDiNucContent==1 &&printfile==1)) && len>0){
	  /*allocate 4*4 matrix*/
	  D=(double**)malloc(4*sizeof(double*));
	  for(i=0;i<4;i++)
      D[i]=(double *)malloc(4*sizeof(double));
	  for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[i][j]=0.00;
	  }
    /*LEFT*/	
    for(i=1;i<numBases;i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i-1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/	
    if(printfile ==1 ){
      /* for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(stderr, "%c%c ", nucleotide[i],nucleotide[j]);*/
      sprintf(dinucFilename, "%s_dinuc_left.sissi", filename);
      if((dinucleft_fv=fopen(dinucFilename,"a"))==NULL)
	      fprintf(stderr, "dinuc.txt could not open\n");
      /* fprintf(dinucleft_fv,"%f \n", sumdinuc);	  
         for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(dinucleft_fv, "%c%c ", nucleotide[i],nucleotide[j]);
         fprintf(dinucleft_fv," \n");	*/
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucleft_fv,"%f ", D[j][i]);
        }	
      }	
      fprintf(dinucleft_fv, "\n");
	  } 
	  /*print*/
	  else{    
      for(i=0;i<4;i++){
		    for(j=0;j<4;j++) {
          fprintf(stderr ,"%f ", D[j][i]);
          fprintf(stderr, "%c%c ", nucleotide[j],nucleotide[i]);
        }   
      }
      fprintf(stderr, "\n"); 
	  }   	 
    /*Right*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]=0.00;
    }
    for(i=0;i<(numBases-1);i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i+1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/
    if(printfile==1){	
      sprintf(dinucRightFilename, "%s_dinuc_right.sissi", filename);
      if((dinucright_fv=fopen(dinucRightFilename,"a"))==NULL)
        fprintf(stderr, "dinuc.txt could not open\n");	  
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucright_fv,"%f ", D[j][i]);
        }	
	    }	 
	    fprintf(dinucright_fv, "\n");	
	  }/*print */
	  else{
      for(i=0;i<4;i++){
        for(j=0;j<4;j++) 
		      fprintf(stderr,"%f ", D[j][i]);
		  }	 
		  fprintf(stderr, "\n\n"); 
	  }	 
    /*free memory and close files*/
    /* free(rat);*/
	  if(printfile==1){
      fclose(dh_fv);
      fclose(dinucleft_fv);
      fclose(dinucright_fv);
    }
    
    free(D);  

  }
	free(ancestorseq);			
  free(proselsite);

}

 
 void MutateSequenceIndqx(char *seq, double len, double *Qij){
  char    *ancestorseq;
  double  *rat;                                  
  double  rsum;
  double  t=0;
  int     secpos,i,j;    
  int     substitutions, osubstitutions;
  double  *probmutij;
  int     mutation;
  double  **D;
  char dinucFilename[256];
  char dinucRightFilename[256];
  char dhFilename[256];
  double sumdinuc;
  double *proselsite;  /*stefan camb*/
     
  proselsite=(double *)malloc(seqlen*sizeof(double));  /*stefan camb*/
  ancestorseq =(char*)malloc(numBases*sizeof(char));
  for(i=0; i<numBases; i++) ancestorseq[i]=seq[i];
  rat=(double *)calloc(seqlen,sizeof(double));
  rate_independent(seq,Qij,rat);
  rsum = 0.0;	   
 /* rsum = rate_sum(rat);*/
  rsum=-seqlen;
  /*************Time*****************************************/
  t    = timeoverall(rsum);
  substitutions=0;
  /**********************************              
	        EVOLUTION PROCESS                        
          ********************************* */
  if(t>len && verboseWarning==1){
		fprintf(stderr, "No substitution! (waiting time:%f > branch length:%f)\n", t, len);
  }
  while(t<=len)
    {
      secpos=-1;
	  secpos=selsite_uniform();
      /*********************************/
      mutation=-1;
	  probmutij=probmutqx(secpos,Qij,seq);
      mutation=choosenuc(probmutij);
	  if(nucleotide[mutation] != seq[secpos]){
	    substitutions++;
		seq=changeseq(seq,secpos,mutation);
		t+=timeoverall(rsum);
	  }	
      free(probmutij);
      /*********************************/
	 /* seq=changeseq(seq,secpos,mutation);*/
    /*  if(fastrates==1){*/
        /* fprintf(stderr, "rsum %lf ",rsum);*/
      /*  rsum=minusrsumInd(secpos,rat, rsum);*/
        /*	fprintf(stderr, "  %lf",rsum);*/
     /*   rat=changerateInd(secpos, rat, seq,Qij);
        rsum=plusrsumInd(secpos,rat, rsum);*/
        /*	fprintf(stderr, "  %lf\n",rsum);*/
    /*  }	
      else{ 
        rat  = rate_independent(seq,Qij,rat);
        rsum= -1;
        rsum = rate_sum(rat);
      }	 */
      /* if(verbosewaiting==1) fprintf(stderr, "\n rsum= %f\n",rsum);*/
      /* fprintf(stderr, "rsum= %f\n",rsum);*/
    /*  t += timeoverall(rsum);*/
      /********************************************************/
    }	
	free(rat);
	/*free(probmuttriple);*/
	osubstitutions=0;
	for(i=0;i<numBases;i++){
		if(ancestorseq[i]!=seq[i]) osubstitutions=osubstitutions+1;
		/*if(P[i]!=seq[i]) osubstitutions=osubstitutions+1;*/
	}
	if(verboseEvolutiontime==1){
	  if(len==0.00) fprintf(stderr, "branch length: %f - no substituion\n", len);
	  else{
      if(printfile==1) {	
        sprintf(dhFilename, "%s_dh.sissi", filename);
        if((dh_fv=fopen(dhFilename,"a"))==NULL)
          fprintf(dh_fv, "dinuc.txt could not open\n");	
        fprintf(dh_fv ,  "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
      }else fprintf(stderr , "branch length: %f \t substitutions: %6d %10f \t Hammingdistance: %6d %10f \n", len,substitutions,(double)substitutions/seqlen,osubstitutions,(double)osubstitutions/seqlen);
    }
	}
	/*print DiNucContent*/
	if((verboseDiNucContent==1 || (verboseDiNucContent==1 &&printfile==1)) && len>0){
	  /*allocate 4*4 matrix*/
	  D=(double**)malloc(4*sizeof(double*));
	  for(i=0;i<4;i++)
      D[i]=(double *)malloc(4*sizeof(double));
	  for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[i][j]=0.00;
	  }
    /*LEFT*/	
    for(i=1;i<numBases;i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i-1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/	
    if(printfile ==1 ){
      /* for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(stderr, "%c%c ", nucleotide[i],nucleotide[j]);*/
      sprintf(dinucFilename, "%s_dinuc_left.sissi", filename);
      if((dinucleft_fv=fopen(dinucFilename,"a"))==NULL)
	      fprintf(stderr, "dinuc.txt could not open\n");
      /* fprintf(dinucleft_fv,"%f \n", sumdinuc);	  
         for(i=0;i<4;i++)
         for(j=0;j<4;j++) fprintf(dinucleft_fv, "%c%c ", nucleotide[i],nucleotide[j]);
         fprintf(dinucleft_fv," \n");	*/
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucleft_fv,"%f ", D[j][i]);
        }	
      }	
      fprintf(dinucleft_fv, "\n");
	  } 
	  /*print*/
	  else{    
      for(i=0;i<4;i++){
		    for(j=0;j<4;j++) {
          fprintf(stderr ,"%f ", D[j][i]);
          fprintf(stderr, "%c%c ", nucleotide[j],nucleotide[i]);
        }   
      }
      fprintf(stderr, "\n"); 
	  }   	 
    /*Right*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]=0.00;
    }
    for(i=0;i<(numBases-1);i++) D[nucleoTO123(seq[i])][nucleoTO123(seq[i+1])]++;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++) D[j][i]/=(numBases-1);}
    sumdinuc=0.00;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++) sumdinuc += D[j][i];
    if(sumdinuc !=1.00){
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++) D[j][i]/=sumdinuc;}
    }	  
    /*print in File*/
    if(printfile==1){	
      sprintf(dinucRightFilename, "%s_dinuc_right.sissi", filename);
      if((dinucright_fv=fopen(dinucRightFilename,"a"))==NULL)
        fprintf(stderr, "dinuc.txt could not open\n");	  
      for(i=0;i<4;i++){
	      for(j=0;j<4;j++){
          fprintf(dinucright_fv,"%f ", D[j][i]);
        }	
	    }	 
	    fprintf(dinucright_fv, "\n");	
	  }/*print */
	  else{
      for(i=0;i<4;i++){
        for(j=0;j<4;j++) 
		      fprintf(stderr,"%f ", D[j][i]);
		  }	 
		  fprintf(stderr, "\n\n"); 
	  }	 
    /*free memory and close files*/
    /* free(rat);*/
	  if(printfile==1){
      fclose(dh_fv);
      fclose(dinucleft_fv);
      fclose(dinucright_fv);
    }
    
    free(D);  

  }
	free(ancestorseq);			
  free(proselsite);

}


 
 
 
 
 
 /****************************************************************************************************************/
 /*                              Choose a site k                                                                 */
 /****************************************************************************************************************/
 /**********************************/
 /*  DEPENDENT = INDEPENDENT CASE  */
 /**********************************/
/* double *selsite(double rat[], double rsum, double *proselsite){
       int pos;
    
       for(pos=0;pos<seqlen;pos++){
	        if(sitefactor==0) proselsite[pos]=rat[pos]/rsum ;
			else  proselsite[pos]= (scalingfactor[pos]*rat[pos])/rsum ;
		}
      
}*/

int selsite_uniform(){
  
  int m, secpos;
  double dummy=0.00, x=0.00;
  double a;
  int upper, lower, n;
	
#  if SPRNG
  if(randomgeneratordkiss == 1)
#  endif    
    x=dkiss();               
#  if SPRNG	   
  else	
    x=sprng(randstream);  
#  endif  

  /*
    m=0;
    secpos=-1;
    do
    {  
      secpos=m;
      if(sitefactor==0) dummy += (1.00/(seqlen));
      else dummy += (1.00/(seqlen))*scalingfactor[m];	 
      m++;
    }
  while(x>dummy && m<seqlen);
  */

   
  x*=seqlen;
  if (sitefactor==0){
    secpos=ceil(x)-1;
  } else {
    secpos=0;
    if (x>=scalingfactor_sum[0]){
      upper=seqlen-1;
      lower=0;
      while (upper-lower>1){
        n=lower+floor((upper-lower+1)/2);
        if (scalingfactor_sum[n]>x) {
          upper=n;
        } else {
          lower=n;
        }
      }
      secpos=upper;
    }
  }
  
  return(secpos);	
}  


inline int selsite(double rat[], double rsum, double *proselsite){
       int pos;
	    double x;
		double dummy=0.00;
		int m;
		int secpos;
		double proselsitesum;
		
		/*for(pos=0;pos<seqlen;pos++){
	        if(sitefactor==0) proselsite[pos]=rat[pos]/rsum ;
			else  proselsite[pos]= (scalingfactor[pos]*rat[pos])/rsum ;
		}*/
		
		for(pos=0;pos<seqlen;pos++){
		    proselsite[pos]=rat[pos]/rsum ;
		/*   fprintf(stderr, "rat %lf rsum %lf = %lf\n", rat[pos], rsum, proselsite[pos]);*/
		}
		

		proselsitesum=0.00;
		for(pos=0;pos<seqlen;pos++){
		/* fprintf(stderr, "pos %d proselsite %lf \n",pos, proselsite[pos]);*/
		 proselsitesum += proselsite[pos];
		 }
		/* fprintf(stderr, "%lf \n", proselsitesum);*/
		 
		 if((proselsitesum < 1.00000 - NEARZERO) || (proselsitesum > 1.00000 + NEARZERO)){
			fprintf(stderr, "\n  Normalize prob. %lf to pick site k with rsum %lf \n", proselsitesum, rsum );
		     for(pos=0;pos<seqlen;pos++){
			   proselsite[pos]= proselsite[pos]/proselsitesum ;
			 }
			  
		 }	 
		 
		 proselsitesum=0.00;
		for(pos=0;pos<seqlen;pos++){
		/* fprintf(stderr, "pos %d proselsite %lf \n",pos, proselsite[pos]);*/
		 proselsitesum += proselsite[pos];
		 }
		/* fprintf(stderr, "\nnew %lf \n", proselsitesum);*/
         		
		
		
		/* x=getrandomdouble();*/
		/* x=(double) rand()/RAND_MAX;*/
		 #  if SPRNG
        if(randomgeneratordkiss == 1)
        #  endif    
	    x=dkiss();               
        #  if SPRNG	   
        else	
        x=sprng(randstream);  
        #  endif  
		
		 m=0;
		 secpos=-1;

         do
         {
             secpos=m;
             dummy += proselsite[m];
             m++;
         }
        /* while(x>dummy && m<seqlen);*/
        while(x>dummy && m<seqlen);
		
		
   return(secpos);		
       
}
  /************************************************************/
  /* RANDOM CHOOSE: Select one positon for the next mutation  */
  /************************************************************/
   int rand_choose(double *proselsite){
         double x;
         double dummy=0;
         int m;
         int secpos;

         x=getrandomdouble();
         m=0;
		 secpos=-1;

         do
         {
             secpos=m;
             dummy += proselsite[m];
             m++;
         }
         while(x>dummy && m<seqlen);
       /* free(proselsite);*/
         return(secpos);
    }
   /****************************************************************************************************************/
   /*                                Replace xk with yo                                                            */
   /****************************************************************************************************************/
    /***************************************/
	/*    INDEPENDENT CASE                 */
	/***************************************/       
    double *probmut( int secpos, double *Qij, char *seq){
         int nucsecpos,j ;
         double *probmutij;
		 
         nucsecpos=-1;
         probmutij =(double *)malloc(4*sizeof(double));
         nucsecpos= nucleoTO123(seq[secpos]);

         for(j=0;j<4;j++) probmutij[j]=0.00;
         for(j=0;j<4;j++){
            if(nucsecpos != j){
		      probmutij[j]= -Qij[nucsecpos*4+j]/Qij[nucsecpos*4+nucsecpos];
            }
         }
         return(probmutij);
     }
	  double *probmutqx( int secpos, double *Qij, char *seq){
         int    nucsecpos,j ;
         double *probmutij;
		 double maxwait=0.00;
		 double probsum=0.00;
		 
         nucsecpos=-1;
         probmutij =(double *)malloc(4*sizeof(double));
         nucsecpos= nucleoTO123(seq[secpos]);
		 
		 for(j=0;j<4;j++) {
		   if(Qij[j*4+j]<maxwait) maxwait=Qij[j*4+j];
		 }
		 maxwait=-maxwait;
		 probsum=0.00;
		 for(j=0;j<4;j++) probmutij[j]=0.00;
		 
		 for(j=0;j<4;j++) {
		   if(nucsecpos != j){
		       probmutij[j]= Qij[nucsecpos*4+j]/maxwait;
		   } else  probmutij[j]= (maxwait+Qij[nucsecpos*4+j])/maxwait;
		   probsum += probmutij[j];
		 }
		 if(verbose==1) fprintf(stderr, "sum %f \n", probsum);  

         
         return(probmutij);
     }

     /*CHOOSE PER RANDOM NUCLEOTID J */
	 int choosenuc(double *probmutij){
           int j, mutation;
           double dummy, x;

           x=getrandomdouble();
           mutation=0;
           j=0;
           dummy=0.00;
           do{
              dummy += probmutij[j];
              mutation=j;
              ++j;
            }while(x>=dummy);

            return(mutation);
	 }
/*******************************************************/
/*              nk=2                                   */
/*******************************************************/
/*******************************************************/
/*
int choosetriplet(double *probmutij_trip){
           int j ;
           double dummy , x;
           int tripletmutation;


           x=getrandomdouble();
           tripletmutation=0;
           j=0;
           dummy=0.00;

           do{
             dummy =  dummy + probmutij_trip[j];
             if(probmutij_trip[j] != 0) tripletmutation=j;
             ++j;
           }while(x>=dummy ) ;
		   free(probmutij_trip);
           return (tripletmutation);
}
*/
/*******************************************************/
/* Triplet  with categorie                            */
/*******************************************************/
double *cprobmutTripletz(int pos, double **QT2ij, char *seq){
           int       j,t=-1 ;
           double    *probmutij_trip;
           double    probsum=0.00;
		  /* int       secposneigh1=-1;
		   int       secposneigh2=-1; */
		   
		    probmutij_trip=(double *)malloc(64*sizeof(double));
			for(j=0 ; j<64 ; j++) probmutij_trip[j]=0.00;

            if(pos==0) t=tripletTO64(seq[seqlen-1],seq[pos],seq[pos+1]);
			else if(pos==(seqlen-1)) t=tripletTO64(seq[pos-1],seq[pos],seq[0]);
            else t=tripletTO64(seq[pos-1],seq[pos],seq[pos+1]);
			for(j=0 ; j<64; j++) {
				if(t != j) {
						probmutij_trip[j] = - QT2ij[t][j]/QT2ij[t][t];
						probsum += probmutij_trip[j] ;
				}
			}

				 /*  if(Nseq.item[pos].cat_pos==0) fprintf(stderr, "\nPos %d Cateforie 0\n", pos); 
                  if(verbosewaiting==1) fprintf(stderr, "triple %d:%s with", t, triplet[t]);
                  if(verbosewaiting==1) fprintf(stderr, "sum %f->", probsum);*/
           return(probmutij_trip);
}


double** computePij(double **QT2ij){

  int i,k;
  int       j,t=-1 ;
  double   *probmutij_trip,probsum_mut,waiting_qjj, norm ;
  double    probsum=0.00;
  double    maxwait=0.00;

  Pij=(double**)malloc(64*sizeof(double*));

  for (t=0;t<64;t++){

    probmutij_trip=(double *)malloc(64*sizeof(double));

    for(j=0 ; j<64 ; j++) probmutij_trip[j]=0.00;
			
    probsum_mut=0.00;
    waiting_qjj=0.00;
    norm=0.00;
			
    for(j=0 ; j<64; j++) {
      if(QT2ij[j][j]<maxwait) maxwait=QT2ij[j][j];
    }	

		maxwait=-maxwait;

    probsum=0.00;
    for(j=0 ; j<64; j++) {
      if(t != j) 
        probmutij_trip[j] =  QT2ij[t][j]/maxwait;
      else	probmutij_trip[j] =  (maxwait + QT2ij[t][j])/maxwait ;
    }

    Pij[t]=probmutij_trip;

    /*

    for (k=0;k<64;k++) {
      if (probmutij_trip[k]<0.0001) continue;
      printf("%.3f ",probmutij_trip[k]);
    }

    printf("\n");

    */

  }

  return Pij;

}		

double *cprobmutTripletzqx(int pos, double **QT2ij, char *seq){
  int       j,t=-1 ;
  double    *probmutij_trip,probsum_mut,waiting_qjj, norm ;
  double    probsum=0.00;
  double    maxwait=0.00;
		   
		   
		    probmutij_trip=(double *)malloc(64*sizeof(double));
			for(j=0 ; j<64 ; j++) probmutij_trip[j]=0.00;

            if(pos==0) t=tripletTO64(seq[seqlen-1],seq[pos],seq[pos+1]);
			else if(pos==(seqlen-1)) t=tripletTO64(seq[pos-1],seq[pos],seq[0]);
            else t=tripletTO64(seq[pos-1],seq[pos],seq[pos+1]);
			
			probsum_mut=0.00;
			waiting_qjj=0.00;
			norm=0.00;
			
			for(j=0 ; j<64; j++) {
			    if(QT2ij[j][j]<maxwait) maxwait=QT2ij[j][j];
			}	
			maxwait=-maxwait;
			/*fprintf(stderr, " max %lf-> ", maxwait);*/
			probsum=0.00;
		 	for(j=0 ; j<64; j++) {
			   if(t != j) 
					probmutij_trip[j] =  QT2ij[t][j]/maxwait;
			   else	probmutij_trip[j] =  (maxwait + QT2ij[t][j])/maxwait ;
			  /* if(probmutij_trip[j]!=0.00) fprintf(stderr, "%d=%lf ",j, probmutij_trip[j] );*/
			 /*  if( t == j) fprintf(stderr, "%lf -> " , probmutij_trip[j] );*/
         /*          probsum += probmutij_trip[j] ;*/
			}

			/*fprintf(stderr, "waiting waiting_qjj %lf, probsum_mut %lf , sum%lf  ",waiting_qjj, probsum_mut, probsum_mut +waiting_qjj);*/
			/*fprintf(stderr, "triple %d:%s with", t, triplet[t]);
			fprintf(stderr, "sum %f->\n", probsum);*/
          
      /*for(j=0 ; j<64; j++) norm = norm+ (tripi[j]*(maxwait+(QT2ij[j][j]/maxwait)));*/
			/*fprintf(stderr, "number of substitutions per site and unit time: %lf \n " ,norm);*/		
				
           return(probmutij_trip);
}



/***************************************/
int cchoosetriplet(double *probmutij_trip){
           int j ;
           double dummy , x;
           int tripletmutation;


           x=getrandomdouble();
           tripletmutation=0;
           j=0;
           dummy=0.00;

           do{
             dummy =  dummy + probmutij_trip[j];
             if(probmutij_trip[j] != 0) tripletmutation=j;
             ++j;
           }while(x>=dummy ) ;
		   
           /*free(probmutij_trip);*/
		   
           return (tripletmutation);
         }


          /*AND GIVE THE NUCLEOTIDCHANGEMUTATION J */
          int ctripletnucleotidmutation(int tripletmutation, int secpos, char *seq)
          {
            int  ii, iii;
            int mutation=-1;
			
			
			ii =(tripletmutation+4)/4;
			iii=(tripletmutation+16)/16;
		
			 mutation=(ii-1)%4; 
			return(mutation);
}
/****************************************************************************************************************/
/* CALCULATE THE SPONTANEOUS SUBSTITUTIONSRATE FOR THE POSITIONS                                                */
/****************************************************************************************************************/
double *rateDIz(char *seq, double **QT2ij, double *rat)
{
	/*double *rat ;*/
	int pos,t;

	/*rat=(double *)calloc(seqlen,sizeof(double));*/
	for(pos=0; pos<seqlen;pos++){
	    rat[pos]=10;
	    /* fprintf(stderr, "%d von %d\n",pos, seqlen);*/
		if(pos==0){
              t=tripletTO64(seq[seqlen-1],seq[0],seq[1]);
		}else if(pos==seqlen-1){	 
		      t=tripletTO64(seq[seqlen-2],seq[seqlen-1],seq[0]);
	    }else {
	          t=tripletTO64(seq[pos-1],seq[pos],seq[pos+1]);
	    }
		
		if(sitefactor==0)  rat[pos]= QT2ij[t][t];
	    else               rat[pos]= QT2ij[t][t]*scalingfactor[pos];
	  
	/*  fprintf(stderr,"%d, %lf\n", pos, rat[pos]);*/
			 

  /* stefan mod.*/
  /*	if(rat[pos]==10) {printf("WRONG RATE1\n") ; 
      fprintf(stderr, "on position %d", pos);
      exit (1);
    }
  */
  }
    return(rat);
}


double *changerateDIz(int secpos, double *rat, char *seq, double **QT2ij){
  int t1,t2,t3;
 /* fprintf(stderr, "pos %d", secpos);*/
  if(secpos==0){
  t1=tripletTO64(seq[seqlen-2],seq[seqlen-1],seq[secpos]);
  t2=tripletTO64(seq[seqlen-1],seq[secpos],seq[secpos+1]);
  t3=tripletTO64(seq[secpos],seq[secpos+1],seq[secpos+2]);
  if(sitefactor==0){  
     rat[seqlen-1]=QT2ij[t1][t1];
     rat[secpos]=QT2ij[t2][t2];
     rat[secpos+1]=QT2ij[t3][t3];
  }else{
     rat[seqlen-1]=QT2ij[t1][t1]*scalingfactor[seqlen-1];
     rat[secpos]=QT2ij[t2][t2]*scalingfactor[secpos];
     rat[secpos+1]=QT2ij[t3][t3]*scalingfactor[secpos+1];
  }	 
}  
else if(secpos==1){
  t1=tripletTO64(seq[seqlen-1],seq[secpos-1],seq[secpos]);
  t2=tripletTO64(seq[secpos-1],seq[secpos],seq[secpos+1]);
  t3=tripletTO64(seq[secpos],seq[secpos+1],seq[secpos+2]);
  if(sitefactor==0){ 
     rat[secpos-1]=QT2ij[t1][t1];
     rat[secpos]=QT2ij[t2][t2];
     rat[secpos+1]=QT2ij[t3][t3];
  }else{
	 rat[secpos-1]=QT2ij[t1][t1]*scalingfactor[secpos-1];
     rat[secpos]=QT2ij[t2][t2]*scalingfactor[secpos];
     rat[secpos+1]=QT2ij[t3][t3]*scalingfactor[secpos+1];
 }	 
}  
else if(secpos==seqlen-1){
  t1=tripletTO64(seq[secpos-2],seq[secpos-1],seq[secpos]);
  t2=tripletTO64(seq[secpos-1],seq[secpos],seq[0]);
  t3=tripletTO64(seq[secpos],seq[0],seq[1]);
  if(sitefactor==0){ 
     rat[secpos-1]=QT2ij[t1][t1];
     rat[secpos]=QT2ij[t2][t2];
     rat[0]=QT2ij[t3][t3];
  }else{
  	 rat[secpos-1]=QT2ij[t1][t1]*scalingfactor[secpos-1];
     rat[secpos]=QT2ij[t2][t2]*scalingfactor[secpos];
     rat[0]=QT2ij[t3][t3]*scalingfactor[0];
  }	 
} 
else if(secpos==seqlen-2){
  t1=tripletTO64(seq[secpos-2],seq[secpos-1],seq[secpos]);
  t2=tripletTO64(seq[secpos-1],seq[secpos],seq[secpos+1]);
  t3=tripletTO64(seq[secpos],seq[secpos+1],seq[0]);
  if(sitefactor==0){ 
     rat[secpos-1]=QT2ij[t1][t1];
     rat[secpos]=QT2ij[t2][t2];
     rat[secpos+1]=QT2ij[t3][t3];
  }else{
      rat[secpos-1]=QT2ij[t1][t1]*scalingfactor[secpos-1];
      rat[secpos]=QT2ij[t2][t2]*scalingfactor[secpos];
      rat[secpos+1]=QT2ij[t3][t3]*scalingfactor[secpos+1];
 }	   
} 
else{
  t1=tripletTO64(seq[secpos-2],seq[secpos-1],seq[secpos]);
  t2=tripletTO64(seq[secpos-1],seq[secpos],seq[secpos+1]);
  t3=tripletTO64(seq[secpos],seq[secpos+1],seq[secpos+2]);
  if(sitefactor==0){ 
     rat[secpos-1]=QT2ij[t1][t1];
     rat[secpos]=QT2ij[t2][t2];
     rat[secpos+1]=QT2ij[t3][t3];
 }else{
     rat[secpos-1]=QT2ij[t1][t1]*scalingfactor[secpos-1];
     rat[secpos]=QT2ij[t2][t2]*scalingfactor[secpos];
     rat[secpos+1]=QT2ij[t3][t3]*scalingfactor[secpos+1];
}	 
 	 
}

   return(rat);
}



double minusrsumDIz(int secpos, double *rat, double rsum){
  if(secpos==0){
    rsum=rsum- rat[seqlen-1];
    rsum=rsum- rat[secpos];
    rsum=rsum- rat[secpos+1];
}  
else if(secpos==seqlen-1){
 rsum=rsum-  rat[secpos-1];
 rsum=rsum-  rat[secpos];
 rsum=rsum-  rat[0];
} 
else{
 rsum=rsum- rat[secpos-1];
 rsum=rsum-  rat[secpos];
 rsum=rsum-  rat[secpos+1];
}

   return(rsum);
}

double plusrsumDIz(int secpos, double *rat, double rsum){
  if(secpos==0){
    rsum=rsum+ rat[seqlen-1];
    rsum=rsum+ rat[secpos];
    rsum=rsum+ rat[secpos+1];
}  
else if(secpos==seqlen-1){
 rsum=rsum+  rat[secpos-1];
 rsum=rsum+  rat[secpos];
 rsum=rsum+  rat[0];
} 
else{
 rsum=rsum+ rat[secpos-1];
 rsum=rsum+  rat[secpos];
 rsum=rsum+  rat[secpos+1];
}



   return(rsum);
}




double *rate_independent(char *seq, double *Qij,double *rat){
  /*           double *rat ;*/
          /* double addrat=0;*/
           int pos,w;

           /*rat=(double *)calloc(seqlen,sizeof(double));*/
		   w=0;
		   
		   for(pos=0; pos<seqlen;pos++) rat[pos]=-10;
		   
		   
		/*   if(sitefactor==0)  rat[pos]= QT2ij[t][t];
	    else               rat[pos]= QT2ij[t][t]*scalingfactor[pos];*/
		   
		   
		   for(pos=0; pos<seqlen;pos++){
		       if(sitefactor==0) rat[pos]=Qij[nucleoTO123(seq[pos])*4+nucleoTO123(seq[pos])];
			   else rat[pos]=(Qij[nucleoTO123(seq[pos])*4+nucleoTO123(seq[pos])])*scalingfactor[pos];

           if(rat[pos]==-10) {printf("WRONG RATE2\n") ; exit (1);}

		   }		 

          /*  if(verbosewaiting==1) addrat += rat[pos];*/
		  /* if(verbosewaiting==1) fprintf(stderr,"summierte Raten: %f\n",addrat);*/
		  /* fprintf(stderr"\n addrat=%f\n", addrat);*/
          return(rat);
      
 }
 /****************************************************************************************************************/
 /*    CALCULATE THE SPONTANEOUS SUBSTITUTIONSRATE OVER ALL POSITIONS  q(x)                                       */
 /****************************************************************************************************************/
  /**********************************/
  /*  DEPENDENT = INDEPENDENT CASE  */
  /**********************************/
  double rate_sum(double rat[]){
         int lpos;
         double rsum;
		 char qxFilename[256];

         rsum=0.0;
		  for(lpos=0; lpos<seqlen; lpos++) rsum += rat[lpos];
		 
		/*if(sitefactor==0)
             for(lpos=0; lpos<seqlen; lpos++){ rsum += rat[lpos];
											    }
		else for(lpos=0; lpos<seqlen; lpos++){ rsum += (rat[lpos]*scalingfactor[lpos]);
											   }*/
			
		if(scalesequencerate==1){
		      fprintf(stderr, " Sequence rate with %lf\n", sqx);
			  rsum=rsum*sqx;
		}	  						
			
		/* rsum = rsum/seqlen;	*/
        if(verbosewaiting==1){
		     if(printfile==1) {
		          sprintf(qxFilename, "%s_qx.sissi", filename);
		          if((qx_fv=fopen(qxFilename,"a"))==NULL)
		    	  fprintf(qx_fv, "qx.txt could not open with option -wp\n");	
				  if(verbosewaitingcount==0) fprintf(qx_fv, "\n") ;  
				   fprintf(qx_fv, " %.4f ", rsum);
		     }else fprintf(stderr, " %.4f ", rsum);
		    
	   }	
	   

	   
		 if(printfile==1 && verbosewaiting==1) fclose(qx_fv);
		 verbosewaitingcount=1;	 /*only one time printed out */
		 
		
         return(rsum);

   }

/***************************************************************************/ 
double *changerateInd(int secpos, double *rat, char *seq, double *Qij){
   if(sitefactor==0)  rat[secpos]=(Qij[nucleoTO123(seq[secpos])*4+nucleoTO123(seq[secpos])]);
   else rat[secpos]=(Qij[nucleoTO123(seq[secpos])*4+nucleoTO123(seq[secpos])])*(scalingfactor[secpos]);
   return(rat);
}
double minusrsumInd(int secpos, double *rat, double rsum){
   rsum=rsum-rat[secpos];
   return(rsum);
}
double plusrsumInd(int secpos, double *rat, double rsum){
   rsum=rsum+  rat[secpos];
   return(rsum);
}
/***************************************************************************/ 

 /*
  * Random Sequence indpendent
  */		
  void RandomSequence_independent(char *seq)
  {
	char *P, *rseq;
	double p[4], x;
    int i, nr, ink, pos;

	rseq=(char*)calloc(seqlen,sizeof(char));
	P=seq;
	
	   p[0]=0;
       p[1]=mpi[0];
       p[2]=mpi[0]+mpi[1];
       p[3]=mpi[0]+mpi[1]+mpi[2];

	for(pos=0; pos<seqlen; pos++)
	{
			do{
			   nr=0;
			   x=getrandomdouble();
			   for(ink=1; ink<4;ink++)
				  if(x>p[ink]) nr=ink;
			}while(nr<0 || nr>3);

			rseq[pos]=nucleotide[nr];
	}
           	
	for (i=0; i<seqlen; i++) {
		*P=rseq[i];
		P++;
	}
	if(verboserandom==1) for(pos=0; pos<seqlen; pos++) fprintf(stderr, "%c", rseq[pos]);

  free(rseq);

 }

 
	

