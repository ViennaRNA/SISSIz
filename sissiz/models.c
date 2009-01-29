 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <ctype.h>

 #include "models.h"
 #include "global.h"
 #include "mutate.h"
 
 #define NEARZERO 1.0e-10



 double mpi[4];             /* stationary nucleotide distribution [agct]*/
 double dipi[16];           /* stationary dinucleotide distribution     */
 double tripi[64];          /*  stationary trinucleotide distribution   */
 int    matrixSH;
 int    balance;
 double Rmat[6];            /*  6 free rate parameters*/
 /*
  double Rmat16[24];
 double Rmat2T[96];
 */
 double trtv;
 double ratio;             
 double quo;  
 int    transitiontransverion,prmat, prmat2t, prmat16;
 int    simind, simindependent;
  double K;
 /*
  int    mononucneigh2, mononucneigh1;
 int    CpG;
 double CpGTpG, CpGCpA;
*/
 
 double ratio;              /*  apha/beta in the model of HKY*/
 double quo;                /*  ratio expected #transitions/#transversions*/
 int    mix;
 int    cpg;
 int     verboseFreqOrder, verboseCluster;
 int    readfreq;
 

 
 
 char   nucleotide[4];       /*  nucleotide:m RNA:{A,C,G,U}, with optin -d: DNA:{A,C,G,T}*/
/* char   nucleotide[]={"A","C","G","U"};*/
 /*char   *doublet[];
 char   *triplet[];*/
 char   *doublet[]={"AA","AC","AG","AU","CA","CC","CG","CU","GA","GC","GG","GU","UA","UC","UG","UU"};
 char   *triplet[]={"AAA","AAC","AAG","AAU","ACA","ACC","ACG","ACU","AGA","AGC","AGG","AGU","AUA","AUC","AUG","AUU",
					"CAA","CAC","CAG","CAU","CCA","CCC","CCG","CCU","CGA","CGC","CGG","CGU","CUA","CUC","CUG","CUU",
					"GAA","GAC","GAG","GAU","GCA","GCC","GCG","GCU","GGA","GGC","GGG","GGU","GUA","GUC","GUG","GUU",
					"UAA","UAC","UAG","UAU","UCA","UCC","UCG","UCU","UGA","UGC","UGG","UGU","UUA","UUC","UUG","UUU",};
  
 char freqFilename[256];
 char matrFilename[256];
 FILE *freq_fv;
 FILE *matr_fv;
 
 
 /*functions*/
 /****************************************************************************************************************/
 /*                COMPUTE FREQUENCIES                                                                        */
 /****************************************************************************************************************/
 /***************************************/
 /*  CHECK FREQUENCIES                  */
 /***************************************/
 void checkFreq(){
    int i;
	double sumdipi, sumpi;
	double sumtripi;
 
 
    /*verboseorder*/
	if(verboseFreqOrder==1){
	  fprintf(stderr,"\n");
	  for(i=0;i<4;i++) fprintf(stderr, "%c ", nucleotide[i]);
	  fprintf(stderr,"\n");
	  for(i=0;i<16;i++) fprintf(stderr, "%s ", doublet[i]);
	  fprintf(stderr,"\n");
	  for(i=0;i<64;i++) fprintf(stderr, "%s ", triplet[i]);
	  fprintf(stderr,"\n");
	}

    /*********************/
	/*MONO*/
	sumpi=0;
	for(i=0;i<4;i++) sumpi += mpi[i];
	if((sumpi < 1.00000 - NEARZERO) || (sumpi > 1.00000 + NEARZERO)){
	  for(i=0;i<4;i++) mpi[i] = mpi[i]/sumpi;
	  if(verboseFrequencies==1){
	    fprintf(stderr,"New\n pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
	    fprintf(stderr,"sum pi(n_k=0)%f\n",sumpi);
	  }	
     }
	/*********************/
	/*DI*/
     sumdipi=0;
	 for(i=0;i<16;i++)
	 sumdipi+= dipi[i];	
	 if(verboseFrequencies==1){
	   fprintf(stderr, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", dipi[0], dipi[1], dipi[2], dipi[3], dipi[4], dipi[5], dipi[6], dipi[7],dipi[8], dipi[9],dipi[10],dipi[11], dipi[12],dipi[13],dipi[14],dipi[15]);
       fprintf(stderr, "sum pi(nk=1) %f\n", sumdipi);
	 }  
	 if((sumdipi < 1.00000 - NEARZERO) || (sumdipi > 1.00000 + NEARZERO)){
	    if(verboseFrequencies==1) fprintf(stderr, "1.New Doublet_freq :\n");
		for(i=0;i<16;i++) dipi[i] = dipi[i]/sumdipi;
		sumdipi=0;	
		for(i=0;i<16;i++)   sumdipi+= dipi[i];	
		if(verboseFrequencies==1 && verboseFreqOrder==1){
		  if(nucleotide[3]=='t' || nucleotide[3]=='T'){
             fprintf(stderr, "frequencies:aa:%f\t ac:%f\t ag:%f\t at:%f\t\n",  dipi[0], dipi[1], dipi[2], dipi[3]);
             fprintf(stderr, "           :ca:%f\t cc:%f\t cg:%f\t ct:%f\t\n",  dipi[4], dipi[5], dipi[6], dipi[7]);
			 fprintf(stderr, "           :ga:%f\t gc:%f\t gg:%f\t gt:%f\t\n",  dipi[8], dipi[9],dipi[10],dipi[11]);
             fprintf(stderr, "           :ta:%f\t tc:%f\t tg:%f\t tt:%f\t\n", dipi[12],dipi[13],dipi[14],dipi[15]);
		  }else {
             fprintf(stderr, "frequencies:aa:%f\t ac:%f\t ag:%f\t au:%f\t\n",  dipi[0], dipi[1], dipi[2], dipi[3]);
		  	 fprintf(stderr, "           :ca:%f\t cc:%f\t cg:%f\t cu:%f\t\n",  dipi[4], dipi[5], dipi[6], dipi[7]);
             fprintf(stderr, "           :ga:%f\t gc:%f\t gg:%f\t gu:%f\t\n",  dipi[8], dipi[9],dipi[10],dipi[11]);
             fprintf(stderr, "           :ua:%f\t uc:%f\t ug:%f\t uu:%f\t\n",  dipi[12],dipi[13],dipi[14],dipi[15]);
			  
		  }
		} 
		if(verboseDiNucContent==1 || verboseFrequencies==1){
		 fprintf(stderr, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", dipi[0], dipi[1], dipi[2], dipi[3], dipi[4], dipi[5], dipi[6], dipi[7],dipi[8], dipi[9],dipi[10],dipi[11], dipi[12],dipi[13],dipi[14],dipi[15]);
         fprintf(stderr, "sum pi(nk=1) %f\n", sumdipi);
		} 
	 }	
	 if(verboseDiNucContentFinal==1)
	 fprintf(stderr, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", dipi[0], dipi[1], dipi[2], dipi[3], dipi[4], dipi[5], dipi[6], dipi[7],dipi[8], dipi[9],dipi[10],dipi[11], dipi[12],dipi[13],dipi[14],dipi[15]);
	/*********************/
	/*TRI*/
	 sumdipi=0;
	 for(i=0;i<64;i++)
	 sumtripi+= tripi[i];	
	 if(verboseFrequencies==1 && quiet !=1){
       fprintf(stderr, "sum pi(nk=1) %f\n", sumtripi);
	 }  
	 if((sumtripi < 1.00000 - NEARZERO) || (sumtripi > 1.00000 + NEARZERO)){
	    if(verboseFrequencies==1 && quiet !=1 ) fprintf(stderr, "1.New Tripletfreq :\n");
		for(i=0;i<64;i++) tripi[i] = tripi[i]/sumtripi;
		sumtripi=0;	
		for(i=0;i<64;i++)   sumtripi+= tripi[i];	
		if((verboseDiNucContent==1 || verboseFrequencies==1) && quiet !=1 ) fprintf(stderr, "sum pi(nk=1) %f\n", sumtripi);
	
	 }	
	
	
}
/***************************************/
/*  FROM DI TO MONONUCLEOTIDES         */
/***************************************/
/*check dinuc, mononuc, etc to one*/
/*verboseCluster*/
void di2mono(){
 int i,j;
 double sum=0.00;
 double sumpi=0.00;
 

 if(verboseFrequencies==1) fprintf(stderr, "\n1.Mono:");
 if(verboseFrequencies==1) for(i=0;i<4;i++) fprintf(stderr, " %c=%lf ",nucleotide[i], mpi[i]);
 for(i=0;i<4;i++) mpi[i]=0.00;
 if(verboseFrequencies==1){
    fprintf(stderr, "\n2.Mono:");
    for(i=0;i<4;i++)  fprintf(stderr, " %c=%lf ",nucleotide[i], mpi[i]);
	fprintf(stderr, "\n Dinuc:");
	for(i=0;i<16;i++) fprintf(stderr, " %s=%lf",doublet[i], dipi[i]);
    fprintf(stderr, "\n");
 }
  for(i=0;i<4;i++){
      if(verboseFrequencies==1) fprintf(stderr, "\n%c=", nucleotide[i] );
      for(j=0;j<4;j++){
         mpi[i]=mpi[i]+dipi[(i*4)+j];
	     if(verboseFrequencies==1) fprintf(stderr, "%s+", doublet[(i*4)+j]);
	  }
/*	   for(j=0;j<4;j++){
	    mpi[i]=mpi[i]+dipi[(i)+(j*4)];
		if(verboseFrequencies==1) fprintf(stderr, "%s+", doublet[(i)+(j*4)]);
	  }	*/
	 /* mpi[i]=mpi[i]/8;*/
	    mpi[i]=mpi[i]/4;
	  if(verboseFrequencies==1) fprintf(stderr, "=%lf", mpi[i]);
	  sum= sum + mpi[i];
  }	   
     if(verboseFrequencies==1) fprintf(stderr, "\n sum=%lf", sum);
	 sumpi=sum;
	 
	  if((sumpi < 1.00000 - NEARZERO) || (sumpi > 1.00000 + NEARZERO)){
	    for(i=0;i<4;i++)
	      mpi[i] = mpi[i]/sumpi;
		sumpi = 0;
		for(i=0;i<4;i++) sumpi += mpi[i];
		
		if(verboseFrequencies==1){  
		  fprintf(stderr, "\n4.Mono:");
		  fprintf(stderr,"%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
		  fprintf(stderr,"sum %lf \n",sumpi);
		}  
	}	
	
	  if(verboseCluster){
	   fprintf(stderr,"%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]); 
	 }
}

/**********************************************************/
/*  FROM MONO & DI TO TRINUCFREQ (2cluster Approach)      */
/**********************************************************/
  void DiMono2Tri()
 {
   int i;
   double sumtripi;

   for(i=0;i<64;i++) tripi[i]=0.015625;
     if(verboseFrequencies==1) fprintf(stderr, "Compute Trinuc. Freq. with 2 Cluster Approach \n");
     if(verboseFrequencies==1) fprintf(stderr, "1:TriNuc\n");
	  for(i=0;i<64;i++) {
	    tripi[i]=(dipi[i/4]*dipi[i%16])/mpi[(i%16)/4] ;
	    if(verboseFrequencies==1)   fprintf(stderr, "f(%s)=f(%s)*f(%s)/f(%c) ", triplet[i], doublet[i/4] , doublet[i%16],  nucleotide[(i%16)/4]);
	    if(verboseFrequencies==1)	fprintf(stderr, "%f = %f*%f / %f\n", tripi[i], dipi[i/4] , dipi[i%16] , mpi[(i%16)/4]);
	  }
	
		
	sumtripi=0.00;
	for(i=0;i<64;i++) sumtripi += tripi[i];
	if(verboseFrequencies==1) fprintf(stderr, "sum pi(nk=1) %f\n", sumtripi);
	if((sumtripi < 1.00000 - NEARZERO) || (sumtripi > 1.00000 + NEARZERO)){	
	   if(verboseFrequencies==1)fprintf(stderr, "pi(nk=1) unequal 1.00: diff=%e :\n", 1.0 - sumtripi);
	   for(i=0;i<64;i++) tripi[i] = tripi[i]/sumtripi;
	   sumtripi=0.00;
	   for(i=0;i<64;i++) sumtripi += tripi[i];	
	   if(verboseFrequencies==1) fprintf(stderr, "sum pi(nk=1) %f\n", sumtripi);
	}
	if(verboseFrequencies==1){
	  fprintf(stderr, "3:TriNuc\n");
	  for(i=0;i<64;i++)  fprintf(stderr, "f(%s) = %lf\n", triplet[i],tripi[i]);
	  fprintf(stderr, "\n");
	}  	
	
    if(verboseCluster==1) fprintf(stderr, "TriNuc: ");	
	if(verboseCluster==1 || verboseFrequencies==1){
	  for(i=0;i<64;i++) fprintf(stderr, "%.16lf ", tripi[i]);		
	  fprintf(stderr, "\n");	
	}			
   
 }
 




 
 
 
 
 /****************************************************************************************************************/
 /*                CALCULATE RATE MATRIX                                                                         */
 /****************************************************************************************************************/
 /***************************************/
 /*  RATE MATRIX ACCORDING TO REV       */
 /***************************************/
   
    double *rateRV(){

    int i, j, k;
    double mr;
    double *Qij;
	double sumpi;

    Qij=(double *)malloc(16*sizeof(double));
 
    if(verboseQMatrix==1 || verboseFrequencies==1){
      fprintf(stderr,"\n");
	  fprintf(stderr,"\t\t\t***MODEL-(nk=0)***\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n"
                           ,nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
    }						   

    sumpi=0;
	for(i=0;i<4;i++)
	    sumpi += mpi[i];
		
	if(verboseQMatrix==1  || verboseFrequencies==1) fprintf(stderr,"sum of frequencies: %f\n\n", sumpi);	/* if not 1 ....*/
	
	/*if(sumpi!=1.00){ */
	 if((sumpi < 1.00000 - NEARZERO) || (sumpi > 1.00000 + NEARZERO)){
	    for(i=0;i<4;i++)
	      mpi[i] = mpi[i]/sumpi;
		sumpi = 0;
		 for(i=0;i<4;i++) sumpi += mpi[i];
		if(verboseQMatrix==1 || verboseFrequencies==1){  
		  fprintf(stderr,"New\n pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
		  fprintf(stderr,"sum pi(n_k=0)%f\n",sumpi);
		}
     }
	 
    if(verboseQMatrix==1) fprintf(stderr,"rate matrix=gamma1:%f, alpha1:%f, beta1:%f,\n beta2:%f, alpha2:%f, gamma2:%f\n\n",Rmat[0],Rmat[1],Rmat[2],Rmat[3],Rmat[4],Rmat[5]);

    k=0;
	
    /*
	 *Initialize the rate parameter
	 */
	 
    for(i=0;i<3;i++)
        for(j=i+1;j<4;j++)
          if(i*4+j !=11)
             Qij[i*4+j]=Qij[j*4+i]=Rmat[k++];     /*reversible*/
    Qij[3*4+2]=Qij[2*4+3]=1.0;                    /*the last one with 1.0*/
	
	
	if(transitiontransverion==1 && prmat==0) {
	     K=2*trtv;
       /*		 fprintf(stderr, "\n HKY Model or Kimura Model with k=%lf trtv=%lf\n", K, trtv);*/
	     Qij[0*4+2]=Qij[0*4+2]*K;
		 Qij[1*4+3]=Qij[1*4+3]*K;
		 Qij[2*4+0]=Qij[2*4+0]*K;
		 Qij[3*4+1]=Qij[3*4+1]*K;
	}
	if(transitiontransverion==1 && prmat==1){
	   fprintf(stderr, "You cannot use TS/TV Parameter and the GTR model!!\n");
	   exit (0);
	}   

	
	 /*
	 * print ratenmatrix(REV)
	 */
	if(verboseQMatrix==1){
      fprintf(stderr,"unnormalise Rk=0 (independent instantaneous substitution matrix)\n:\n");
     
      for(i=0;i<4;i++){
         fprintf(stderr,"%c ||", nucleotide[i]);
         for(j=0;j<4;j++)
            fprintf(stderr, " %8.5f |",Qij[i*4+j]);
            fprintf(stderr, "\n");
	  }
	
	}
	
	
	/*
	 * multiplied by stationary distribution
	 */
    for(i=0;i<4;i++)
       for(j=0;j<4;j++)
         Qij[i*4+j]*=mpi[j];

	/*
	 * the diagonal values
	 */
    mr=0;
    for (i=0; i<4; i++){
        Qij[i*4+i]=0;                                              /*diagonal entries zero*/
        Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);   /*mathematical requiremnet: sum of the row*/
        mr-=mpi[i]*Qij[i*4+i];                                     /*sum of the diagonal values*pi */
    }

    /*
	 * print ratenmatrix(REV)
	 */
	if(verboseQMatrix==1){
      fprintf(stderr,"unnormalise Qk=0 (independent instantaneous substitution matrix)\n:\n");
     
      for(i=0;i<4;i++){
         fprintf(stderr,"%c ||", nucleotide[i]);
         for(j=0;j<4;j++)
            fprintf(stderr, " %8.5f |",Qij[i*4+j]);
            fprintf(stderr, "\n");
	  }
	
         fprintf(stderr, "sum(diagonal.values*pi):expected number of substitutions = %f\n", mr);
         fprintf(stderr, "\n");
	  	 
    }
   
    /*
	 * normalize the rate matrix
	 */
      for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        Qij[i*4+j] /= mr;        /*mr:sum of the diagonal values*pi*/


    if(verboseQMatrix==1){
      fprintf(stderr,"  normalise Qk=0 (independent instantaneous substitution matrix)\n");
      for(i=0;i<4;i++){
            fprintf(stderr,"%c ||", nucleotide[i]);
			for(j=0;j<4;j++)
               fprintf(stderr," %8.5f |",Qij[i*4+j]);
            fprintf(stderr,"\n");
	   }
     }			

    /*
	 * diagonal values
	 */
    mr=0;
    for (i=0; i<4; i++){
       mr-=mpi[i]*Qij[i*4+i];      /*sum of the diagonal values: must be one!*/

    }
	
	if(verboseRate==1 || verboseQMatrix==1)
     fprintf(stderr, " Q(nk=0): sum(diagonal.values*pi):expected number of substitutions = %f\n", mr);

    return(Qij);
}





    double *rateRV_old(){

    int i, j, k;
    double mr;
    double *Qij;
	double sumpi;
	
	/*fprintf(stderr, "check model.c %lf", trtv);*/

    Qij=(double *)malloc(16*sizeof(double));
 
    if(verboseQMatrix==1 || verboseFrequencies==1 ||(verboseQMatrix==1 &&printfile==1) ||(verboseFrequencies==1 && printfile==1)){
      if(printfile==1 && verboseFrequencies==1){
	   fprintf(freq_fv,"\n");
	   fprintf(freq_fv,"\t\t\t***MODEL-(nk=0)***\n");
       fprintf(freq_fv,"\n");
       fprintf(freq_fv,"pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n"
                           ,nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
	 }else 
	   if(printfile==1 && verboseQMatrix==1){
	   fprintf(matr_fv,"\n");
	   fprintf(matr_fv,"\t\t\t***MODEL-(nk=0)***\n");
       fprintf(matr_fv,"\n");
       fprintf(matr_fv,"pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n"
                           ,nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);	
      }else				
	  fprintf(stderr,"\n");
	  fprintf(stderr,"\t\t\t***MODEL-(nk=0)***\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n"
                           ,nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
    }						   

    sumpi=0;
	for(i=0;i<4;i++)
	    sumpi += mpi[i];
		
	if((verboseQMatrix==1  || verboseFrequencies==1) && printfile != 1) fprintf(stderr,"sum of frequencies: %f\n\n", sumpi);	/* if not 1 ....*/

	/*if(sumpi!=1.00){ */
	 if((sumpi < 1.00000 - NEARZERO) || (sumpi > 1.00000 + NEARZERO)){
	    for(i=0;i<4;i++)
	      mpi[i] = mpi[i]/sumpi;
		sumpi = 0;
		 for(i=0;i<4;i++) sumpi += mpi[i];
		if(verboseQMatrix==1 || verboseFrequencies==1 ||(verboseQMatrix==1 &&printfile==1) ||(verboseFrequencies==1 && printfile==1)){  
		 if(printfile==1 && verboseQMatrix==1){
		  fprintf(matr_fv,"New\n pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
		  fprintf(matr_fv,"sum pi(n_k=0)%f\n",sumpi);
		 }else if(printfile==1 && verboseFrequencies==1){
		  fprintf(freq_fv,"New\n pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
		  fprintf(freq_fv,"sum pi(n_k=0)%f\n",sumpi);
		 }else
		  fprintf(stderr,"New\n pi(nk=0):%c:%f, %c:%f, %c:%f, %c:%f\n",nucleotide[0],mpi[0],nucleotide[1],mpi[1],nucleotide[2],mpi[2],nucleotide[3],mpi[3]);
		  fprintf(stderr,"sum pi(n_k=0)%f\n",sumpi);
		}
     }
	 
    if(verboseQMatrix==1 || (verboseQMatrix==1 && printfile==1)){
	  if(printfile==1) 
	    fprintf(stderr,"rate matrix=gamma1:%f, alpha1:%f, beta1:%f,\n beta2:%f, alpha2:%f, gamma2:%f\n\n",Rmat[0],Rmat[1],Rmat[2],Rmat[3],Rmat[4],Rmat[5]);
	  else
	   fprintf(stderr,"rate matrix=gamma1:%f, alpha1:%f, beta1:%f,\n beta2:%f, alpha2:%f, gamma2:%f\n\n",Rmat[0],Rmat[1],Rmat[2],Rmat[3],Rmat[4],Rmat[5]);
    }	  
	  
      k=0;
	
    /*
	 *Initialize the rate parameter
	 */
	 
    for(i=0;i<3;i++)
        for(j=i+1;j<4;j++)
          if(i*4+j !=11)
             Qij[i*4+j]=Qij[j*4+i]=Rmat[k++];     /*reversible*/
    Qij[3*4+2]=Qij[2*4+3]=1.0;                    /*the last one with 1.0*/
	
	
	/*
	 * multiplied by stationary distribution
	 */
    for(i=0;i<4;i++)
       for(j=0;j<4;j++)
         Qij[i*4+j]*=mpi[j];
		 
		 
		 
   if(transitiontransverion==1) {
	    Qij[0*4+1]=Qij[0*4+1]*trtv;
		Qij[1*4+0]=Qij[1*4+0]*trtv;
		Qij[0*4+3]=Qij[0*4+3]*trtv;
		Qij[3*4+0]=Qij[3*4+0]*trtv;
		Qij[1*4+2]=Qij[1*4+2]*trtv;
		Qij[2*4+1]=Qij[2*4+1]*trtv;
		Qij[2*4+3]=Qij[2*4+3]*trtv;
		Qij[3*4+2]=Qij[3*4+2]*trtv;
     }
		
		
	  	 

	/*
	 * the diagonal values
	 */
    mr=0;
    for (i=0; i<4; i++){
        Qij[i*4+i]=0;                                              /*diagonal entries zero*/
        Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);   /*mathematical requiremnet: sum of the row*/
        mr-=mpi[i]*Qij[i*4+i];                                     /*sum of the diagonal values*pi */
    }

    /*
	 * print ratenmatrix(REV)
	 */
	if(verboseQMatrix==1 ){
	/*  if(printfile==1){
      fprintf(matr_fv,"unnormalise Qk=0 (independent instantaneous substitution matrix)\n:\n");
      for(i=0;i<4;i++){
         fprintf(matr_fv,"%c ||", nucleotide[i]);
         for(j=0;j<4;j++)
            fprintf(matr_fv, " %8.5f |",Qij[i*4+j]);
            fprintf(matr_fv, "\n");
	  } else{*/
	  fprintf(stderr,"unnormalise Qk=0 (independent instantaneous substitution matrix)\n:\n");
      for(i=0;i<4;i++){
         fprintf(stderr,"%c ||", nucleotide[i]);
         for(j=0;j<4;j++)
            fprintf(stderr, " %8.5f |",Qij[i*4+j]);
            fprintf(stderr, "\n");
	 /* }*/		
	  }
	
         fprintf(stderr, "sum(diagonal.values*pi):expected number of substitutions = %f\n", mr);
         fprintf(stderr, "\n");
	  	 
    }
   
    /*
	 * normalize the rate matrix
	 */
      for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        Qij[i*4+j] /= mr;        /*mr:sum of the diagonal values*pi*/


    if(verboseQMatrix==1){
      fprintf(stderr,"  normalise Qk=0 (independent instantaneous substitution matrix)\n");
      for(i=0;i<4;i++){
            fprintf(stderr,"%c ||", nucleotide[i]);
			for(j=0;j<4;j++)
               fprintf(stderr," %8.5f |",Qij[i*4+j]);
            fprintf(stderr,"\n");
	   }
     }			

    /*
	 * diagonal values
	 */
    mr=0;
    for (i=0; i<4; i++){
       mr-=mpi[i]*Qij[i*4+i];      /*sum of the diagonal values: must be one!*/

    }
	
	if(verboseRate==1 || verboseQMatrix==1)
     fprintf(stderr, " Q(nk=0): sum(diagonal.values*pi):expected number of substitutions = %f\n", mr);

    return(Qij);
}





  /******************************************************************/
  /*  TRIPLET RATE MATRIX FOR THE DINUCLEOTIDCONTENT                */
  /******************************************************************/

/*tanja:trtv*/
double **rateT2(){
    int i,j, ii,jj ,iii,jjj;
    double dig, rowsum[64];
    double **QT2ij, **RT2ij;
	/*double sumtripi;*/
	
	
	/********************************************
	 *Speicher allozieren
	 ********************************************/
    QT2ij=(double**)malloc(64*sizeof(double*));
    for(i=0;i<64;i++)
       QT2ij[i]=(double *)malloc(64*sizeof(double));
	   
	  for(i=0;i<64;i++)
	   for(j=0;j<64;j++) QT2ij[i][j]=0.00;   
	   
	RT2ij=(double**)malloc(64*sizeof(double*));
    for(i=0;i<64;i++)
       RT2ij[i]=(double *)malloc(64*sizeof(double));

    for(i=0;i<64;i++)
	   for(j=0;j<64;j++) RT2ij[i][j]=1.00;
   

    dig=0.00;
	
	/********************************************
	 * Rmatrix
	 ********************************************/
		if(prmat==1 && prmat2t==0){			  
	  for(i=0;i<64;i++)
      {
        ii=(i+4)/4;
		iii=(ii+3)/4;
        for(j=0;j<64;j++)
        {
        jj=(j+4)/4;
		jjj=(jj+3)/4;
		 if(((i)%16)/4==0 && ((j)%16)/4==1) RT2ij[i][j]=Rmat[0];
		 if(((i)%16)/4==0 && ((j)%16)/4==2) RT2ij[i][j]=Rmat[1];
	 	 if(((i)%16)/4==0 && ((j)%16)/4==3) RT2ij[i][j]=Rmat[2];
		 if(((i)%16)/4==1 && ((j)%16)/4==2) RT2ij[i][j]=Rmat[3];
		 if(((i)%16)/4==1 && ((j)%16)/4==3) RT2ij[i][j]=Rmat[4];
		 if(((i)%16)/4==2 && ((j)%16)/4==3) RT2ij[i][j]=Rmat[5];
		 if(((i)%16)/4==1 && ((j)%16)/4==0) RT2ij[i][j]=Rmat[0];
		 if(((i)%16)/4==2 && ((j)%16)/4==0) RT2ij[i][j]=Rmat[1];
	 	 if(((i)%16)/4==3 && ((j)%16)/4==0) RT2ij[i][j]=Rmat[2];
		 if(((i)%16)/4==2 && ((j)%16)/4==1) RT2ij[i][j]=Rmat[3];
		 if(((i)%16)/4==3 && ((j)%16)/4==1) RT2ij[i][j]=Rmat[4];
		 if(((i)%16)/4==3 && ((j)%16)/4==2) RT2ij[i][j]=Rmat[5];
		 
		  if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           RT2ij[i][j]=0;
          if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           RT2ij[i][j]=0;
          if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) RT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
		  if(i==j) RT2ij[i][j]=0.00;                                        /*diagonale*/
		 }
	   }
	   
}  

/*if(prmat2t==1){			  
   for(j=0;j<96;j++) fprintf(stderr, " rmat%d Rmat2T %lf\n", j,Rmat2T[j]);	 
	for(k=0;k<16;k++){
	    for(i=0;i<64;i++)
        {
        ii=(i+4)/4;
		iii=(ii+3)/4;
	    for(j=0;j<64;j++)
        {
        jj=(j+4)/4;
		jjj=(jj+3)/4;
		if( (k==(j%4)+(jjj-1)*4) && (k== i%4 + (iii-1)*4)){

		 if(((i)%16)/4==0 && ((j)%16)/4==1) RT2ij[i][j]=Rmat2T[0+(k*6)];
		 if(((i)%16)/4==0 && ((j)%16)/4==2) RT2ij[i][j]=Rmat2T[1+(k*6)];
	 	 if(((i)%16)/4==0 && ((j)%16)/4==3) RT2ij[i][j]=Rmat2T[2+(k*6)];
		 if(((i)%16)/4==1 && ((j)%16)/4==2) RT2ij[i][j]=Rmat2T[3+(k*6)];
		 if(((i)%16)/4==1 && ((j)%16)/4==3) RT2ij[i][j]=Rmat2T[4+(k*6)];
		 if(((i)%16)/4==2 && ((j)%16)/4==3) RT2ij[i][j]=Rmat2T[5+(k*6)];
		 if(((i)%16)/4==1 && ((j)%16)/4==0) RT2ij[i][j]=Rmat2T[0+(k*6)];
		 if(((i)%16)/4==2 && ((j)%16)/4==0) RT2ij[i][j]=Rmat2T[1+(k*6)];
	 	 if(((i)%16)/4==3 && ((j)%16)/4==0) RT2ij[i][j]=Rmat2T[2+(k*6)];
		 if(((i)%16)/4==2 && ((j)%16)/4==1) RT2ij[i][j]=Rmat2T[3+(k*6)];
		 if(((i)%16)/4==3 && ((j)%16)/4==1) RT2ij[i][j]=Rmat2T[4+(k*6)];
		 if(((i)%16)/4==3 && ((j)%16)/4==2) RT2ij[i][j]=Rmat2T[5+(k*6)];
         }
		  if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )              
		             RT2ij[i][j]=0;
          if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     
           RT2ij[i][j]=0;
          if(matrixSH==0){
		      if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) RT2ij[i][j]=0;}            
			  if(i==j) RT2ij[i][j]=0.00;                                       
			}	
		}
	   } 
	 }   */
	 
			if(transitiontransverion==1) {
	     K=2*trtv;
       /*		 fprintf(stderr, "\nHKY Model or Kimura Model with k=%lf trtv=%lf\n", K, trtv);*/
		 for(i=0;i<64;i++){ /* wieso vorher 16?? 0312*/
		   ii=(i+4)/4;
		   iii=(ii+3)/4;
		       for(j=0;j<64;j++){
		       jj=(j+4)/4;
		       jjj=(jj+3)/4;
	     if(((i)%16)/4==0 && ((j)%16)/4==2) RT2ij[i][j]=K*RT2ij[i][j];
		 if(((i)%16)/4==1 && ((j)%16)/4==3) RT2ij[i][j]=K*RT2ij[i][j];
		 if(((i)%16)/4==2 && ((j)%16)/4==0) RT2ij[i][j]=K*RT2ij[i][j];
		 if(((i)%16)/4==3 && ((j)%16)/4==1) RT2ij[i][j]=K*RT2ij[i][j];
   
		 if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           RT2ij[i][j]=0;
		 if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           RT2ij[i][j]=0;
          if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) RT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
		  if(i==j) RT2ij[i][j]=0.00;                 
	}}}
	if(transitiontransverion==1 && prmat==1){
	   fprintf(stderr, "You cannot use TS/TV Parameter and the GTR model!!\n");
	   exit (0);
	}  
	
	
/*	  if(CpG==1){
		   
			 RT2ij[6][14]=CpGCpA;
			 RT2ij[22][30]=CpGCpA;
			 RT2ij[38][46]=CpGCpA;
			 RT2ij[54][62]=CpGCpA;
			
			 RT2ij[24][16]=CpGTpG;
			 RT2ij[25][17]=CpGTpG;
			 RT2ij[26][18]=CpGTpG;
			 RT2ij[27][19]=CpGTpG;
		
			  fprintf(stderr, "\n%s %s %lf\n",triplet[6],  triplet[14],CpGCpA);
			  fprintf(stderr, "%s %s %lf\n",  triplet[22], triplet[30],CpGCpA);
			  fprintf(stderr, "%s %s %lf\n",  triplet[38], triplet[46],CpGCpA);
			  fprintf(stderr, "%s %s %lf\n",  triplet[54], triplet[62],CpGCpA);
			  fprintf(stderr, "%s %s %lf\n",  triplet[24], triplet[16],CpGTpG);
			  fprintf(stderr, "%s %s %lf\n",  triplet[25], triplet[17],CpGTpG);
			  fprintf(stderr, "%s %s %lf\n",  triplet[26], triplet[18],CpGTpG);
			  fprintf(stderr, "%s %s %lf\n",  triplet[27], triplet[19],CpGTpG);
   }*/
   
   
     if(verboseQMatrix==1){
     fprintf(stderr,"\n RT2 \n");
     fprintf(stderr,"Triplet2:dig=%f\n",dig);
     fprintf(stderr,"\n");
     fprintf(stderr,"        ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");
	  fprintf(stderr,"        ");
	 for(i=0; i<64; i++) fprintf(stderr,"%6c|",nucleotide[((i)%16)/4]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
	   fprintf(stderr,"%c||",nucleotide[((i)%16)/4]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",RT2ij[i][j]);
       fprintf(stderr,"\n");
     }
	 fprintf(stderr,"Triplet2: dig=%f\n",dig);
  }
	 
	
	/********************************************
	 *old
	 ********************************************/


	for(i=0;i<64;i++)
    {
      ii=(i+4)/4;
      iii=(ii+3)/4;
      for(j=0;j<64;j++)
      {
        jj=(j+4)/4;
        jjj=(jj+3)/4;
        QT2ij[i][j]=tripi[j]*RT2ij[i][j]; /*tanja:trtv*/
        if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) QT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
        if(i==j) QT2ij[i][j]=0.00;                                        /*diagonale*/
      }
    }

   dig=0.00;
   for(i=0; i<64; i++){
      rowsum[i]=0.00;
      for(j=0; j<64; j++){
         if(i!=j)  rowsum[i]=rowsum[i]+ QT2ij[i][j];}
	  QT2ij[i][i]= -1*rowsum[i];                         /*Diagonalentries*/ 
      dig=dig + tripi[i]*rowsum[i];                        /*Diagonalentries multipli. with the rowsum together*/
    }                                                   

   if(verboseQMatrix==1){
     fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"Triplet2:dig=%f\n",dig);
     fprintf(stderr,"\n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");
	 
	 for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
	 fprintf(stderr,"Triplet2: dig=%f\n",dig);
  }
    
  for(i=0;i<64;i++)
        for(j=0;j<64;j++)
                QT2ij[i][j] /= dig ;                     /*AB HIER PI*RATE -> NORMIERUNG*/

  if(verboseQMatrix==1){
     fprintf(stderr,"\n");
	 fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
   }

    dig=0;
    for (i=0; i<64; i++)
       dig -= tripi[i]*QT2ij[i][i];      /*sum of the diagonal values*/

    if(verboseRate==1 || verboseQMatrix==1)
        fprintf(stderr, "Q2(nk=2): sum(diagonal.values*pi):expected number of substitutions = %f\n", dig);
		
	free(RT2ij);
	return(QT2ij);
}
   
double **rateT2_old(){
    int i,j, ii,jj ,iii,jjj;
    double dig, rowsum[64];
    double **QT2ij;
	/*double sumtripi;*/

    QT2ij=(double**)malloc(64*sizeof(double*));
    for(i=0;i<64;i++)
       QT2ij[i]=(double *)malloc(64*sizeof(double));

	for(i=0;i<64;i++)
    {
      ii=(i+4)/4;
      iii=(ii+3)/4;
      for(j=0;j<64;j++)
      {
        jj=(j+4)/4;
        jjj=(jj+3)/4;
        QT2ij[i][j]=tripi[j];
        if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) QT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
        if(i==j) QT2ij[i][j]=0.00;                                        /*diagonale*/
      }
    }

   dig=0.00;
   for(i=0; i<64; i++){
      rowsum[i]=0.00;
      for(j=0; j<64; j++){
         if(i!=j)  rowsum[i]=rowsum[i]+ QT2ij[i][j];}
	  QT2ij[i][i]= -1*rowsum[i];                         /*Diagonalentries*/ 
      dig=dig + tripi[i]*rowsum[i];                        /*Diagonalentries multipli. with the rowsum together*/
    }                                                   

   if(verboseQMatrix==1){
     fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"Triplet2:dig=%f\n",dig);
     fprintf(stderr,"\n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");
	 
	 for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
	 fprintf(stderr,"Triplet2: dig=%f\n",dig);
  }
    
  for(i=0;i<64;i++)
        for(j=0;j<64;j++)
                QT2ij[i][j] /= dig ;                     /*AB HIER PI*RATE -> NORMIERUNG*/

  if(verboseQMatrix==1){
     fprintf(stderr,"\n");
	 fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
   }

    dig=0;
    for (i=0; i<64; i++)
       dig -= tripi[i]*QT2ij[i][i];      /*sum of the diagonal values*/

    if(verboseRate==1 || verboseQMatrix==1)
        fprintf(stderr, "Q2(nk=2): sum(diagonal.values*pi):expected number of substitutions = %f\n", dig);
	return(QT2ij);
}
   



  double **rateT2sum(){
    int i,j, ii,jj ,iii,jjj;
    double dig, rowsum[64];
    double **QT2ij;
	int dinuc_13, nuc_2;
	int first, second, third, l;
	double sum123;

    QT2ij=(double**)malloc(64*sizeof(double*));
    for(i=0;i<64;i++)
       QT2ij[i]=(double *)malloc(64*sizeof(double));

	for(i=0;i<64;i++)
    {
      ii=(i+4)/4;
      iii=(ii+3)/4;
	  dinuc_13= (i%4)+(4*(iii-1))  ;
	  nuc_2= (ii-1)%4;
	 /* printf("%s->%d,%c\n",triplet[i],nuc_2,nucleotide[nuc_2]);*/
      for(j=0;j<64;j++)
      {
        jj=(j+4)/4;
        jjj=(jj+3)/4;
		if(i!=j){
		   for(l=0;l<4;l++){
		   first= dinuc_13+(l*16);
		   third= (dinuc_13*4) +l;
		   second=dinuc_13%4+(l*4) + ((dinuc_13/4)*16);
		   sum123= (tripi[first]+tripi[second]+tripi[third])/3;
		  /* printf("%s-%s: %s %s %s \n",triplet[j], doublet[dinuc_13],triplet[first],triplet[second],triplet[third]);*/
           QT2ij[i][j]=sum123;
		   }
		}
        if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) QT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
        if(i==j) QT2ij[i][j]=0.00;                                        /*diagonale*/
      }
    }

   dig=0.00;
   for(i=0; i<64; i++){
      rowsum[i]=0.00;
      for(j=0; j<64; j++){
         if(i!=j)  rowsum[i]=rowsum[i]+ QT2ij[i][j];}
	  QT2ij[i][i]= -1*rowsum[i];                         /*Diagonalentries*/ 
      dig=dig + tripi[i]*rowsum[i];                        /*Diagonalentries multipli. with the rowsum together*/
    }                                                   

   if(verboseQMatrix==1){
     fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"Triplet2:dig=%f\n",dig);
     fprintf(stderr,"\n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");
	 
	 for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
	 fprintf(stderr,"Triplet2: dig=%f\n",dig);
  }
    
  for(i=0;i<64;i++)
        for(j=0;j<64;j++)
                QT2ij[i][j] /= dig ;                     /*AB HIER PI*RATE -> NORMIERUNG*/

  if(verboseQMatrix==1){
     fprintf(stderr,"\n");
	 fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
   }

    dig=0;
    for (i=0; i<64; i++)
       dig -= tripi[i]*QT2ij[i][i];      /*sum of the diagonal values*/

    if(verboseRate==1 || verboseQMatrix==1)
        fprintf(stderr, "Q2(nk=2): sum(diagonal.values*pi):expected number of substitutions = %f\n", dig);
	return(QT2ij);
}
 

   
double **rateT2CpG(){
    int i,j, ii,jj ,iii,jjj;
    double dig, rowsum[64];
    double **QT2ij;
	/*double sumtripi;*/
	int dinuc_j_left, dinuc_j_right;
	int dinuc_i_left, dinuc_i_right;
	
	fprintf(stderr,"\n r(CpG)=%d\n", cpg);

    QT2ij=(double**)malloc(64*sizeof(double*));
    for(i=0;i<64;i++)
       QT2ij[i]=(double *)malloc(64*sizeof(double));

	for(i=0;i<64;i++)
    {
      ii=(i+4)/4;
      iii=(ii+3)/4;
	  dinuc_i_left=ii-1;
	  dinuc_i_right= (i%16);
      for(j=0;j<64;j++)
      {
        jj=(j+4)/4;
        jjj=(jj+3)/4;
		dinuc_j_left=jj-1;
	    dinuc_j_right=(j%16);
        QT2ij[i][j]=tripi[j];
		if(dinuc_i_left  == 6 && dinuc_j_left==4)  {fprintf(stderr,"%s->%s plus %d\n", triplet[i],triplet[j],cpg); QT2ij[i][j]=QT2ij[i][j]+cpg;}
		if(dinuc_i_left  == 6 && dinuc_j_left==14) {fprintf(stderr,"%s->%s plus %d\n", triplet[i],triplet[j],cpg); QT2ij[i][j]=QT2ij[i][j]+cpg;}
		if(dinuc_i_right == 6 && dinuc_j_right==4) {fprintf(stderr,"%s->%s plus %d\n", triplet[i],triplet[j],cpg); QT2ij[i][j]=QT2ij[i][j]+cpg;}
		if(dinuc_i_right == 6 && dinuc_j_right==4) {fprintf(stderr,"%s->%s plus %d\n", triplet[i],triplet[j],cpg);QT2ij[i][j]=QT2ij[i][j]+cpg;}
		if(dinuc_i_right == 6 && dinuc_j_right==14){fprintf(stderr,"%s->%s plus %d\n", triplet[i],triplet[j],cpg); QT2ij[i][j]=QT2ij[i][j]+cpg;}
		
        if( (ii != jj) &&  (i-4*(ii-1) != j-4*(jj-1)) )               /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if( (iii != jjj) && ( (ii-4*(iii-1))!=(jj -4*(jjj-1)) ) )     /*Hammingdistance==2*/
           QT2ij[i][j]=0;
        if(matrixSH==0){
		     if( ((ii-1)-4*(iii-1) == (jj-1)-4*(jjj-1))  ) QT2ij[i][j]=0;}               /*expand model, only substitution on one position*/
        if(i==j) QT2ij[i][j]=0.00;                                        /*diagonale*/
		
      }
    }

   dig=0.00;
   for(i=0; i<64; i++){
      rowsum[i]=0.00;
      for(j=0; j<64; j++){
         if(i!=j)  rowsum[i]=rowsum[i]+ QT2ij[i][j];}
	  QT2ij[i][i]= -1*rowsum[i];                         /*Diagonalentries*/ 
      dig=dig + tripi[i]*rowsum[i];                        /*Diagonalentries multipli. with the rowsum together*/
    }                                                   

   if(verboseQMatrix==1){
     fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"Triplet2:dig=%f\n",dig);
     fprintf(stderr,"\n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");
	 
	 for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
	 fprintf(stderr,"Triplet2: dig=%f\n",dig);
  }
    
  for(i=0;i<64;i++)
        for(j=0;j<64;j++)
                QT2ij[i][j] /= dig ;                     /*AB HIER PI*RATE -> NORMIERUNG*/

  if(verboseQMatrix==1){
     fprintf(stderr,"\n");
	 fprintf(stderr,"\n QT2 \n");
     fprintf(stderr,"     ");
     for(i=0; i<64; i++) fprintf(stderr,"%6s|",triplet[i]);
     fprintf(stderr,"\n");

     for(i=0; i<64; i++){
       fprintf(stderr,"%s||",triplet[i]);
       for(j=0; j<64; j++)
       fprintf(stderr,"% .3f|",QT2ij[i][j]);
       fprintf(stderr,"\n");
     }
   }

    dig=0;
    for (i=0; i<64; i++)
       dig -= tripi[i]*QT2ij[i][i];      /*sum of the diagonal values*/

    if(verboseRate==1 || verboseQMatrix==1)
        fprintf(stderr, "Q2(nk=2): sum(diagonal.values*pi):expected number of substitutions = %f\n", dig);
	return(QT2ij);
}   




