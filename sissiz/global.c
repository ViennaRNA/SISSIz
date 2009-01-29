#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "global.h"



int verbose=0, verboseMemory=0, quiet=0, relaxedPhylip, randomgeneratordkiss, FastaAln, Clustal;
int verboseQMatrix, verboseWarning, verboseStep, verboseNeighborMatrix, verboseNeighborCt, verboseDiNucContent,verboseDiNucContentFinal, verboseMarkov, verboseEvolutiontime, verboseRate,verboseFrequencies, verboserandom, verbosewaiting, verbosewaitingcount,userSeed, verboseRandomSeed;
unsigned long randomSeed=-1; 
long totalMem=0;
int *randstream; /* current SPRNG random numberstream */
int seqlen;	
double minimumvalue;	  

double** Pij;


  /*************************************/
  /*change 'string' to 'int'*/
    int getintNumber (char *sNumber) {
	 char *endPtr_;
	 int number_ = strtol (sNumber, &endPtr_, 10);
	 return number_;
   }

  /*********************************************************************************************************/
  /*change 'string" to 'double'*/
   double getdoubleNumber(char *sNumber){
            char *endPtr_;
            double number_;
            number_=strtod(sNumber, &endPtr_);
            return number_;
    }
  /*********************************************************************************************************/
   int GetDoubleParams(int argc, char **argv, int *argn, char *pos, int numParams, double *params)
   {
	   int i;
	   char *st, buf[256];

	   i=0;
	   strcpy(buf, pos);
	   st=strtok(buf, "\t,/");
	   do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}

		if (sscanf(st, "%lf", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	   } while (i<numParams);

      return 0;
   }

/*************************************/
int GetIntParams(int argc, char **argv, int *argn, char *pos, int numParams, int *params)
{
	int i;
	char *st, buf[256];

	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}

		if (sscanf(st, "%d", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	} while (i<numParams);

	return 0;
}


/*************************************/
int GetUnsignedLongParams(int argc, char **argv, int *argn, char *pos, int numParams, unsigned long *params)
{
	int i;
	char *st, buf[256];
	
	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}
		
		if (sscanf(st, "%lu", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	} while (i<numParams);

	return 0;
}


/*************************************/
/*seq-gen*/
void *AllocMem(long n, char *name, char *func, int showInfo)
{
	void *P;
	
	if ( (P=malloc(n))==NULL ) {
		fprintf(stderr, "Out of memory allocating '%s': %s()\n", name, func);
		exit(0);
	}
	
	totalMem+=n;

	if (showInfo && verboseMemory)
		fprintf(stderr, "%s in %s() - %ld bytes\n", name, func, n);
	
	return P;
}

/*************************************/
/*seq-gen*/
void *CAllocMem(long n, char *name, char *func, int showInfo)
{
	void *P;
	
	if ( (P=calloc(n, 1))==NULL ) {
		fprintf(stderr, "Out of memory allocating '%s': %s()\n", name, func);
		exit(0);
	}
	
	totalMem+=n;
	if (showInfo && verboseMemory)
		fprintf(stderr, "%s in %s() - %ld bytes\n", name, func, n);
	
	return P;
}






 /***************translate nucleotide in int****************************************************/
   int nucleoTO123(char letter){
   int  zahl=-1;

   switch(letter){

    case 'A': case 'a':                     zahl = 0; break;
    case 'C': case 'c':                     zahl = 1; break;
    case 'G': case 'g':                     zahl = 2; break;
    case 'T': case 't':                     zahl = 3; break;
    case 'U': case 'u':                     zahl = 3; break;
    case '-':                               zahl = 4; break;
    default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter);break;
   }

   return(zahl);
   }
  /***************translate doublets in int****************************************************/
   int  doubletsTO15(char letter1, char letter2){
   int zahl=-1;


      switch(letter1){

      case 'A': case 'a':
         switch(letter2){
                case 'A': case 'a':                     zahl =  0; break;
                case 'C': case 'c':                     zahl =  1; break;
                case 'G': case 'g':                     zahl =  2; break;
                case 'T': case 't':
                case 'U': case 'u':                     zahl =  3; break;
                default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter2);break;
      }
      break;
      case 'C': case 'c':
          switch(letter2){
                case 'A': case 'a':                     zahl =  4; break;
                case 'C': case 'c':                     zahl =  5; break;
                case 'G': case 'g':                     zahl =  6; break;
                case 'T': case 't':
                case 'U': case 'u':                     zahl =  7; break;
				default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter2);break;
      }
      break;
      case 'G': case 'g':
            switch(letter2){
                case 'A': case 'a':                     zahl =  8; break;
                case 'C': case 'c':                     zahl =  9; break;
                case 'G': case 'g':                     zahl = 10; break;
                case 'T': case 't':
                case 'U': case 'u':                     zahl = 11; break;
                default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter2);break;
      }
      break;
      case 'T': case 't':
      case 'U': case 'u':
             switch(letter2){
                case 'A': case 'a':                     zahl = 12; break;
                case 'C': case 'c':                     zahl = 13; break;
                case 'G': case 'g':                     zahl = 14; break;
                case 'T': case 't':
                case 'U': case 'u':                     zahl = 15; break;
                default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter2);break;
      }
    }
    return(zahl);
    }

     /***************translate doublets in int****************************************************/
   int  tripletTO64(char letter1, char letter2, char letter3){
   int zahl=-1;


      switch(letter1){
	  case 'A': case 'a':
         switch(letter2){
                case 'A': case 'a':                     
				         switch(letter3){
                         case 'A': case 'a':                     zahl =  0; break;
                         case 'C': case 'c':                     zahl =  1; break;
                         case 'G': case 'g':                     zahl =  2; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl =  3; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'C': case 'c':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl =  4; break;
                         case 'C': case 'c':                     zahl =  5; break;
                         case 'G': case 'g':                     zahl =  6; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl =  7; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'G': case 'g':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl =  8; break;
                         case 'C': case 'c':                     zahl =  9; break;
                         case 'G': case 'g':                     zahl = 10; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 11; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'U': case 'u':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 12; break;
                         case 'C': case 'c':                     zahl = 13; break;
                         case 'G': case 'g':                     zahl = 14; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 15; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
               break;
		   }
	  break;
		 case 'C': case 'c':
         switch(letter2){
                case 'A': case 'a':                     
				         switch(letter3){
                         case 'A': case 'a':                     zahl = 16; break;
                         case 'C': case 'c':                     zahl = 17; break;
                         case 'G': case 'g':                     zahl = 18; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 19; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'C': case 'c':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 20; break;
                         case 'C': case 'c':                     zahl = 21; break;
                         case 'G': case 'g':                     zahl = 22; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 23; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'G': case 'g':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 24; break;
                         case 'C': case 'c':                     zahl = 25; break;
                         case 'G': case 'g':                     zahl = 26; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 27; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'U': case 'u':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 28; break;
                         case 'C': case 'c':                     zahl = 29; break;
                         case 'G': case 'g':                     zahl = 30; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 31; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
               break;
		   }
		break;
		case 'G': case 'g':
         switch(letter2){
                case 'A': case 'a':                     
				         switch(letter3){
                         case 'A': case 'a':                     zahl = 32; break;
                         case 'C': case 'c':                     zahl = 33; break;
                         case 'G': case 'g':                     zahl = 34; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 35; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'C': case 'c':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 36; break;
                         case 'C': case 'c':                     zahl = 37; break;
                         case 'G': case 'g':                     zahl = 38; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 39; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'G': case 'g':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 40; break;
                         case 'C': case 'c':                     zahl = 41; break;
                         case 'G': case 'g':                     zahl = 41; break;
                         case 'T': case 't':                     
                         case 'U': case 'u':                     zahl = 43; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'U': case 'u':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 44; break;
                         case 'C': case 'c':                     zahl = 45; break;
                         case 'G': case 'g':                     zahl = 46; break;
                         case 'T': case 't':                     
                         case 'U': case 'u':                     zahl = 47; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
               break;
		   }
		break;
		case 'U': case 'u':
         switch(letter2){
                case 'A': case 'a':                     
				         switch(letter3){
                         case 'A': case 'a':                     zahl = 48; break;
                         case 'C': case 'c':                     zahl = 49; break;
                         case 'G': case 'g':                     zahl = 50; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 51; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'C': case 'c':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 52; break;
                         case 'C': case 'c':                     zahl = 53; break;
                         case 'G': case 'g':                     zahl = 54; break;
                         case 'T': case 't':                     
                         case 'U': case 'u':                     zahl = 55; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'G': case 'g':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 56; break;
                         case 'C': case 'c':                     zahl = 57; break;
                         case 'G': case 'g':                     zahl = 58; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 59; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
			   break;			 
			   case 'U': case 'u':   			 
						 switch(letter3){
                         case 'A': case 'a':                     zahl = 60; break;
                         case 'C': case 'c':                     zahl = 61; break;
                         case 'G': case 'g':                     zahl = 62; break;
                         case 'T': case 't':                    
                         case 'U': case 'u':                     zahl = 63; break;
                         default: fprintf(stderr, "Illegal nucleotide symbol! %c\n",letter3);break;
                         }
               break;
		   }
	  break;
      }
      return(zahl);
}












