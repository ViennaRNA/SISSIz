
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

     extern int verbose, verboseMemory, quiet, relaxedPhylip, FastaAln, Clustal;
	 extern int verboseQMatrix, verboseWarning, verboseStep , verboseNeighborMatrix,verboseDiNucContent, verboseDiNucContentFinal,verboseMarkov,verboseNeighborCt, verboseEvolutiontime, verboseRate, verboseFrequencies, verboserandom, verbosewaiting, verbosewaitingcount,userSeed, verboseRandomSeed;
     extern unsigned long randomSeed; 
     extern long totalMem;
	 extern int randomgeneratordkiss;
	 extern int *randstream; /* current SPRNG random numberstream */
	 extern int seqlen;

extern double** Pij;


	/* extern double minimumvalue;*/
	 

    int    getintNumber (char *sNumber);                /*string to int*/
    double getdoubleNumber(char *sNumber) ;             /*string to double*/
    int    GetDoubleParams(int argc, char **argv, int *argn, char *pos, int numParams, double *params);
    int    doubletsTO15(char letter1, char letter2);
	int    tripletTO64(char letter1, char letter2, char letter3);
    int    nucleoTO123(char letter);
    int    GetIntParams(int argc, char **argv, int *argn, char *pos, int numParams, int *params);
	int    GetUnsignedLongParams(int argc, char **argv, int *argn, char *pos, int numParams, unsigned long *params);

	
	void *AllocMem(long n, char *name, char *func, int showInfo);
    void *CAllocMem(long n, char *name, char *func, int showInfo);
#endif /* _GLOBAL_H_ */
