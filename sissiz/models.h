  #ifndef MODELS_H
  #define MODELS_H

	extern   double mpi[4];
	extern   double dipi[16];
    extern   double tripi[64];
	extern   double Rmat[6];
	extern   double ratio;              /*  apha/beta in the model of HKY*/
	extern   double quo;                /*  ratio expected #transitions/#transversions*/  
	extern   int mix; 
	extern   int readfreq;
	
	
	
	extern   char   nucleotide[4];  
	/*extern   char *nucleotide[];*/
    extern   char   *doublet[];
	extern   char   *triplet[];
	extern   double trtv;
	extern   int transitiontransverion, prmat, prmat2t, prmat16;
    extern   int   cpg, matrixSH, balance;      /*only to see SH matrix: don't generate alignment with that!*/
    extern char freqFilename[256];
    extern char matrFilename[256];
    extern FILE *freq_fv;
    extern FILE *matr_fv;
	extern int verboseFreqOrder;
	extern int verboseCluster;
	extern int simind, simindependent;
	
	

	
 

    /*prototypes*/
	void checkFreq();
	void di2mono();
	void DiMono2Tri();
	
	
    double *rateRV();
	double **rateT2();    
	double **rateT2sum();
	double **rateT2CpG();     

								            

  #endif

