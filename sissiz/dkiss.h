#ifndef DKISS_H
#define DKISS_H

#define ranDum dkiss
#define starTup start_kiss

double dkiss(void);
void start_kiss(int seed);
void restart_kiss(unsigned int *vals);
void kiss_state(unsigned int *vals);

#endif


