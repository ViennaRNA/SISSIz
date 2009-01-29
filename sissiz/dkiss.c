#include "dkiss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/******************************************************************/
/*  variables for the kiss random number generator                */
/******************************************************************/
unsigned int k,m,x,y,z,w,carry,r;
/******************************************************************/
/*  kickstart the kiss random number generator                    */
/******************************************************************/
void start_kiss(int seed)
{
     x=seed;y=102;z=12;w=34535;
     x = x * 69069 + 1;
     y ^= y << 13;
     y ^= y >> 17;
     y ^= y << 5;
     k = (z >> 2) + (w >> 3) + (carry >> 2);
     m = w + w + z + carry;
     z = w;
     w = m;
     carry = k >> 30;
}
/******************************************************************/
void restart_kiss(unsigned int *vals)
{
        k=vals[0];
        m=vals[1];
        x=vals[2];
        y=vals[3];
        z=vals[4];
        w=vals[5];
        carry=vals[6];
        r=vals[7];
}
/******************************************************************/
void kiss_state(unsigned int *vals)
{
        vals[0]=k;
        vals[1]=m;
        vals[2]=x;
        vals[3]=y;
        vals[4]=z;
        vals[5]=w;
        vals[6]=carry;
        vals[7]=r;
}
/******************************************************************/
/*   Keep random number generator from George
     Marsaglia's DIEHARD cdrom                                    */
/******************************************************************/
unsigned int kiss(void)
{
     x = x * 69069 + 1;
     y ^= y << 13;
     y ^= y >> 17;
     y ^= y << 5;
     k = (z >> 2) + (w >> 3) + (carry >> 2);
     m = w + w + z + carry;
     z = w;
     w = m;
     carry = k >> 30;
     return x+y+z;
}
/******************************************************************/
double dkiss(void)
{
    return ((double)kiss()+0.5)/4294967296.0;
}
/******************************************************************/


