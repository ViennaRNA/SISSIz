#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rnaz_utils.h"

#define PREC 8
#define PRC  100
#define LEN  1000

/* Olivier Gascuel's BIONJ */

typedef struct word{
  char name[LEN];
  struct word *suiv;
} WORD;

typedef struct pointers{
  WORD *head;
  WORD *tail;
} POINTERS;


void IInitialize(float **delta, FILE *input, int n, POINTERS *trees);
void CCompute_sums_Sx(float **delta, int n);
void BBest_pair(float **delta, int r, int *a, int *b, int n);
void FFinish(float **delta, int n, POINTERS *trees, char *treeString);
void CConcatenate(char chain1[LEN], int ind, POINTERS *trees, int post);
void PPrint_output(int i, POINTERS *trees, char *output);
float DDistance(int i, int j, float **delta);
float VVariance(int i, int j, float **delta);
float SSum_S(int i, float **delta);
float AAgglomerative_criterion(int i, int j, float **delta, int r);
float BBranch_length(int a, int b, float **delta, int r);
float RReduction4(int a, float la, int b, float lb, int i, float lamda,float **delta);
float RReduction10(int a, int b, int i, float lamda, float vab, float **delta);
float LLamda(int a, int b, float vab, float **delta, int n, int r);
float FFinish_branch_length(int i, int j, int k, float **delta);
int EEmptied(int i, float **delta);
int SSymmetrize(float **delta, int n);


void bbionj(char* treeString,const struct aln *alignment[],double* distMatrix){
  
  POINTERS *trees;
  char *chain1;
  char *chain2;
  int *a, *b;
  float **delta;
  float la;
  float lb;
  float vab;
  float lamda;
  int i;
  int ok;
  int r;
  int n;
  int x, y;
  double t;

  int lig;
  int col;
  float distance;
  char name_taxon[LEN];
  WORD *name;

  a=(int*)calloc(1,sizeof(int));
  b=(int*)calloc(1,sizeof(int));
  chain1=(char *)calloc(LEN,sizeof(char));
  chain2=(char *)calloc(LEN,sizeof(char));

  for (n=0; alignment[n]!=NULL; n++);

  delta=(float **)calloc(n+1,sizeof(float*));

  for(i=1; i<= n; i++){
    delta[i]=(float *)calloc(n+1, sizeof(float));
    if(delta[i] == NULL){
      printf("Out of memories!!");
      exit(0);
    }
  }

  trees=(POINTERS *)calloc(n+1,sizeof(POINTERS));
  if(trees == NULL){
    printf("Out of memories!!");
    exit(0);
  }
  
  r=n;
  *a=0;
  *b=0;
  
  for(lig=1; lig <= n; lig++){
    strcpy(name_taxon,alignment[lig-1]->name);
    name=(WORD *)calloc(1,sizeof(WORD));      
    if(name == NULL){                          
      printf("Out of memories !!");
      exit(0);
    } else {
      strcpy(name->name,name_taxon);
      name->suiv=NULL;
      trees[lig].head=name;
      trees[lig].tail=name;
      for(col= 1; col <= n; col++){
        distance=(float)distMatrix[(lig-1)*n+(col-1)];
        delta[lig][col]=distance;
      }
    }
  }


  ok=SSymmetrize(delta, n);
  if(!ok)
    printf("\n The matrix  is not symmetric.\n ");
  while (r > 3){                          
	  CCompute_sums_Sx(delta, n);            
	  BBest_pair(delta, r, a, b, n);         
	  vab=VVariance(*a, *b, delta);          
	  la=BBranch_length(*a, *b, delta, r);   
	  lb=BBranch_length(*b, *a, delta, r);   
	  lamda=LLamda(*a, *b, vab, delta, n, r);
	  for(i=1; i <= n; i++){
      if(!EEmptied(i,delta) && (i != *a) && (i != *b)){
        if(*a > i){
          x=*a;
          y=i;
        }
        else {
		      x=i;
		      y=*a;                           
		    }                                 
        delta[x][y]=RReduction4(*a, la, *b, lb, i, lamda, delta);
        delta[y][x]=RReduction10(*a, *b, i, lamda, vab, delta);
      }
    }
	  strcpy(chain1,"");
	  strcat(chain1,"(");
	  CConcatenate(chain1, *a, trees, 0);    
	  strcpy(chain1,"");
	  strcat(chain1,":");
	  sprintf(chain1+strlen(chain1),"%f",la);
	  strcat(chain1,",");
	  CConcatenate(chain1,*a, trees, 1);
	  trees[*a].tail->suiv=trees[*b].head;
	  trees[*a].tail=trees[*b].tail;
	  strcpy(chain1,"");
	  strcat(chain1,":");
	  sprintf(chain1+strlen(chain1),"%f",lb);
	  strcat(chain1,")");
	  CConcatenate(chain1, *a, trees, 1);
	  delta[*b][0]=1.0;                     /* make the b line empty     */
	  trees[*b].head=NULL;
	  trees[*b].tail=NULL;
	  r=r-1;                                /* decrease r                */
	}
  
  FFinish(delta, n, trees, treeString);

  for(i=1; i<=n; i++){
    delta[i][0]=0.0;
	  trees[i].head=NULL;
	  trees[i].tail=NULL;
	}

 
  for(i=1; i<= n; i++){
    free(delta[i]);
  }

  free(delta);


  free(trees);
  free(a);
  free(b);
  free(chain1);
  free(chain2);
}


void IInitialize(float **delta, FILE *input, int n, POINTERS *trees)
{
  int lig;                                          /* matrix line       */
  int col;                                          /* matrix column     */
  float distance;
  char name_taxon[LEN];                             /* taxon’s name      */
  WORD *name;
  
  for(lig=1; lig <= n; lig++)
    {
      fscanf(input,"%s",name_taxon);                  /* read taxon’s name */
      name=(WORD *)calloc(1,sizeof(WORD));            /* taxon’s name is   */
      if(name == NULL)                                /* put in trees      */
	{
	  printf("Out of memories !!");
	  exit(0);
	}
      else
	{
	  strcpy(name->name,name_taxon);
	  name->suiv=NULL;
	  trees[lig].head=name;
	  trees[lig].tail=name;
	  for(col= 1; col <= n; col++)
	    {
	      fscanf(input,"%f",&distance);             /* read the distance  */
	      delta[lig][col]=distance;
	    }
	}
    }
}


void PPrint_output(int i, POINTERS *trees, char *output)
{
  WORD *parcour;
  char tmpString[1000];
  parcour=trees[i].head;
  while(parcour != NULL)
    {
      sprintf(tmpString,"%s",parcour->name);
      strcat(output,tmpString);
      parcour=parcour->suiv;
    }
  
}

int SSymmetrize(float **delta, int n)
{
  int lig;                                         /* matrix line        */
  int col;                                         /* matrix column      */
  float value;                                     /* symmetrized value  */
  int symmetric;
  
  symmetric=1;
  for(lig=1; lig  <=  n; lig++)
    {
      for(col=1; col< lig; col++)
	{
	  if(delta[lig][col] != delta[col][lig])
	    {
	      value= (delta[lig][col]+delta[col][lig])/2;
	      delta[lig][col]=value;
	      delta[col][lig]=value;
	      symmetric=0;
	    }
        }
    }
  if(!symmetric)
    printf("The matrix is not symmetric");
  return(symmetric);
}

void CConcatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
{
  WORD *bran;
  
  bran=(WORD *)calloc(1,sizeof(WORD));
  if(bran == NULL)
    {
      printf("Out of memories");
      exit(0);
    }
  else
    {
      strcpy(bran->name,chain1);
      bran->suiv=NULL;
    }
  if(post == 0)
    {
      bran->suiv=trees[ind].head;
      trees[ind].head=bran;
    }
  else
    {
      trees[ind].tail->suiv=bran;
      trees[ind].tail=trees[ind].tail->suiv;
    }
}

float DDistance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[i][j]);
  else
    return(delta[j][i]);
}


float VVariance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[j][i]);
  else
    return(delta[i][j]);
}


int EEmptied(int i, float **delta)      /* test if the ith line is emptied */
{
  return((int)delta[i][0]);
}


float SSum_S(int i, float **delta)          /* get sum Si form the diagonal */
{
  return(delta[i][i]);
}


void CCompute_sums_Sx(float **delta, int n)
{
  float sum;
  int i;
  int j;
  
  for(i= 1; i <= n ; i++)
    {
      if(!EEmptied(i,delta))
	{
	  sum=0;
	  for(j=1; j <=n; j++)
	    {
	      if(i != j && !EEmptied(j,delta))           /* compute the sum Si */
		sum=sum + DDistance(i,j,delta);
	    }
	}
      delta[i][i]=sum;                           /* store the sum Si in */
    }                                               /* delta’s diagonal    */
}


void BBest_pair(float **delta, int r, int *a, int *b, int n)
{
  float Qxy;                         /* value of the criterion calculated*/
  int x,y;                           /* the pair which is tested         */
  float Qmin;                        /* current minimun of the criterion */
  
  Qmin=1.0e300;
  for(x=1; x <= n; x++)
    {
      if(!EEmptied(x,delta))
        {
	  for(y=1; y < x; y++)
	    {
	      if(!EEmptied(y,delta))
		{
		  Qxy=AAgglomerative_criterion(x,y,delta,r);
		  if(Qxy < Qmin-0.000001)
		    {
		      Qmin=Qxy;
		      *a=x;          
		      *b=y;
		    }
		}  
	    }
        }
    }
}

float FFinish_branch_length(int i, int j, int k, float **delta)
{
  float length;
  length=0.5*(DDistance(i,j,delta) + DDistance(i,k,delta)
	      -DDistance(j,k,delta));
  return(length);
}


void FFinish(float **delta, int n, POINTERS *trees, char *output)
{
  int l=1;
  int i=0;
  float length;
  char *str;
  WORD *bidon;
  WORD *ele;
  int last[3];                            /* the last three subtrees     */
  char tmpString[1000];
  
  str=(char *)calloc(LEN,sizeof(char));
  
  if(str == NULL)
    {
      printf("Out of memories !!");
      exit(0);
    }
  while(l <= n)
    {                                       /* find the last tree subtree  */
      if(!EEmptied(l, delta))
	{
	  last[i]=l;
	  i++;
	}
      l++;
    }
  
  length=FFinish_branch_length(last[0],last[1],last[2],delta);
  strcat(output,"(");
  PPrint_output(last[0],trees,output);
  strcat(output,":");
  sprintf(tmpString,"%f,",length);
  strcat(output,tmpString);

  length=FFinish_branch_length(last[1],last[0],last[2],delta);
  PPrint_output(last[1],trees,output);
  strcat(output,":");
  sprintf(tmpString,"%f,",length);
  strcat(output,tmpString);
    
  length=FFinish_branch_length(last[2],last[1],last[0],delta);
  PPrint_output(last[2],trees,output);
  strcat(output,":");

  sprintf(tmpString,"%f",length);
  strcat(output,tmpString);
  strcat(output,");");
  strcat(output,"\n");
  
  for(i=0; i < 3; i++)
    {
      bidon=trees[last[i]].head;
      ele=bidon;
      while(bidon!=NULL)
	{
	  ele=ele->suiv;
	  free(bidon);
	  bidon=ele;
	}
    }
  free(str);
}


float AAgglomerative_criterion(int i, int j, float **delta, int r)
{
  float Qij;
  Qij=(r-2)*DDistance(i,j,delta)                           /* Formula (1) */
    -SSum_S(i,delta)
    -SSum_S(j,delta); 
  
  return(Qij);                       
}


float BBranch_length(int a, int b, float **delta, int r)
{
  float length;
  length=0.5*(DDistance(a,b,delta)                         /* Formula (2) */
	      +(SSum_S(a,delta)
		-SSum_S(b,delta))/(r-2)); 
  return(length);                                   
}


float RReduction4(int a, float la, int b, float lb, int i, float lamda,
		 float **delta)
{
  float Dui;
  Dui=lamda*(DDistance(a,i,delta)-la)
    +(1-lamda)*(DDistance(b,i,delta)-lb);                /* Formula (4) */
  return(Dui);
}


float RReduction10(int a, int b, int i, float lamda, float vab,
		  float **delta)
{
  float Vci;
  Vci=lamda*VVariance(a,i,delta)+(1-lamda)*VVariance(b,i,delta)
    -lamda*(1-lamda)*vab;                              /*Formula (10)  */
  return(Vci);
}

float LLamda(int a, int b, float vab, float **delta, int n, int r)
{
  float lamda=0.0;
  int i;
  
  if(vab==0.0)
    lamda=0.5;
  else
    {
      for(i=1; i <= n ; i++)
	{
          if(a != i && b != i && !EEmptied(i,delta))
            lamda=lamda + (VVariance(b,i,delta) - VVariance(a,i,delta));
	}
      lamda=0.5 + lamda/(2*(r-2)*vab);             
    }                                              /* Formula (9) and the  */
  if(lamda > 1.0)                                /* constraint that lamda*/
    lamda = 1.0;                             /* belongs to [0,1]     */
  if(lamda < 0.0)
    lamda=0.0;
  return(lamda);
}





