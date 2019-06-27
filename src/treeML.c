#include "phyml/spr.h"
#include "phyml/utilities.h"
#include "phyml/lk.h"
#include "phyml/optimiz.h"
#include "phyml/bionj.h"
#include "phyml/models.h"
#include "phyml/free.h"
#include "phyml/options.h"
#include "phyml/simu.h"
#include "phyml/eigen.h"
#include "phyml/pars.h"
#include "phyml/alrt.h"
#include "rnaz_utils.h"

seq **Get_Seq_local(seq**, option *io,  int rw);

void treeML(const struct aln *alignment[], int catGamma, char** treeString, double* kappa){

  int i;
  seq **data;
  allseq *alldata;
  option *io;
  arbre *tree;
  int n_otu;
  matrix *mat;
  model *mod;
  time_t t_beg,t_end;
  div_t hour,min;
  phydbl best_lnL;
  int L,r_seed;
  int pairwiseAln=0;

  /* Initialize data structures */

  tree = NULL;
  mod  = NULL;
  data = NULL;
  best_lnL = UNLIKELY;

  io = (option *)Make_Input();
  Set_Defaults_Input(io);
  Set_Defaults_Model(io->mod);
  Set_Defaults_Optimiz(io->mod->s_opt);
  r_seed = time(NULL);
  srand(r_seed);

  /* Set options */

  io->mod->datatype=NT; /* Nucleotides not amino acids */
  io->mod->whichmodel = HKY85;
  io->mod->s_opt->print=0; /* Shut of verbose output */
  io->mod->s_opt->opt_topo = 0; /* Do not optimize topology */
  io->mod->s_opt->opt_bl   = 1; /* Optimize branch lengths*/
  /* Optimize kappa, initialize with 4.0*/
  io->mod->s_opt->opt_kappa = 1;
  io->mod->kappa = 4.0;
  
  if (catGamma>0){
    io->mod->s_opt->opt_alpha     = 1;
    io->mod->n_catg = catGamma;
  }


  io->mod->s_opt->opt_num_param = 1;

  mod = io->mod;

  Make_Model_Complete(io->mod);
		   
  /* Manually fill data structure with the alignment */
  
  n_otu = 0;

  L=strlen(alignment[0]->seq);
  for (n_otu=0; alignment[n_otu]!=NULL; n_otu++);

  /* Currently PHYML fails for pairwise alignments. As a simple
   workaround we duplicate one sequence to get three. This does not
   change the calculation of the kappa which we are mainly interested
   here.. */

  if (n_otu==2){
    pairwiseAln=1;
  }


  io->mod->n_otu=n_otu;

  if (pairwiseAln){
    data = (seq **)mCalloc(3,sizeof(seq *));
  } else {
    data = (seq **)mCalloc(n_otu,sizeof(seq *));
  }

  for (i=0;i<n_otu;i++){
    data[i] = (seq *)mCalloc(1,sizeof(seq));
    data[i]->len = L;
    data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    strcpy(data[i]->name,alignment[i]->name);
    data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
    strcpy(data[i]->state,alignment[i]->seq);
    data[i]->is_ambigu = NULL;
  }

  if (pairwiseAln){

    data[2] = (seq *)mCalloc(1,sizeof(seq));
    data[2]->len = L;
    data[2]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    strcpy(data[2]->name,"duplicated_dummy");
    data[2]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
    strcpy(data[2]->state,alignment[1]->seq);
    data[2]->is_ambigu = NULL;
    
    n_otu++;
    io->mod->n_otu=n_otu;

  }

  /* Call modified version of Get_Seq, that does some processing */
  data=Get_Seq_local(data, io, 0);
  alldata = Compact_Seq(data,io);
  Free_Seq(data,alldata->n_otu);
  Check_Ambiguities(alldata,io->mod->datatype,io->mod->stepsize);
  Init_Model(alldata,mod);


  /* Calculate pairwise distances and make BIONJ tree*/
  mat = ML_Dist(alldata,mod);
  Fill_Missing_Dist(mat);
  mat->tree = Make_Tree_From_Scratch(alldata->n_otu,alldata);
  Bionj(mat);
  tree = mat->tree;
  tree->mat = mat;
  tree->mod = mod;
  tree->io = io;
  tree->data = alldata;
  tree->both_sides = 1;
  tree->n_pattern = tree->data->crunch_len/tree->mod->stepsize;
  time(&t_beg);
  time(&(tree->t_beg));

  /* Print_Settings(io);*/

  /* Prepare for optimization and optimize */

  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,alldata,alldata->init_len);
  Make_Tree_4_Lk(tree,alldata,alldata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  Br_Len_Not_Involving_Invar(tree);
  Order_Tree_CSeq(tree,alldata);
  Round_Optimize(tree,tree->data);

  //treeString=(char *)mCalloc(T_MAX_LINE,sizeof(char));
  
  *treeString=Write_Tree(tree);
  *kappa=io->mod->kappa;

}

seq **Get_Seq_local(seq** data, option *io,  int rw)
{
  int i,j;
  char **buff;
  int n_unkn,n_removed,pos;
  int *remove;


/*   rewind(fp_seq); */


  if(data)
    {
      buff = (char **)mCalloc(io->mod->n_otu,sizeof(char *));
      For(i,io->mod->n_otu) buff[i] = (char *)mCalloc(data[0]->len,sizeof(char));
      remove = (int *)mCalloc(data[0]->len,sizeof(int));

      n_removed = 0;

      For(i,data[0]->len)
	{
	  For(j,io->mod->n_otu)
	    {
	      if((data[j]->state[i] == '?') || (data[j]->state[i] == '-')) data[j]->state[i] = 'X';
	      if((io->mod->datatype == NT) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
	      if(data[j]->state[i] == 'U') data[j]->state[i] = 'T';
	    }

	  n_unkn = 0;
	  For(j,io->mod->n_otu) if(data[j]->state[i] == 'X') n_unkn++;

	  if(n_unkn == io->mod->n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }

	  For(j,io->mod->n_otu) buff[j][i] = data[j]->state[i];
	}

      if(n_removed > 0)
	{
	  if(io->mod->datatype == NT)
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n",n_removed);
	  else
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?')...\n",n_removed);
	}

      pos = 0;
      For(i,data[0]->len)
	{
/* 	  if(!remove[i]) */
/* 	    { */
	      For(j,io->mod->n_otu) data[j]->state[pos] = buff[j][i];
	      pos++;
/* 	    } */
	}

      For(i,io->mod->n_otu) data[i]->len = pos;
      For(i,io->mod->n_otu) Free(buff[i]);
      Free(buff);
      Free(remove);
    }
  return data;
}
