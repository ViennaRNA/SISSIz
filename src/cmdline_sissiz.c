/*
  File autogenerated by gengetopt version 2.19.1
  generated with the following command:
  gengetopt -F cmdline_sissiz --unamed-opts --no-handle-version --no-handle-help --no-handle-error 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getopt.h"

#include "cmdline_sissiz.h"

const char *gengetopt_args_info_purpose = "";

const char *gengetopt_args_info_usage = "Usage: " CMDLINE_PARSER_PACKAGE " [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "  -r, --nossr                   No site specific rates  (default=off)",
  "  -f, --flanks=INT              Flanking sites",
  "  -o, --outfile=STRING          Output filename",
  "  -n, --num-samples=INT         Number of samples",
  "  -m, --num-samples-regression=INT\n                                Number of samples in regression",
  "  -p, --precision=DOUBLE        Cut-off for mononucleotide content (Euclidean \n                                  distance of frequency vector)",
  "  -s, --simulate                Simulate only  (default=off)",
  "  -v, --verbose                 verbose  (default=off)",
  "  -i, --mono                    Mononucleotide  (default=off)",
  "  -d, --di                      Mononucleotide  (default=on)",
  "  -q, --dna                     Print Us or Ts  (default=on)",
  "  -u, --rna                     Print Us or Ts  (default=off)",
  "  -t, --tstv                    Transition/Transversion ratio  (default=on)",
  "  -k, --kappa=FLOAT             kappa parameter",
  "  -g, --gamma=INT               Categories for gamma for kappa ML estimation",
  "  -x, --print-rates             Print rates to file (debugging)  (default=off)",
  "      --clustal                 Output format CLUSTAL W  (default=on)",
  "      --maf                     Output format MAF  (default=off)",
  "      --fasta                   Output format FASTA  (default=off)",
    0
};

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->nossr_given = 0 ;
  args_info->flanks_given = 0 ;
  args_info->outfile_given = 0 ;
  args_info->num_samples_given = 0 ;
  args_info->num_samples_regression_given = 0 ;
  args_info->precision_given = 0 ;
  args_info->simulate_given = 0 ;
  args_info->verbose_given = 0 ;
  args_info->mono_given = 0 ;
  args_info->di_given = 0 ;
  args_info->dna_given = 0 ;
  args_info->rna_given = 0 ;
  args_info->tstv_given = 0 ;
  args_info->kappa_given = 0 ;
  args_info->gamma_given = 0 ;
  args_info->print_rates_given = 0 ;
  args_info->clustal_given = 0 ;
  args_info->maf_given = 0 ;
  args_info->fasta_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->nossr_flag = 0;
  args_info->flanks_orig = NULL;
  args_info->outfile_arg = NULL;
  args_info->outfile_orig = NULL;
  args_info->num_samples_orig = NULL;
  args_info->num_samples_regression_orig = NULL;
  args_info->precision_orig = NULL;
  args_info->simulate_flag = 0;
  args_info->verbose_flag = 0;
  args_info->mono_flag = 0;
  args_info->di_flag = 1;
  args_info->dna_flag = 1;
  args_info->rna_flag = 0;
  args_info->tstv_flag = 1;
  args_info->kappa_orig = NULL;
  args_info->gamma_orig = NULL;
  args_info->print_rates_flag = 0;
  args_info->clustal_flag = 1;
  args_info->maf_flag = 0;
  args_info->fasta_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{
  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->nossr_help = gengetopt_args_info_help[2] ;
  args_info->flanks_help = gengetopt_args_info_help[3] ;
  args_info->outfile_help = gengetopt_args_info_help[4] ;
  args_info->num_samples_help = gengetopt_args_info_help[5] ;
  args_info->num_samples_regression_help = gengetopt_args_info_help[6] ;
  args_info->precision_help = gengetopt_args_info_help[7] ;
  args_info->simulate_help = gengetopt_args_info_help[8] ;
  args_info->verbose_help = gengetopt_args_info_help[9] ;
  args_info->mono_help = gengetopt_args_info_help[10] ;
  args_info->di_help = gengetopt_args_info_help[11] ;
  args_info->dna_help = gengetopt_args_info_help[12] ;
  args_info->rna_help = gengetopt_args_info_help[13] ;
  args_info->tstv_help = gengetopt_args_info_help[14] ;
  args_info->kappa_help = gengetopt_args_info_help[15] ;
  args_info->gamma_help = gengetopt_args_info_help[16] ;
  args_info->print_rates_help = gengetopt_args_info_help[17] ;
  args_info->clustal_help = gengetopt_args_info_help[18] ;
  args_info->maf_help = gengetopt_args_info_help[19] ;
  args_info->fasta_help = gengetopt_args_info_help[20] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n", CMDLINE_PARSER_PACKAGE, CMDLINE_PARSER_VERSION);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  printf("\n%s\n\n", gengetopt_args_info_usage);

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n", gengetopt_args_info_description);

  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = NULL;
  args_info->inputs_num = 0;
}

static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  
  unsigned int i;
  if (args_info->flanks_orig)
    {
      free (args_info->flanks_orig); /* free previous argument */
      args_info->flanks_orig = 0;
    }
  if (args_info->outfile_arg)
    {
      free (args_info->outfile_arg); /* free previous argument */
      args_info->outfile_arg = 0;
    }
  if (args_info->outfile_orig)
    {
      free (args_info->outfile_orig); /* free previous argument */
      args_info->outfile_orig = 0;
    }
  if (args_info->num_samples_orig)
    {
      free (args_info->num_samples_orig); /* free previous argument */
      args_info->num_samples_orig = 0;
    }
  if (args_info->num_samples_regression_orig)
    {
      free (args_info->num_samples_regression_orig); /* free previous argument */
      args_info->num_samples_regression_orig = 0;
    }
  if (args_info->precision_orig)
    {
      free (args_info->precision_orig); /* free previous argument */
      args_info->precision_orig = 0;
    }
  if (args_info->kappa_orig)
    {
      free (args_info->kappa_orig); /* free previous argument */
      args_info->kappa_orig = 0;
    }
  if (args_info->gamma_orig)
    {
      free (args_info->gamma_orig); /* free previous argument */
      args_info->gamma_orig = 0;
    }
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);
  
  if (args_info->inputs_num)
    free (args_info->inputs);
  
  clear_given (args_info);
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  if (args_info->help_given) {
    fprintf(outfile, "%s\n", "help");
  }
  if (args_info->version_given) {
    fprintf(outfile, "%s\n", "version");
  }
  if (args_info->nossr_given) {
    fprintf(outfile, "%s\n", "nossr");
  }
  if (args_info->flanks_given) {
    if (args_info->flanks_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "flanks", args_info->flanks_orig);
    } else {
      fprintf(outfile, "%s\n", "flanks");
    }
  }
  if (args_info->outfile_given) {
    if (args_info->outfile_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "outfile", args_info->outfile_orig);
    } else {
      fprintf(outfile, "%s\n", "outfile");
    }
  }
  if (args_info->num_samples_given) {
    if (args_info->num_samples_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "num-samples", args_info->num_samples_orig);
    } else {
      fprintf(outfile, "%s\n", "num-samples");
    }
  }
  if (args_info->num_samples_regression_given) {
    if (args_info->num_samples_regression_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "num-samples-regression", args_info->num_samples_regression_orig);
    } else {
      fprintf(outfile, "%s\n", "num-samples-regression");
    }
  }
  if (args_info->precision_given) {
    if (args_info->precision_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "precision", args_info->precision_orig);
    } else {
      fprintf(outfile, "%s\n", "precision");
    }
  }
  if (args_info->simulate_given) {
    fprintf(outfile, "%s\n", "simulate");
  }
  if (args_info->verbose_given) {
    fprintf(outfile, "%s\n", "verbose");
  }
  if (args_info->mono_given) {
    fprintf(outfile, "%s\n", "mono");
  }
  if (args_info->di_given) {
    fprintf(outfile, "%s\n", "di");
  }
  if (args_info->dna_given) {
    fprintf(outfile, "%s\n", "dna");
  }
  if (args_info->rna_given) {
    fprintf(outfile, "%s\n", "rna");
  }
  if (args_info->tstv_given) {
    fprintf(outfile, "%s\n", "tstv");
  }
  if (args_info->kappa_given) {
    if (args_info->kappa_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "kappa", args_info->kappa_orig);
    } else {
      fprintf(outfile, "%s\n", "kappa");
    }
  }
  if (args_info->gamma_given) {
    if (args_info->gamma_orig) {
      fprintf(outfile, "%s=\"%s\"\n", "gamma", args_info->gamma_orig);
    } else {
      fprintf(outfile, "%s\n", "gamma");
    }
  }
  if (args_info->print_rates_given) {
    fprintf(outfile, "%s\n", "print-rates");
  }
  if (args_info->clustal_given) {
    fprintf(outfile, "%s\n", "clustal");
  }
  if (args_info->maf_given) {
    fprintf(outfile, "%s\n", "maf");
  }
  if (args_info->fasta_given) {
    fprintf(outfile, "%s\n", "fasta");
  }
  
  fclose (outfile);

  i = EXIT_SUCCESS;
  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}


/* gengetopt_strdup() */
/* strdup.c replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = NULL;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;

  result = cmdline_parser_internal (argc, argv, args_info, override, initialize, check_required, NULL);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  return EXIT_SUCCESS;
}

int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = 1;
  optopt = '?';

  while (1)
    {
      int option_index = 0;
      char *stop_char;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "nossr",	0, NULL, 'r' },
        { "flanks",	1, NULL, 'f' },
        { "outfile",	1, NULL, 'o' },
        { "num-samples",	1, NULL, 'n' },
        { "num-samples-regression",	1, NULL, 'm' },
        { "precision",	1, NULL, 'p' },
        { "simulate",	0, NULL, 's' },
        { "verbose",	0, NULL, 'v' },
        { "mono",	0, NULL, 'i' },
        { "di",	0, NULL, 'd' },
        { "dna",	0, NULL, 'q' },
        { "rna",	0, NULL, 'u' },
        { "tstv",	0, NULL, 't' },
        { "kappa",	1, NULL, 'k' },
        { "gamma",	1, NULL, 'g' },
        { "print-rates",	0, NULL, 'x' },
        { "clustal",	0, NULL, 0 },
        { "maf",	0, NULL, 0 },
        { "fasta",	0, NULL, 0 },
        { NULL,	0, NULL, 0 }
      };

      stop_char = 0;
      c = getopt_long (argc, argv, "hVrf:o:n:m:p:svidqutk:g:x", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          if (local_args_info.help_given)
            {
              fprintf (stderr, "%s: `--help' (`-h') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->help_given && ! override)
            continue;
          local_args_info.help_given = 1;
          args_info->help_given = 1;
          cmdline_parser_free (&local_args_info);
          return 0;

        case 'V':	/* Print version and exit.  */
          if (local_args_info.version_given)
            {
              fprintf (stderr, "%s: `--version' (`-V') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->version_given && ! override)
            continue;
          local_args_info.version_given = 1;
          args_info->version_given = 1;
          cmdline_parser_free (&local_args_info);
          return 0;

        case 'r':	/* No site specific rates.  */
          if (local_args_info.nossr_given)
            {
              fprintf (stderr, "%s: `--nossr' (`-r') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->nossr_given && ! override)
            continue;
          local_args_info.nossr_given = 1;
          args_info->nossr_given = 1;
          args_info->nossr_flag = !(args_info->nossr_flag);
          break;

        case 'f':	/* Flanking sites.  */
          if (local_args_info.flanks_given)
            {
              fprintf (stderr, "%s: `--flanks' (`-f') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->flanks_given && ! override)
            continue;
          local_args_info.flanks_given = 1;
          args_info->flanks_given = 1;
          args_info->flanks_arg = strtol (optarg, &stop_char, 0);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->flanks_orig)
            free (args_info->flanks_orig); /* free previous string */
          args_info->flanks_orig = gengetopt_strdup (optarg);
          break;

        case 'o':	/* Output filename.  */
          if (local_args_info.outfile_given)
            {
              fprintf (stderr, "%s: `--outfile' (`-o') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->outfile_given && ! override)
            continue;
          local_args_info.outfile_given = 1;
          args_info->outfile_given = 1;
          if (args_info->outfile_arg)
            free (args_info->outfile_arg); /* free previous string */
          args_info->outfile_arg = gengetopt_strdup (optarg);
          if (args_info->outfile_orig)
            free (args_info->outfile_orig); /* free previous string */
          args_info->outfile_orig = gengetopt_strdup (optarg);
          break;

        case 'n':	/* Number of samples.  */
          if (local_args_info.num_samples_given)
            {
              fprintf (stderr, "%s: `--num-samples' (`-n') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->num_samples_given && ! override)
            continue;
          local_args_info.num_samples_given = 1;
          args_info->num_samples_given = 1;
          args_info->num_samples_arg = strtol (optarg, &stop_char, 0);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->num_samples_orig)
            free (args_info->num_samples_orig); /* free previous string */
          args_info->num_samples_orig = gengetopt_strdup (optarg);
          break;

        case 'm':	/* Number of samples in regression.  */
          if (local_args_info.num_samples_regression_given)
            {
              fprintf (stderr, "%s: `--num-samples-regression' (`-m') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->num_samples_regression_given && ! override)
            continue;
          local_args_info.num_samples_regression_given = 1;
          args_info->num_samples_regression_given = 1;
          args_info->num_samples_regression_arg = strtol (optarg, &stop_char, 0);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->num_samples_regression_orig)
            free (args_info->num_samples_regression_orig); /* free previous string */
          args_info->num_samples_regression_orig = gengetopt_strdup (optarg);
          break;

        case 'p':	/* Cut-off for mononucleotide content (Euclidean distance of frequency vector).  */
          if (local_args_info.precision_given)
            {
              fprintf (stderr, "%s: `--precision' (`-p') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->precision_given && ! override)
            continue;
          local_args_info.precision_given = 1;
          args_info->precision_given = 1;
          args_info->precision_arg = strtod (optarg, &stop_char);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->precision_orig)
            free (args_info->precision_orig); /* free previous string */
          args_info->precision_orig = gengetopt_strdup (optarg);
          break;

        case 's':	/* Simulate only.  */
          if (local_args_info.simulate_given)
            {
              fprintf (stderr, "%s: `--simulate' (`-s') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->simulate_given && ! override)
            continue;
          local_args_info.simulate_given = 1;
          args_info->simulate_given = 1;
          args_info->simulate_flag = !(args_info->simulate_flag);
          break;

        case 'v':	/* verbose.  */
          if (local_args_info.verbose_given)
            {
              fprintf (stderr, "%s: `--verbose' (`-v') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->verbose_given && ! override)
            continue;
          local_args_info.verbose_given = 1;
          args_info->verbose_given = 1;
          args_info->verbose_flag = !(args_info->verbose_flag);
          break;

        case 'i':	/* Mononucleotide.  */
          if (local_args_info.mono_given)
            {
              fprintf (stderr, "%s: `--mono' (`-i') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->mono_given && ! override)
            continue;
          local_args_info.mono_given = 1;
          args_info->mono_given = 1;
          args_info->mono_flag = !(args_info->mono_flag);
          break;

        case 'd':	/* Mononucleotide.  */
          if (local_args_info.di_given)
            {
              fprintf (stderr, "%s: `--di' (`-d') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->di_given && ! override)
            continue;
          local_args_info.di_given = 1;
          args_info->di_given = 1;
          args_info->di_flag = !(args_info->di_flag);
          break;

        case 'q':	/* Print Us or Ts.  */
          if (local_args_info.dna_given)
            {
              fprintf (stderr, "%s: `--dna' (`-q') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->dna_given && ! override)
            continue;
          local_args_info.dna_given = 1;
          args_info->dna_given = 1;
          args_info->dna_flag = !(args_info->dna_flag);
          break;

        case 'u':	/* Print Us or Ts.  */
          if (local_args_info.rna_given)
            {
              fprintf (stderr, "%s: `--rna' (`-u') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->rna_given && ! override)
            continue;
          local_args_info.rna_given = 1;
          args_info->rna_given = 1;
          args_info->rna_flag = !(args_info->rna_flag);
          break;

        case 't':	/* Transition/Transversion ratio.  */
          if (local_args_info.tstv_given)
            {
              fprintf (stderr, "%s: `--tstv' (`-t') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->tstv_given && ! override)
            continue;
          local_args_info.tstv_given = 1;
          args_info->tstv_given = 1;
          args_info->tstv_flag = !(args_info->tstv_flag);
          break;

        case 'k':	/* kappa parameter.  */
          if (local_args_info.kappa_given)
            {
              fprintf (stderr, "%s: `--kappa' (`-k') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->kappa_given && ! override)
            continue;
          local_args_info.kappa_given = 1;
          args_info->kappa_given = 1;
          args_info->kappa_arg = (float)strtod (optarg, &stop_char);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->kappa_orig)
            free (args_info->kappa_orig); /* free previous string */
          args_info->kappa_orig = gengetopt_strdup (optarg);
          break;

        case 'g':	/* Categories for gamma for kappa ML estimation.  */
          if (local_args_info.gamma_given)
            {
              fprintf (stderr, "%s: `--gamma' (`-g') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->gamma_given && ! override)
            continue;
          local_args_info.gamma_given = 1;
          args_info->gamma_given = 1;
          args_info->gamma_arg = strtol (optarg, &stop_char, 0);
          if (!(stop_char && *stop_char == '\0')) {
            fprintf(stderr, "%s: invalid numeric value: %s\n", argv[0], optarg);
            goto failure;
          }
          if (args_info->gamma_orig)
            free (args_info->gamma_orig); /* free previous string */
          args_info->gamma_orig = gengetopt_strdup (optarg);
          break;

        case 'x':	/* Print rates to file (debugging).  */
          if (local_args_info.print_rates_given)
            {
              fprintf (stderr, "%s: `--print-rates' (`-x') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->print_rates_given && ! override)
            continue;
          local_args_info.print_rates_given = 1;
          args_info->print_rates_given = 1;
          args_info->print_rates_flag = !(args_info->print_rates_flag);
          break;


        case 0:	/* Long option with no short option */
          /* Output format CLUSTAL W.  */
          if (strcmp (long_options[option_index].name, "clustal") == 0)
          {
            if (local_args_info.clustal_given)
              {
                fprintf (stderr, "%s: `--clustal' option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
                goto failure;
              }
            if (args_info->clustal_given && ! override)
              continue;
            local_args_info.clustal_given = 1;
            args_info->clustal_given = 1;
            args_info->clustal_flag = !(args_info->clustal_flag);
          }
          /* Output format MAF.  */
          else if (strcmp (long_options[option_index].name, "maf") == 0)
          {
            if (local_args_info.maf_given)
              {
                fprintf (stderr, "%s: `--maf' option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
                goto failure;
              }
            if (args_info->maf_given && ! override)
              continue;
            local_args_info.maf_given = 1;
            args_info->maf_given = 1;
            args_info->maf_flag = !(args_info->maf_flag);
          }
          /* Output format FASTA.  */
          else if (strcmp (long_options[option_index].name, "fasta") == 0)
          {
            if (local_args_info.fasta_given)
              {
                fprintf (stderr, "%s: `--fasta' option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
                goto failure;
              }
            if (args_info->fasta_given && ! override)
              continue;
            local_args_info.fasta_given = 1;
            args_info->fasta_given = 1;
            args_info->fasta_flag = !(args_info->fasta_flag);
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}