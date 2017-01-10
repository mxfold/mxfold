/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.6
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "ngsfold"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "ngsfold"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.0.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *full_help_help; /**< @brief Print help, including hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int noncomplementary_flag;	/**< @brief Allow non-canonical base pairs (default=off).  */
  const char *noncomplementary_help; /**< @brief Allow non-canonical base pairs help description.  */
  char * param_arg;	/**< @brief Load parameters from parameter-file.  */
  char * param_orig;	/**< @brief Load parameters from parameter-file original value given at command line.  */
  const char *param_help; /**< @brief Load parameters from parameter-file help description.  */
  int random_seed_arg;	/**< @brief Specify the seed of the random number generator (default='-1').  */
  char * random_seed_orig;	/**< @brief Specify the seed of the random number generator original value given at command line.  */
  const char *random_seed_help; /**< @brief Specify the seed of the random number generator help description.  */
  int verbose_arg;	/**< @brief Verbose output (default='0').  */
  char * verbose_orig;	/**< @brief Verbose output original value given at command line.  */
  const char *verbose_help; /**< @brief Verbose output help description.  */
  int predict_flag;	/**< @brief Prediction mode (default=on).  */
  const char *predict_help; /**< @brief Prediction mode help description.  */
  float* mea_arg;	/**< @brief MEA decoding with gamma (default='6.0').  */
  char ** mea_orig;	/**< @brief MEA decoding with gamma original value given at command line.  */
  unsigned int mea_min; /**< @brief MEA decoding with gamma's minimum occurreces */
  unsigned int mea_max; /**< @brief MEA decoding with gamma's maximum occurreces */
  const char *mea_help; /**< @brief MEA decoding with gamma help description.  */
  float* gce_arg;	/**< @brief Generalized centroid decoding with gamma (default='4.0').  */
  char ** gce_orig;	/**< @brief Generalized centroid decoding with gamma original value given at command line.  */
  unsigned int gce_min; /**< @brief Generalized centroid decoding with gamma's minimum occurreces */
  unsigned int gce_max; /**< @brief Generalized centroid decoding with gamma's maximum occurreces */
  const char *gce_help; /**< @brief Generalized centroid decoding with gamma help description.  */
  int bpseq_flag;	/**< @brief Output predicted results as the BPSEQ format (default=off).  */
  const char *bpseq_help; /**< @brief Output predicted results as the BPSEQ format help description.  */
  int constraints_flag;	/**< @brief Use contraints (default=off).  */
  const char *constraints_help; /**< @brief Use contraints help description.  */
  char * train_arg;	/**< @brief Trainining mode (write the trained parameters into output-file).  */
  char * train_orig;	/**< @brief Trainining mode (write the trained parameters into output-file) original value given at command line.  */
  const char *train_help; /**< @brief Trainining mode (write the trained parameters into output-file) help description.  */
  int max_iter_arg;	/**< @brief The maximum number of iterations for training (default='100').  */
  char * max_iter_orig;	/**< @brief The maximum number of iterations for training original value given at command line.  */
  const char *max_iter_help; /**< @brief The maximum number of iterations for training help description.  */
  int burn_in_arg;	/**< @brief The number of iterations for initial training from labeled data (default='10').  */
  char * burn_in_orig;	/**< @brief The number of iterations for initial training from labeled data original value given at command line.  */
  const char *burn_in_help; /**< @brief The number of iterations for initial training from labeled data help description.  */
  float weight_weak_label_arg;	/**< @brief The weight for weak labeled data (default='0.1').  */
  char * weight_weak_label_orig;	/**< @brief The weight for weak labeled data original value given at command line.  */
  const char *weight_weak_label_help; /**< @brief The weight for weak labeled data help description.  */
  char ** structure_arg;	/**< @brief The lists of training data with full structures.  */
  char ** structure_orig;	/**< @brief The lists of training data with full structures original value given at command line.  */
  unsigned int structure_min; /**< @brief The lists of training data with full structures's minimum occurreces */
  unsigned int structure_max; /**< @brief The lists of training data with full structures's maximum occurreces */
  const char *structure_help; /**< @brief The lists of training data with full structures help description.  */
  char ** unpaired_reactivity_arg;	/**< @brief The lists of training data with unpaired reactivity.  */
  char ** unpaired_reactivity_orig;	/**< @brief The lists of training data with unpaired reactivity original value given at command line.  */
  unsigned int unpaired_reactivity_min; /**< @brief The lists of training data with unpaired reactivity's minimum occurreces */
  unsigned int unpaired_reactivity_max; /**< @brief The lists of training data with unpaired reactivity's maximum occurreces */
  const char *unpaired_reactivity_help; /**< @brief The lists of training data with unpaired reactivity help description.  */
  char ** paired_reactivity_arg;	/**< @brief The lists of training data with paired reactivity.  */
  char ** paired_reactivity_orig;	/**< @brief The lists of training data with paired reactivity original value given at command line.  */
  unsigned int paired_reactivity_min; /**< @brief The lists of training data with paired reactivity's minimum occurreces */
  unsigned int paired_reactivity_max; /**< @brief The lists of training data with paired reactivity's maximum occurreces */
  const char *paired_reactivity_help; /**< @brief The lists of training data with paired reactivity help description.  */
  float eta_arg;	/**< @brief Initial step width for the subgradient optimization (default='0.5').  */
  char * eta_orig;	/**< @brief Initial step width for the subgradient optimization original value given at command line.  */
  const char *eta_help; /**< @brief Initial step width for the subgradient optimization help description.  */
  float pos_w_arg;	/**< @brief The weight for positive base-pairs (default='4').  */
  char * pos_w_orig;	/**< @brief The weight for positive base-pairs original value given at command line.  */
  const char *pos_w_help; /**< @brief The weight for positive base-pairs help description.  */
  float neg_w_arg;	/**< @brief The weight for negative base-pairs (default='1').  */
  char * neg_w_orig;	/**< @brief The weight for negative base-pairs original value given at command line.  */
  const char *neg_w_help; /**< @brief The weight for negative base-pairs help description.  */
  float lambda_arg;	/**< @brief The weight for the L1 regularization term (default='0.125').  */
  char * lambda_orig;	/**< @brief The weight for the L1 regularization term original value given at command line.  */
  const char *lambda_help; /**< @brief The weight for the L1 regularization term help description.  */
  float scale_reactivity_arg;	/**< @brief The scale of reactivity (default='0.1').  */
  char * scale_reactivity_orig;	/**< @brief The scale of reactivity original value given at command line.  */
  const char *scale_reactivity_help; /**< @brief The scale of reactivity help description.  */
  float threshold_unpaired_reactivity_arg;	/**< @brief The threshold of reactiviy for unpaired bases (default='0.7').  */
  char * threshold_unpaired_reactivity_orig;	/**< @brief The threshold of reactiviy for unpaired bases original value given at command line.  */
  const char *threshold_unpaired_reactivity_help; /**< @brief The threshold of reactiviy for unpaired bases help description.  */
  float threshold_paired_reactivity_arg;	/**< @brief The threshold of reactiviy for paired bases (default='0.7').  */
  char * threshold_paired_reactivity_orig;	/**< @brief The threshold of reactiviy for paired bases original value given at command line.  */
  const char *threshold_paired_reactivity_help; /**< @brief The threshold of reactiviy for paired bases help description.  */
  int discretize_reactivity_flag;	/**< @brief Discretize reactivity with reactivity thresholds (default=off).  */
  const char *discretize_reactivity_help; /**< @brief Discretize reactivity with reactivity thresholds help description.  */
  char * out_param_arg;	/**< @brief Output parameter sets for each step.  */
  char * out_param_orig;	/**< @brief Output parameter sets for each step original value given at command line.  */
  const char *out_param_help; /**< @brief Output parameter sets for each step help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int full_help_given ;	/**< @brief Whether full-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int noncomplementary_given ;	/**< @brief Whether noncomplementary was given.  */
  unsigned int param_given ;	/**< @brief Whether param was given.  */
  unsigned int random_seed_given ;	/**< @brief Whether random-seed was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int predict_given ;	/**< @brief Whether predict was given.  */
  unsigned int mea_given ;	/**< @brief Whether mea was given.  */
  unsigned int gce_given ;	/**< @brief Whether gce was given.  */
  unsigned int bpseq_given ;	/**< @brief Whether bpseq was given.  */
  unsigned int constraints_given ;	/**< @brief Whether constraints was given.  */
  unsigned int train_given ;	/**< @brief Whether train was given.  */
  unsigned int max_iter_given ;	/**< @brief Whether max-iter was given.  */
  unsigned int burn_in_given ;	/**< @brief Whether burn-in was given.  */
  unsigned int weight_weak_label_given ;	/**< @brief Whether weight-weak-label was given.  */
  unsigned int structure_given ;	/**< @brief Whether structure was given.  */
  unsigned int unpaired_reactivity_given ;	/**< @brief Whether unpaired-reactivity was given.  */
  unsigned int paired_reactivity_given ;	/**< @brief Whether paired-reactivity was given.  */
  unsigned int eta_given ;	/**< @brief Whether eta was given.  */
  unsigned int pos_w_given ;	/**< @brief Whether pos-w was given.  */
  unsigned int neg_w_given ;	/**< @brief Whether neg-w was given.  */
  unsigned int lambda_given ;	/**< @brief Whether lambda was given.  */
  unsigned int scale_reactivity_given ;	/**< @brief Whether scale-reactivity was given.  */
  unsigned int threshold_unpaired_reactivity_given ;	/**< @brief Whether threshold-unpaired-reactivity was given.  */
  unsigned int threshold_paired_reactivity_given ;	/**< @brief Whether threshold-paired-reactivity was given.  */
  unsigned int discretize_reactivity_given ;	/**< @brief Whether discretize-reactivity was given.  */
  unsigned int out_param_given ;	/**< @brief Whether out-param was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];
/** @brief all the lines making the full help output (including hidden options) */
extern const char *gengetopt_args_info_full_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the full help (including hidden options)
 */
void cmdline_parser_print_full_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
