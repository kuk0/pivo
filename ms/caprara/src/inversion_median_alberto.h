#ifdef __cplusplus
  extern "C" {
#endif


/* this is the inversion median solver implemented by Alberto
   and integrated into GRAPPA by Tao Liu.
   There are many global variables used and are initialized in
   the main function
*/
#ifndef INVERSION_MEDIAN_H
#define INVERSION_MEDIAN_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "invdist.h"

#define TWO 2
#define ONE_GOOD 1
#define ONE_BAD 0
#define TLIM 3600.0
#define MILL 1000000
#define NONE -1
#define SOME 0
#define TRUE 1
#define FALSE 0
#define VIO_EPS 0.0001


/* constant definitions */
#define MAXQ 3                  /* maximum number of permutations in the instance */
#define MAX_CYCLES 100000       /* maximum number of active cycle variables */




/* global branch-and-bound variables */



/* data structures */

typedef struct
{
    int perm;                   /* permutation associated with the half cycle */
    int len;                    /* length of the half cycle */
    int *edge;                  /* indices of the edges in the half cycle */
    double value;               /* value in the current LP */
    unsigned active;            /* flag telling whether the ass. variable is active in LP */
    unsigned recovered;         /* flag telling whether the ass. variable was rec. in LP */
    int nitslack;               /* number of consecutive iterations with variable slack in LP */
} half_cycle;



/* ********** functions   ************ */

void init_global_variables ( int, distmem_t * distmem );
void genome_init ( int i );
void median_init (  );

/* calculate the median score */
int median_distance ( int *med );

/* elaps:    total time elapsed */
/* optimal:  termination condition: TRUE = optimal sol. found */
void termin ( float elaps, int optimal );

/* mate:  permutation matching on input */
/* per:   permutation on output */
void mat_2_per ( int ngenes, int vc, int *mate, int *per );

/* pmi:  given permutation matching */
/* k:    index of permutation       */
void find_half_cycles ( int *pmi, int k, int mc, int *ncyc );

/*  lncyc:      number of cycles already formed by the partial matching and the original permutation matchings 
    ln:         local n referred to reduced problem 
    lub;        local upper bound associated with the problem 
    lcurr;      current last node visited by Hamiltonian matching         */
void msbrbb ( int lncyc, int ln, int lub, int lcurr );

/* Given 3 permutations and their size, return the median */
int *albert_inversion_median_circular ( struct genome_struct **gen,
                                        int ngenes, int *genes );
int *albert_inversion_median_noncircular ( struct genome_struct **gen,
                                           int ngenes, int *genes );

int *Find_circular_identity ( int *perm, int k );

void alloc_error ( char *s );
void stop_error ( char *s );
void free_variables (  );

int median_distance_mb ( int *med, int conserved, int common, int pair_pd, int sig);

void genome_init_mb ( int N, int conserved, int common, int sig);

int va_check(int lcurr, int lnext, int **va);

int sign_check(int lnext, int *determined_signs);

void set_determined_signs(int lnext, int **dep, int *dep_len, int *determined_signs, int **loc_determined_signs);

void unset_determined_signs(int lnext, int **dep, int *dep_len, int *determined_signs, int **loc_determined_signs);

int range_check(int element_curr, int element_next,  
		int **element_range_map, int **range_size, int *range_max, int overlap_cmp_cnt);

void set_rangemap(int element_next, int **element_range_map, int **range_size, int overlap_cmp_cnt);

void unset_rangemap(int element_next, int **element_range_map, int **range_size, int overlap_cmp_cnt);

void msbrbb_mb ( int lncyc, int ln, int lub, int lcurr, 
		int *median_cnt, int ***medians, int dynamic_restrictions,
		int conserved, int common, int pair_pd, int sign, int allmed);


/**
 * adaption of the main function of capraras median solver, 
 * - to respect conserved / common intervals
 * - get all medians
 *.
 *@param[in] passgenome input genomes
 *@param[in] ngenes number of genes per genome
 *@param[in,out] genes solution ??? i dont use it, see medians
 *@param[out] median_cnt number of found medians
 *@param[out] medians the found medians
 *@param[in] lower_bound lower bound of the median score, only solutions with a score >= lower_bound are returned,
 *	if the value is 0  the best solutions are returned
 *@param[in] upper_bound upper bound of the median score, only solutions with a score <= upper_bound are returned,
 *	if the value is 0  the best solutions are returned
 *@param[in] allmed return all medians 
 *@param[in] lazy if in the one-median-version the trivial median score is equal to the lower bound
 *	simply take this one
 *@return a median
 */
int *albert_inversion_median_noncircular_mb ( struct genome_struct **passgenome,
		int ngenes, int *genes, int *median_cnt, int ***medians, 
		int lower_bound, int upper_bound, int allmed, int lazy);


#endif

#ifdef __cplusplus
 }
#endif

