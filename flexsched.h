#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef FLEXSCHED_H
#define FLEXSCHED_H

/****************************/
/**** PARAMETER CONSTANTS ***/
/****************************/
#define EPSILON 0.00001

#define RANDOM_SEED 6337

#define GLPK_TIME_LIMIT 10*60*1000 // 10 minutes in milliseconds


/****************************/
/****** VARIOUS MARCOS ******/
/****************************/

#ifndef REALLOC
#define REALLOC(ptr,size) \
  (((size) > 0) ? realloc(ptr, size) : (free(ptr), NULL))
#endif

#ifndef MIN
#define MIN(x,y) (((x) > (y)) ? (y) : (x))
#endif

#ifndef MAX   
#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#endif

#ifndef CMP
#define CMP(x, y) (((x) == (y)) ? 0 : (((x) < (y)) ? 1 : -1))
#endif

#define INTEGER  0
#define RATIONAL 1

/*************************************/
/****** STRUCTURES AND TYPEDEFS ******/
/*************************************/

typedef struct flexsched_problem_struct {
    int num_rigid;
    int num_fluid; 
    int num_servers; 
    int num_services;
    float **rigid_capacities;
    float **fluid_capacities;
    float *slas;
    float **rigid_needs;
    float **fluid_needs;
    // property of the problem, but output...
    double lpbound;
} *flexsched_problem;

typedef struct flexsched_solution_struct {
    flexsched_problem prob;
    // output
    char algorithm[25];
    int success;
    int *mapping;
    float *scaled_yields;
    char misc_output[100];
    /* String printed out at the end of the output that contains
     * misc characters that we may care about:
     *   M: The Maruyama optimization (dominance degrees) helped
     *        possibly multiple occurences for a single run
     *   T: The MILP solver has reached its time limit
     */
} *flexsched_solution;

/* A vector that needs to be packed */
// FIXME: do we need the extraneous info here?
struct vp_vector {
  int service;  // the service to which this vector corresponds
  int num_dims;
  float *x;
  float misc; // used to attached whatever value to a vector
};

/* A Vector Packing instance */
struct vp_instance {
  // INPUT
  int num_dims;
  int num_vectors;
  int num_bins;
  struct vp_vector *vectors;
  struct vp_vector *bins;
  // OUTPUT
  int *mapping;
}; 

/* A binary counter */
struct binary_counter_t {
  int size;
  char *bits;
  int sum;	// keeps the sum current
};

/* An int and a float */
struct int_float {
  int i;
  float f;
};

/* Data structure describing a scheduling algorithm implementation */
struct scheduler_t {
    char *name;
    flexsched_solution (*func)(char *, char *, char*);
    char *arg1;
    char *arg2;
    char *arg3;
};

/******************************/
/****** GLOBAL VARIABLES ******/
/******************************/

/* Global variable that keeps track of the current problem instance */
extern flexsched_problem flex_prob;

/*********************************/
/****** FUNCTION PROTOTYPES ******/
/*********************************/

/* utility function prototypes */
flexsched_solution new_flexsched_solution(char *);
void free_flexsched_solution(flexsched_solution);
void maximize_average_yield(flexsched_solution);
float compute_minimum_yield(flexsched_solution);
float compute_average_yield(flexsched_solution);

/* Scheduler function prototypes */
flexsched_solution LPBOUND_scheduler(char*, char*, char*);
flexsched_solution MILP_scheduler(char*, char*, char*);
flexsched_solution RRND_scheduler(char*, char*, char*);
flexsched_solution RRNZ_scheduler(char*, char*, char*);
flexsched_solution GREEDY_scheduler(char*, char*, char*);
flexsched_solution METAGREEDY_scheduler(char*, char*, char*);
flexsched_solution METAGREEDYLIGHT_scheduler(char*, char*, char*);
flexsched_solution VP_scheduler(char*, char*, char*);
flexsched_solution GA_scheduler(char*, char*, char*);

#endif
