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
#define CMP(x, y) (((x) == (y)) ? 0 : (((x) < (y)) ? -1 : 1))
#endif

#ifndef RCMP // NOTE: revorse compare, NOT Royal Canadian Mounted Police!
#define RCMP(x, y) (((x) == (y)) ? 0 : (((x) < (y)) ? 1 : -1))
#endif


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

#if 0
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
#endif

/******************************/
/****** GLOBAL VARIABLES ******/
/******************************/

/* Global variable that keeps track of the current problem instance */
extern flexsched_problem flex_prob;

/*********************************/
/****** FUNCTION PROTOTYPES ******/
/*********************************/

/* utility function prototypes */
void initialize_global_server_loads();
void free_global_server_loads();
void add_service_load_to_server(int, int);
double compute_LP_bound();
float array_sum(float *, int);
float array_max(float *, int);
float array_min(float *, int);
int array_argmax(float *, int);
float compute_sum_server_load(flexsched_solution, int, const char *);
float compute_sum_server_load_fast(int, const char *);
float compute_server_load_in_dimension(
    flexsched_solution, int, const char*, int);
float compute_server_load_in_dimension_fast(int, const char*, int);
int service_can_fit_on_server(flexsched_solution, int, int);
int service_can_fit_on_server_fast(int, int);
void maximize_minimum_yield_on_server(flexsched_solution, int);
float compute_minimum_yield(flexsched_solution);
void maximize_minimum_yield(flexsched_solution);
void maximize_average_yield_on_server_given_minimum(
    flexsched_solution, int, float);
void maximize_minimum_then_average_yield(flexsched_solution);
float compute_average_yield(flexsched_solution);
float compute_server_sum_alloc(flexsched_solution, int);
float compute_utilization(flexsched_solution);
int sanity_check(flexsched_solution);
flexsched_solution new_flexsched_solution();
void free_flexsched_solution(flexsched_solution);

/*
void increment_binary_counter(struct binary_counter_t *);
void print_binary_counter(struct binary_counter_t *);
int *find_subset_of_size(vp_problem, int *, int, int);
int find_maximum_subset(vp_problem, int *, int, int **, int);
*/


/* Scheduler function prototypes */
flexsched_solution GREEDY_scheduler(char *, char **);
flexsched_solution METAGREEDY_scheduler(char *, char **);
flexsched_solution METAGREEDYLIGHT_scheduler(char *, char **);
flexsched_solution LPBOUND_scheduler(char*, char **);
flexsched_solution MILP_scheduler(char *, char **);
flexsched_solution LPROUNDING_scheduler(char *, char **);
flexsched_solution VP_scheduler(char *, char **);

flexsched_solution HVP_scheduler(char*, char*, char*);

#endif
