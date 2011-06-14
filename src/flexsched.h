#ifndef FLEXSCHED_H
#define FLEXSCHED_H

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/****************************/
/**** PARAMETER CONSTANTS ***/
/****************************/
#define EPSILON 0.00001
#define RANDOM_SEED 6337

#define BEST_FIT 1
#define FIRST_FIT 0
#define WORST_FIT -1

/****************************/
/****** VARIOUS MARCOS ******/
/****************************/

#ifndef MIN
#define MIN(x,y) (((x) > (y)) ? (y) : (x))
#endif

#ifndef MAX   
#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#endif

#ifndef CMP
#define CMP(x, y) (((x) == (y)) ? 0 : (((x) < (y)) ? -1 : 1))
#endif

/*************************************/
/****** STRUCTURES AND TYPEDEFS ******/
/*************************************/

typedef struct server_s {
    double *unit_capacities;
    double *total_capacities;
} *server_t;

typedef struct service_s {
    double *unit_rigid_requirements;
    double *unit_fluid_needs;
    double *total_rigid_requirements;
    double *total_fluid_needs;
} *service_t;

typedef struct flexsched_problem_s {
    int num_resources;
    int num_servers; 
    int num_services;
    server_t *servers;
    service_t *services;
} *flexsched_problem_t;

typedef struct flexsched_solution_s {
    flexsched_problem_t prob;
    int success;
    int *mapping;
    double *yields;
    char misc_output[100];
} *flexsched_solution_t;

/*********************************/
/****** FUNCTION PROTOTYPES ******/
/*********************************/

/* utility function prototypes */
flexsched_solution_t new_flexsched_solution(flexsched_problem_t);
void free_flexsched_solution(flexsched_solution_t);

int service_can_fit_on_server(flexsched_solution_t, int, int);

void initialize_global_resource_availabilities_and_loads(flexsched_problem_t);
void free_global_resource_availabilities_and_loads(flexsched_problem_t);
void put_service_on_server(flexsched_solution_t, int, int);
int service_can_fit_on_server_fast(flexsched_solution_t, int, int);


double compute_available_resource(flexsched_solution_t, int, int);
double compute_minimum_yield(flexsched_solution_t);
double compute_average_yield(flexsched_solution_t);
double compute_utilization(flexsched_solution_t);
int sanity_check(flexsched_solution_t);

/*
double compute_LP_bound();
double array_sum(float *, int);
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
void maximize_minimum_yield(flexsched_solution);
void maximize_average_yield_on_server_given_minimum(
    flexsched_solution, int, float);
void maximize_minimum_then_average_yield(flexsched_solution);
float compute_server_sum_alloc(flexsched_solution, int);
*/

/* Scheduler function prototypes */
flexsched_solution_t GREEDY_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t METAGREEDY_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t METAGREEDYLIGHT_scheduler(flexsched_problem_t, char *, 
    char **);
flexsched_solution_t LPBOUND_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t MILP_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t LPROUNDING_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t VP_scheduler(flexsched_problem_t, char *, char **);
flexsched_solution_t HVP_scheduler(flexsched_problem_t, char *, char **);

/*
flexsched_solution METAHVP_scheduler(char *, char **);
*/

// QSORT_R
#ifdef NO_QSORT_R
#define QSORT_CMP_FUNC_DECL(name) 
    int name(const void *xvp, const void *yvp)
typedef (int)(const void *, const void *) qsort_cmp_func;
extern void *qsort_thunk_vp;
#define qsort_r(base, nel, width, thunk, compar) \
    (qsort_thunk_vp = thunk; qsort(base, nel, width, compar))
#define QSORT_CMP_CALL(cmp_items, thunk, x, y) \
    (qsort_thunk_vp = thunk, cmp_items(x, y))

#else
#define QSORT_CMP_FUNC_DECL(name) \
    int name(void *qsort_thunk_vp, const void *xvp, const void *yvp)
typedef int(qsort_cmp_func)(void *, const void *, const void *); 
#define QSORT_CMP_CALL(cmp_items, thunk, x, y) \
    cmp_items(thunk, x, y)

#endif

#endif
