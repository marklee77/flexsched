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
#define EPSILON 0.0001
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

void put_service_on_server(flexsched_solution_t, int, int);
int service_can_fit_on_server(flexsched_solution_t, int, int);
double compute_allocated_resource(flexsched_solution_t, int, int);
double compute_available_resource(flexsched_solution_t, int, int);
double compute_fluid_load(flexsched_solution_t, int, int);

void initialize_global_resource_availabilities_and_loads(flexsched_problem_t);
void free_global_resource_availabilities_and_loads(flexsched_problem_t);
int service_can_fit_on_server_fast(flexsched_solution_t, int, int);
void put_service_on_server_fast(flexsched_solution_t, int, int);
double compute_available_resource_fast(flexsched_solution_t, int, int);
double compute_fluid_load_fast(flexsched_solution_t, int, int);

double compute_minimum_yield(flexsched_solution_t);
double compute_average_yield(flexsched_solution_t);
double compute_utilization(flexsched_solution_t);
int sanity_check(flexsched_solution_t);

double double_array_sum(double *, int);
double double_array_max(double *, int);
double double_array_min(double *, int);
int double_array_argmax(double *, int);
int double_array_argmin(double *, int);

void maximize_minimum_yield_on_server(flexsched_solution_t, int);
void maximize_average_yield_on_given_minimum(flexsched_solution_t, double);
void maximize_minimum_then_average_yield(flexsched_solution_t);

/*
double compute_LP_bound();
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
#define QSORT_CMP_FUNC_DECL(name) \
    int name(const void *xvp, const void *yvp)
typedef int(qsort_cmp_func)(const void *, const void *); 
extern void *qsort_thunk_vp;
#define QSORT_R(base, nel, width, thunk, compar) \
    ((compar) ? (qsort_thunk_vp = thunk, qsort(base, nel, width, compar)) : 0)
#define QSORT_CMP_CALL(compar, thunk, x, y) \
    ((compar) ? (qsort_thunk_vp = thunk, compar(x, y)) : 0)

#else
#define QSORT_CMP_FUNC_DECL(name) \
    int name(void *qsort_thunk_vp, const void *xvp, const void *yvp)
typedef int(qsort_cmp_func)(void *, const void *, const void *); 
#define QSORT_R(base, nel, width, thunk, compar) \
    ((compar) ? (qsort_r(base, nel, width, thunk, compar)) : 0)
#define QSORT_CMP_CALL(compar, thunk, x, y) \
    ((compar) ? compar(thunk, x, y) : 0)

#endif

#endif
