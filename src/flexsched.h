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
    float *unit_capacities;
    float *total_capacities;
} *server_t;

typedef struct service_s {
    float *unit_rigid_requirements;
    float *unit_fluid_needs;
    float *total_rigid_requirements;
    float *total_fluid_needs;
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
    float *yields;
    char misc_output[100];
} *flexsched_solution_t;

/*********************************/
/****** FUNCTION PROTOTYPES ******/
/*********************************/

/* utility function prototypes */
float compute_minimum_yield(flexsched_solution_t);
float compute_average_yield(flexsched_solution_t);
float compute_utilization(flexsched_solution_t);

/*
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
void maximize_minimum_yield(flexsched_solution);
void maximize_average_yield_on_server_given_minimum(
    flexsched_solution, int, float);
void maximize_minimum_then_average_yield(flexsched_solution);
float compute_server_sum_alloc(flexsched_solution, int);
int sanity_check(flexsched_solution);
flexsched_solution new_flexsched_solution();
void free_flexsched_solution(flexsched_solution);
*/

/* Scheduler function prototypes */
/*
flexsched_solution GREEDY_scheduler(flexsched_problem, char *, char **);
flexsched_solution METAGREEDY_scheduler(flexsched_problem, char *, char **);
flexsched_solution METAGREEDYLIGHT_scheduler(char *, char **);
flexsched_solution LPBOUND_scheduler(char*, char **);
flexsched_solution MILP_scheduler(char *, char **);
flexsched_solution LPROUNDING_scheduler(char *, char **);
flexsched_solution VP_scheduler(char *, char **);
flexsched_solution HVP_scheduler(char *, char **);
flexsched_solution METAHVP_scheduler(char *, char **);
*/

#endif
