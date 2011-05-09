#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/****************************/
/****** VARIOUS MARCOS ******/
/****************************/


#define EPSILON 0.00001

#define RANDOM_SEED 6337

#ifndef REALLOC
#define REALLOC(ptr,size) ((size > 0) ? \
          (ptr ? realloc(ptr, size) : malloc(size)) : \
          (ptr ? (free(ptr), NULL) : NULL))
#endif

#ifndef MIN
#define MIN(x,y) ((x) > (y) ? (y) : (x))
#endif

#ifndef MAX   
#define MAX(x,y) ((x) < (y) ? (y) : (x))
#endif

#define RESOURCE_ALLOCATION_FAILURE 1
#define RESOURCE_ALLOCATION_SUCCESS 0

#define INTEGER  0
#define RATIONAL 1

#define NOT_VERBOSE 0
#define VERBOSE 1

#define SLOW_DIVE 0
#define FAST_DIVE 1

#define VP_SOLVED 1
#define VP_NOT_SOLVED 0

/******************************/
/****** GLOBAL VARIABLES ******/
/******************************/

struct flexsched_instance {
  int numrigid;
  int numfluid; 
  int numservers; 
  int numservices;
  float *slas;
  float **rigidneeds;
  float **fluidneeds;
  float **rigidcapacities;
  float **fluidcapacities;
  // output
  int *mapping;
  float *allocation;
  double lpbound;
  char misc_output[2048];
    /* String printed out at the end of the output that contains
     * misc characters that we may care about:
     *   M: The Maruyama optimization (dominance degrees) helped
     *        possibly multiple occurences for a single run
     *   T: The MILP solver has reached its time limit
     */
};

/* Global variable that keeps track of the current problem instance */
extern struct flexsched_instance INS;

/* This is a global variable that keeps track of server load
   and that is used to speed up a bunch of the naively implemented
   algorithms (especially greedy ones) */
extern float **global_server_fluid_loads;
extern float **global_server_rigid_loads;
extern float **global_server_fluidmin_loads;

/* Global for some VP stuff acceleration*/
extern float **global_vp_bin_loads;

/************************/
/****** STRUCTURES ******/
/************************/

/* A vector that needs to be packed */
struct vp_vector {
  int service;  // the service to which this vector corresponds
  int num_dims;
  float *x;
  float misc; // used to attached whatever value to a vector
};

/* A Vector Packing instance */
struct vp_instance {
  // INPUT
  int num_vectors;
  int num_dims;
  struct vp_vector *vectors;
  // OUTPUT
  int *mapping;
  int num_bins;
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
  int (*func)(char *, char *, char*);
  char *arg1;
  char *arg2;
  char *arg3;
};

extern struct scheduler_t implemented_schedulers[];

