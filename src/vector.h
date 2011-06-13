#ifndef VECTOR_H
#define VECTOR_H

typedef struct vp_vector_s {
    double *units;
    double *totals;
} *vp_vector_t;

typedef struct vp_vector_array_s {
    int num_vectors;
    int num_dims;
    vp_vector_t *vectors;
} *vp_vector_array_t;

typedef struct vp_problem_s {
    vp_vector_array_t items;
    vp_vector_array_t bins;
} *vp_problem_t; 

typedef struct vp_solution_s {
    vp_problem_t prob;
    int *mapping;
    double **loads; // useful scratch variable
    char misc_output[100];
} *vp_solution_t;

// function prototypes
vp_problem_t new_vp_problem(flexsched_problem_t, double);
void free_vp_problem(vp_problem_t);
vp_solution_t new_vp_solution(vp_problem_t);
void free_vp_solution(vp_solution_t);

int vp_vector_can_fit_in_bin(vp_solution_t, int, int);
void vp_put_vector_in_bin(vp_solution_t, int, int);
int vp_put_vector_in_bin_safe(vp_solution_t, int, int);
double vp_compute_sum_load(vp_solution_t, int);

int cmp_int_arrays_lex(int, int [], int []);

// qsort comparison functions
int cmp_ints(const void *, const void *);

qsort_cmp_func cmp_vector_array_idxs_lex;
qsort_cmp_func rcmp_vector_array_idxs_lex;

qsort_cmp_func cmp_vector_array_idxs_max;
qsort_cmp_func rcmp_vector_array_idxs_max;

qsort_cmp_func cmp_vector_array_idxs_sum;
qsort_cmp_func rcmp_vector_array_idxs_sum;

qsort_cmp_func cmp_vector_array_idxs_maxratio;
qsort_cmp_func rcmp_vector_array_idxs_maxratio;

qsort_cmp_func cmp_vector_array_idxs_maxdiff;
qsort_cmp_func rcmp_vector_array_idxs_maxdiff;

qsort_cmp_func cmp_double_array_idxs;
qsort_cmp_func rcmp_double_array_idxs;

qsort_cmp_func *get_vp_cmp_func(char *);

#endif
