#ifndef VECTOR_H
#define VECTOR_H

typedef struct vector_array_t {
    int num_dims;
    float **vectors;
} *vector_array;

typedef struct vp_problem_struct {
    int num_dims;
    int num_items;
    int num_bins;
    float **items;
    float **bins;
    // unlike with flex_prob, we just store the solution in the problem
    int *mapping;
    float **loads; // useful scratch variable
    float **capacities; // useful scratch variable
    float *misc; // scratch -- inelegant but useful...
} *vp_problem; 

// function prototypes
vp_problem new_vp_problem(float);
void free_vp_problem(vp_problem);

int vp_vector_can_fit_in_bin(vp_problem, int, int);
void vp_put_vector_in_bin(vp_problem, int, int);
int vp_put_vector_in_bin_safe(vp_problem, int, int);
float vp_compute_sum_load(vp_problem, int);

int cmp_int_arrays_lex(int, int [], int []);

// qsort comparison functions
int cmp_ints(const void *, const void *);

#ifdef NO_QSORT_R
extern void *global_qsort_vptr;
#define qsort_r(base, nel, width, thunk, compar) global_qsort_vptr = thunk; \
    qsort(base, nel, width, compar);
typedef (int)(const void *, const void *) qsort_cmp_func;

#else
typedef int(qsort_cmp_func)(void *, const void *, const void *);

#endif

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

qsort_cmp_func cmp_float_array_idxs;
qsort_cmp_func rcmp_float_array_idxs;

qsort_cmp_func *get_vp_cmp_func(char *);

#endif
