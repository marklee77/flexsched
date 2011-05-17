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

int cmp_vector_array_idxs_lex(const void *, const void *);
int rcmp_vector_array_idxs_lex(const void *, const void *);

int cmp_vector_array_idxs_max(const void *, const void *);
int rcmp_vector_array_idxs_max(const void *, const void *);

int cmp_vector_array_idxs_sum(const void *, const void *);
int rcmp_vector_array_idxs_sum(const void *, const void *);

int cmp_vector_array_idxs_maxratio(const void *, const void *);
int rcmp_vector_array_idxs_maxratio(const void *, const void *);

int cmp_vector_array_idxs_maxdiff(const void *, const void *);
int rcmp_vector_array_idxs_maxdiff(const void *, const void *);

int cmp_float_array_idxs(const void *, const void *);
int rcmp_float_array_idxs(const void *, const void *);

#else

int cmp_vector_array_idxs_lex(void *, const void *, const void *);
int rcmp_vector_array_idxs_lex(void *, const void *, const void *);

int cmp_vector_array_idxs_max(void *, const void *, const void *);
int rcmp_vector_array_idxs_max(void *, const void *, const void *);

int cmp_vector_array_idxs_sum(void *, const void *, const void *);
int rcmp_vector_array_idxs_sum(void *, const void *, const void *);

int cmp_vector_array_idxs_maxratio(void *, const void *, const void *);
int rcmp_vector_array_idxs_maxratio(void *, const void *, const void *);

int cmp_vector_array_idxs_maxdiff(void *, const void *, const void *);
int rcmp_vector_array_idxs_maxdiff(void *, const void *, const void *);

int cmp_vector_array_idxs_misc(void *, const void *, const void *);
int rcmp_vector_array_idxs_misc(void *, const void *, const void *);

int cmp_float_array_idxs(void *, const void *, const void *);
int rcmp_float_array_idxs(void *, const void *, const void *);

#endif


#endif
