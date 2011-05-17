
#ifndef VECTOR_H
#define VECTOR_H

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
int rcmp_vp_vector_idxs_lex(const void *, const void *);
int rcmp_vp_vector_idxs_max(const void *, const void *);
int rcmp_vp_vector_idxs_sum(const void *, const void *);
int rcmp_vp_vector_idxs_maxratio(const void *, const void *);
int rcmp_vp_vector_idxs_maxdiff(const void *, const void *);
int rcmp_vp_vector_idxs_misc(const void *, const void *);
int cmp_float_array_idxs(const void *, const void *);
int rcmp_float_array_idxs(const void *, const void *);
#else
int rcmp_vp_vector_idxs_lex(void *, const void *, const void *);
int rcmp_vp_vector_idxs_max(void *, const void *, const void *);
int rcmp_vp_vector_idxs_sum(void *, const void *, const void *);
int rcmp_vp_vector_idxs_maxratio(void *, const void *, const void *);
int rcmp_vp_vector_idxs_maxdiff(void *, const void *, const void *);
int rcmp_vp_vector_idxs_misc(void *, const void *, const void *);
int cmp_float_array_idxs(void *, const void *, const void *);
int rcmp_float_array_idxs(void *, const void *, const void *);
#endif


#endif