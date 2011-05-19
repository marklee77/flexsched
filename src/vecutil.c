#include "flexsched.h"
#include "vector.h"

#ifdef NO_QSORT_R
void *global_qsort_vptr;
#endif

vp_problem new_vp_problem(float yield) 
{
    int i, j;
    vp_problem vp_prob;

    vp_prob = (vp_problem)calloc(1, sizeof(struct vp_problem_struct));
    vp_prob->num_dims = flex_prob->num_rigid + flex_prob->num_fluid;
    vp_prob->num_items = flex_prob->num_services;
    vp_prob->num_bins = flex_prob->num_servers;
    vp_prob->items = (float **)calloc(vp_prob->num_items, sizeof(float *));
    vp_prob->bins = 
        (float **)calloc(vp_prob->num_bins, sizeof(float *));
    vp_prob->mapping = (int *)calloc(vp_prob->num_items, sizeof(int));
    vp_prob->loads = (float **)calloc(vp_prob->num_bins, sizeof(float *));
    vp_prob->misc = NULL;
    for (i = 0; i< vp_prob->num_items; i++) {
        vp_prob->items[i] = 
            (float *)calloc(vp_prob->num_dims, sizeof(float));
        vp_prob->mapping[i] = -1;
        for (j = 0; j < flex_prob->num_rigid; j++) {
            vp_prob->items[i][j] = flex_prob->rigid_needs[i][j]; 
        }
        for (j = 0; j < flex_prob->num_fluid; j++) {
            vp_prob->items[i][flex_prob->num_rigid + j] = 
                (flex_prob->slas[i] + yield * (1.0 - flex_prob->slas[i])) * 
                flex_prob->fluid_needs[i][j];
        }
    }
    for (i = 0; i < vp_prob->num_bins; i++) {
        vp_prob->bins[i] = 
            (float *)calloc(vp_prob->num_dims, sizeof(float));
        vp_prob->loads[i] = (float *)calloc(vp_prob->num_dims, sizeof(float));
        for (j = 0; j < flex_prob->num_rigid; j++) {
            vp_prob->bins[i][j] = flex_prob->rigid_capacities[i][j];
            vp_prob->loads[i][j] = 0.0;
        }
        for (j = 0; j < flex_prob->num_fluid; j++) {
            vp_prob->bins[i][flex_prob->num_rigid + j] = 
                flex_prob->fluid_capacities[i][j];
            vp_prob->loads[i][flex_prob->num_rigid + j] = 0.0;
        }
    }
    return vp_prob;
}

void free_vp_problem(vp_problem vp_prob)
{
    int i;
    free(vp_prob->mapping);
    for (i = 0; i < vp_prob->num_items; i++) {
        free(vp_prob->items[i]);
    }
    free(vp_prob->items);
    for (i = 0; i < vp_prob->num_bins; i++) {
        free(vp_prob->bins[i]);
        free(vp_prob->loads[i]);
    }
    free(vp_prob->loads);
    free(vp_prob->misc);
    free(vp_prob->bins);
    free(vp_prob);
    return;
}

int vp_vector_can_fit_in_bin(vp_problem vp_prob, int v, int b)
{
    int i;
    for (i = 0; i < vp_prob->num_dims; i++) {
        if (vp_prob->loads[b][i] + vp_prob->items[v][i] > 
            vp_prob->bins[b][i])
            return 0;
    }
    return 1;
}

void vp_put_vector_in_bin(vp_problem vp_prob, int v, int b)
{
    int i;
    vp_prob->mapping[v] = b;
    for (i = 0; i < vp_prob->num_dims; i++) {
        vp_prob->loads[b][i] += vp_prob->items[v][i];
    }

}

int vp_put_vector_in_bin_safe(vp_problem vp_prob, int v, int b)
{
    if (!vp_vector_can_fit_in_bin(vp_prob, v, b)) return 1;
    vp_put_vector_in_bin(vp_prob, v, b);
    return 0;
}

// load computed as the sum of all dims of all objects
float vp_compute_sum_load(vp_problem vp_prob, int b)
{
  return array_sum(vp_prob->loads[b], vp_prob->num_dims);
}

// NOT sort function
// but could be..
int cmp_int_arrays_lex(int length, int a1[], int a2[]) {
    int i;
    for (i = 0; i < length; i++) {
        if (a1[i] < a2[i]) return -1;
        if (a1[i] > a2[i]) return 1;
    }
    return 0;
}

// compare lexicographically
#ifdef NO_QSORT_R
int cmp_vector_array_idxs_lex(const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)global_qsort_vptr;
#else
int cmp_vector_array_idxs_lex(void *v_ptr, const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)v_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    int i;

    // FIXME: could call cmp int arrays lex here...
    for (i = 0; i < va->num_dims; i++) {
        if (va->vectors[x][i] < va->vectors[y][i]) return -1;
        if (va->vectors[x][i] > va->vectors[y][i]) return 1;
    }
    return 0;
}

#ifdef NO_QSORT_R
int rcmp_vector_array_idxs_lex(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_lex(x_ptr, y_ptr);
}
#else
int rcmp_vector_array_idxs_lex(void *v_ptr, const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_lex(v_ptr, x_ptr, y_ptr);
}
#endif

/* max comparison of vectors */
// for descinding order sort...
#ifdef NO_QSORT_R
int cmp_vector_array_idxs_max(const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)global_qsort_vptr;
#else
int cmp_vector_array_idxs_max(void *v_ptr, const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)v_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(array_max(va->vectors[x], va->num_dims),
        array_max(va->vectors[y], va->num_dims));
}

#ifdef NO_QSORT_R
int rcmp_vector_array_idxs_max(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_max(x_ptr, y_ptr);
}
#else
int rcmp_vector_array_idxs_max(void *v_ptr, const void *x_ptr, 
    const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_max(v_ptr, x_ptr, y_ptr);
}
#endif

/* Sum comparison of vectors */
// for descinding order sort...
#ifdef NO_QSORT_R
int cmp_vector_array_idxs_sum(const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)global_qsort_vptr;
#else
int cmp_vector_array_idxs_sum(void *v_ptr, const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)v_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(array_sum(va->vectors[x], va->num_dims), 
        array_sum(va->vectors[y], va->num_dims));
}

#ifdef NO_QSORT_R
int rcmp_vector_array_idxs_sum(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_sum(x_ptr, y_ptr);
}
#else
int rcmp_vector_array_idxs_sum(void * v_ptr, const void *x_ptr, 
    const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_sum(v_ptr, x_ptr, y_ptr);
}
#endif

#ifdef NO_QSORT_R
int cmp_vector_array_idxs_maxratio(const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)global_qsort_vptr;
#else
int cmp_vector_array_idxs_maxratio(void *v_ptr, const void *x_ptr, 
    const void *y_ptr)
{
    vector_array va = (vector_array)v_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(array_max(va->vectors[x], va->num_dims) / 
        array_min(va->vectors[x], va->num_dims), 
        array_max(va->vectors[y], va->num_dims) / 
        array_min(va->vectors[y], va->num_dims));
}

#ifdef NO_QSORT_R
int rcmp_vector_array_idxs_maxratio(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_maxratio(x_ptr, y_ptr);
}
#else
int rcmp_vector_array_idxs_maxratio(void * v_ptr, const void *x_ptr, 
    const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_maxratio(v_ptr, x_ptr, y_ptr);
}
#endif

#ifdef NO_QSORT_R
int cmp_vector_array_idxs_maxdiff(const void *x_ptr, const void *y_ptr)
{
    vector_array va = (vector_array)global_qsort_vptr;
#else
int cmp_vector_array_idxs_maxdiff(void *v_ptr, const void *x_ptr, 
        const void *y_ptr)
{
    vector_array va = (vector_array)v_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(array_max(va->vectors[x], va->num_dims) - 
        array_min(va->vectors[x], va->num_dims),
        array_max(va->vectors[y], va->num_dims) - 
        array_min(va->vectors[y], va->num_dims));
}

#ifdef NO_QSORT_R
int rcmp_vector_array_idxs_maxdiff(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_maxdiff(x_ptr, y_ptr);
}
#else
int rcmp_vector_array_idxs_maxdiff(void * v_ptr, const void *x_ptr, 
    const void *y_ptr)
{
    return -1 * cmp_vector_array_idxs_maxdiff(v_ptr, x_ptr, y_ptr);
}
#endif

/* Misc comparison of vectors  */
/* By decreasing order of MISC */
// FIXME: only really used by Maruyama, which is not currently implemented
#ifdef NO_QSORT_R
int rcmp_vp_vector_idxs_misc(const void *x_ptr, const void *y_ptr)
{
    vp_problem vp_prob = (vp_problem)global_qsort_vptr;
#else
int rcmp_vp_vector_idxs_misc(void *vp_ptr, const void *x_ptr, const void *y_ptr)
{
    vp_problem vp_prob = (vp_problem)vp_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(vp_prob->misc[x], vp_prob->misc[y]);
}

int cmp_ints(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);
    return CMP(x, y);
}

#ifdef NO_QSORT_R
int cmp_float_array_idxs(const void *x_ptr, const void *y_ptr)
{
    float *vals = (float *)global_qsort_vptr;
#else
int cmp_float_array_idxs(void *vals_ptr, const void *x_ptr, const void *y_ptr)
{
    float *vals = (float *)vals_ptr;
#endif
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(vals[x], vals[y]);
}

#ifdef NO_QSORT_R
int rcmp_float_array_idxs(const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_float_array_idxs(x_ptr, y_ptr);
}
#else
int rcmp_float_array_idxs(void *vals_ptr, const void *x_ptr, const void *y_ptr)
{
    return -1 * cmp_float_array_idxs(vals_ptr, x_ptr, y_ptr);
}
#endif

