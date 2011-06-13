#include "flexsched.h"
#include "vector.h"

vp_problem_t new_vp_problem(flexsched_problem_t flex_prob, double yield) 
{
    int i, j;
    vp_problem_t vp_prob = (vp_problem_t)malloc(sizeof(struct vp_problem_s));
    vp_prob->num_items = flex_prob->num_services;
    vp_prob->items = 
        (vp_vector_t *)calloc(vp_prob->num_items, sizeof(vp_vector_t));
    for (i = 0; i< vp_prob->num_items; i++) {
        vp_prob->items[i] = (vp_vector_t)malloc(sizeof(struct vp_vector_s));
        vp_prob->items[i]->num_dims = flex_prob->num_resources;
        vp_prob->items[i]->units = 
            (double *)calloc(vp_prob->items[i]->num_dims, sizeof(double));
        vp_prob->items[i]->totals = 
            (double *)calloc(vp_prob->items[i]->num_dims, sizeof(double));
        for (j = 0; j < vp_prob->items[i]->num_dims; j++) {
            vp_prob->items[i]->units[j] = 
                flex_prob->services[i]->unit_rigid_requirements[j] +
                yield*flex_prob->services[i]->unit_fluid_needs[j];
            vp_prob->items[i]->totals[j] =
                flex_prob->services[i]->total_rigid_requirements[j] +
                yield*flex_prob->services[i]->total_fluid_needs[j];
        }
    }
    vp_prob->num_bins = flex_prob->num_servers;
    vp_prob->bins = 
        (vp_vector_t *)calloc(vp_prob->num_bins, sizeof(vp_vector_t));
    for (i = 0; i< vp_prob->num_bins; i++) {
        vp_prob->bins[i] = (vp_vector_t)malloc(sizeof(struct vp_vector_s));
        vp_prob->bins[i]->num_dims = flex_prob->num_resources;
        vp_prob->bins[i]->units = 
            (double *)calloc(vp_prob->bins[i]->num_dims, sizeof(double));
        vp_prob->bins[i]->totals = 
            (double *)calloc(vp_prob->bins[i]->num_dims, sizeof(double));
        for (j = 0; j < vp_prob->bins[i]->num_dims; j++) {
            vp_prob->bins[i]->units[j] = 
                flex_prob->servers[i]->unit_capacities[j];
            vp_prob->bins[i]->totals[j] = 
                flex_prob->servers[i]->total_capacities[j];
        }
    }

    return vp_prob;
}

void free_vp_problem(vp_problem_t vp_prob)
{
    int i;
    for (i = 0; i < vp_prob->num_items; i++) {
        free(vp_prob->items[i]->units);
        free(vp_prob->items[i]->totals);
        free(vp_prob->items[i]);
    }
    free(vp_prob->items);
    for (i = 0; i < vp_prob->num_bins; i++) {
        free(vp_prob->bins[i]->units);
        free(vp_prob->bins[i]->totals);
        free(vp_prob->bins[i]);
    }
    free(vp_prob->bins);
    free(vp_prob);
    return;
}

vp_solution_t new_vp_solution(vp_problem_t vp_prob)
{
    int i, j;
    vp_solution_t vp_soln = (vp_solution_t)malloc(sizeof(struct vp_solution_s));
    vp_soln->prob = vp_prob;
    vp_soln->success = 0;
    vp_soln->mapping = (int *)calloc(vp_prob->num_items, sizeof(int));
    for (i = 0; i< vp_prob->num_items; i++) {
        vp_soln->mapping[i] = -1;
    }
    vp_soln->loads = (double **)calloc(vp_prob->num_bins, sizeof(double *));
    for (i = 0; i < vp_prob->num_bins; i++) {
        vp_soln->loads[i] = (double *)calloc(vp_prob->bins[i]->num_dims, 
            sizeof(double));
        for (j = 0; j < vp_prob->bins[i]->num_dims; j++) {
            vp_soln->loads[i][j] = 0.0;
        }
    }

    return vp_soln;
}

void free_vp_solution(vp_solution_t vp_soln)
{
    int i;
    for (i = 0; i < vp_soln->prob->num_bins; i++) {
        free(vp_soln->loads[i]);
    }
    free(vp_soln->loads);
    free(vp_soln->mapping);
    free(vp_soln);
    return;
}

int vp_item_can_fit_in_bin(vp_solution_t vp_soln, int i, int b)
{
    int j;
    for (j = 0; j < vp_soln->prob->bins[b]->num_dims; j++) {
        if (vp_soln->prob->bins[b]->units[j] < 
            vp_soln->prob->items[i]->units[j] || 
            vp_soln->prob->bins[b]->totals[j] <
            vp_soln->loads[b][j] + vp_soln->prob->items[i]->totals[j]) return 0;
    }
    return 1;
}

void vp_put_item_in_bin(vp_solution_t vp_soln, int i, int b)
{
    int j;
    vp_soln->mapping[i] = b;
    for (j = 0; j < vp_soln->prob->bins[b]->num_dims; j++) {
        vp_soln->loads[b][j] += vp_soln->prob->items[i]->totals[j];
    }
}

int vp_put_item_in_bin_safe(vp_solution_t vp_soln, int i, int b)
{
    if (!vp_item_can_fit_in_bin(vp_soln, i, b)) return 1;
    vp_put_item_in_bin(vp_soln, i, b);
    return 0;
}

// load computed as the sum of all dims of all objects
double vp_compute_sum_load(vp_solution_t vp_soln, int b)
{
  return double_array_sum(vp_soln->loads[b], vp_soln->prob->bins[b]->num_dims);
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
QSORT_CMP_FUNC_DECL(cmp_vector_array_idxs_lex)
{
    vp_vector_t *va = (vp_vector_t *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    int i;

    for (i = 0; i < MIN(va[x]->num_dims, va[y]->num_dims); i++) {
        if (va[x]->totals[i] < va[y]->totals[i]) return -1;
        if (va[x]->totals[i] > va[y]->totals[i]) return 1;
    }
    return CMP(va[x]->num_dims, va[y]->num_dims);
}

QSORT_CMP_FUNC_DECL(rcmp_vector_array_idxs_lex)
{
#ifdef NO_QSORT_R
    return -1 * cmp_vector_array_idxs_lex(xvp, yvp);
#else
    return -1 * cmp_vector_array_idxs_lex(qsort_thunk_vp, xvp, yvp);
#endif
}

/* max comparison of vectors */
// for descinding order sort...
QSORT_CMP_FUNC_DECL(cmp_vector_array_idxs_max)
{
    vp_vector_t *va = (vp_vector_t *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(double_array_max(va[x]->totals, va[x]->num_dims),
        double_array_max(va[y]->totals, va[y]->num_dims));
}

QSORT_CMP_FUNC_DECL(rcmp_vector_array_idxs_max)
{
#ifdef NO_QSORT_R
    return -1 * cmp_vector_array_idxs_max(xvp, yvp);
#else
    return -1 * cmp_vector_array_idxs_max(qsort_thunk_vp, xvp, yvp);
#endif
}


/* Sum comparison of vectors */
QSORT_CMP_FUNC_DECL(cmp_vector_array_idxs_sum)
{
    vp_vector_t *va = (vp_vector_t *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(double_array_sum(va[x]->totals, va[x]->num_dims),
        double_array_sum(va[y]->totals, va[y]->num_dims));
}

QSORT_CMP_FUNC_DECL(rcmp_vector_array_idxs_sum)
{
#ifdef NO_QSORT_R
    return -1 * cmp_vector_array_idxs_sum(xvp, yvp);
#else
    return -1 * cmp_vector_array_idxs_sum(qsort_thunk_vp, xvp, yvp);
#endif
}

QSORT_CMP_FUNC_DECL(cmp_vector_array_idxs_maxratio)
{
    vp_vector_t *va = (vp_vector_t *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(double_array_max(va[x]->totals, va[x]->num_dims) / 
        double_array_min(va[x]->totals, va[x]->num_dims), 
        double_array_max(va[y]->totals, va[y]->num_dims) / 
        double_array_min(va[y]->totals, va[y]->num_dims));
}

QSORT_CMP_FUNC_DECL(rcmp_vector_array_idxs_maxratio)
{
#ifdef NO_QSORT_R
    return -1 * cmp_vector_array_idxs_maxratio(xvp, yvp);
#else
    return -1 * cmp_vector_array_idxs_maxratio(qsort_thunk_vp, xvp, yvp);
#endif
}

QSORT_CMP_FUNC_DECL(cmp_vector_array_idxs_maxdiff)
{
    vp_vector_t *va = (vp_vector_t *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(double_array_max(va[x]->totals, va[x]->num_dims) -
        double_array_min(va[x]->totals, va[x]->num_dims), 
        double_array_max(va[y]->totals, va[y]->num_dims) - 
        double_array_min(va[y]->totals, va[y]->num_dims));
}

QSORT_CMP_FUNC_DECL(rcmp_vector_array_idxs_maxdiff)
{
#ifdef NO_QSORT_R
    return -1 * cmp_vector_array_idxs_maxdiff(xvp, yvp);
#else
    return -1 * cmp_vector_array_idxs_maxdiff(qsort_thunk_vp, xvp, yvp);
#endif
}

int cmp_ints(const void *xvp, const void *yvp)
{
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(x, y);
}

QSORT_CMP_FUNC_DECL(cmp_double_array_idxs)
{
    double *vals = (double *)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    return CMP(vals[x], vals[y]);
}

QSORT_CMP_FUNC_DECL(rcmp_double_array_idxs)
{
#ifdef NO_QSORT_R
    return -1 * cmp_double_array_idxs(xvp, yvp);
#else
    return -1 * cmp_double_array_idxs(qsort_thunk_vp, xvp, yvp);
#endif
}

qsort_cmp_func *get_vp_cmp_func(char *cmp_name) {
    if (NULL == cmp_name) {
        return NULL;
    } else if (!strcmp(cmp_name, "ALEX")) {
        return cmp_vector_array_idxs_lex;
    } else if (!strcmp(cmp_name, "AMAX")) {
        return cmp_vector_array_idxs_max;
    } else if (!strcmp(cmp_name, "ASUM")) {
        return cmp_vector_array_idxs_sum;
    } else if (!strcmp(cmp_name, "AMAXRATIO")) {
        return cmp_vector_array_idxs_maxratio;
    } else if (!strcmp(cmp_name, "AMAXDIFF")) {
        return cmp_vector_array_idxs_maxdiff;
    } else if (!strcmp(cmp_name, "DLEX")) {
        return rcmp_vector_array_idxs_lex;
    } else if (!strcmp(cmp_name, "DMAX")) {
        return rcmp_vector_array_idxs_max;
    } else if (!strcmp(cmp_name, "DSUM")) {
        return rcmp_vector_array_idxs_sum;
    } else if (!strcmp(cmp_name, "DMAXRATIO")) {
        return rcmp_vector_array_idxs_maxratio;
    } else if (!strcmp(cmp_name, "DMAXDIFF")) {
        return rcmp_vector_array_idxs_maxdiff;
    } else {
        fprintf(stderr, "Unknown comparison '%s'\n", cmp_name);
        exit(1);
    }   

    return NULL;
}
