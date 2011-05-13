#include "flexsched.h"

vp_problem new_vp_problem(float yield) 
{
    int i, j;
    vp_problem vp_prob;

    vp_prob = (vp_problem)calloc(1, sizeof(struct vp_problem_struct));
    vp_prob->num_dims = flex_prob->num_rigid + flex_prob->num_fluid;
    vp_prob->num_vectors = flex_prob->num_services;
    vp_prob->num_bins = flex_prob->num_servers;
    vp_prob->vectors = (float **)calloc(vp_prob->num_vectors, sizeof(float *));
    vp_prob->bin_capacities = 
        (float **)calloc(vp_prob->num_bins, sizeof(float *));
    vp_prob->mapping = (int *)calloc(vp_prob->num_vectors, sizeof(int));
    vp_prob->loads = (float **)calloc(vp_prob->num_bins, sizeof(float *));
    vp_prob->misc = NULL;
    for (i = 0; i< vp_prob->num_vectors; i++) {
        vp_prob->vectors[i] = 
            (float *)calloc(vp_prob->num_dims, sizeof(float));
        vp_prob->mapping[i] = -1;
        for (j = 0; j < flex_prob->num_rigid; j++) {
            vp_prob->vectors[i][j] = flex_prob->rigid_needs[i][j]; 
        }
        for (j = 0; j < flex_prob->num_fluid; j++) {
            vp_prob->vectors[i][flex_prob->num_rigid + j] = 
                (flex_prob->slas[i] + yield * (1.0 - flex_prob->slas[i])) * 
                flex_prob->fluid_needs[i][j];
        }
    }
    for (i = 0; i < vp_prob->num_bins; i++) {
        vp_prob->bin_capacities[i] = 
            (float *)calloc(vp_prob->num_dims, sizeof(float));
        vp_prob->loads[i] = (float *)calloc(vp_prob->num_dims, sizeof(float));
        for (j = 0; j < flex_prob->num_rigid; j++) {
            vp_prob->bin_capacities[i][j] = flex_prob->rigid_capacities[i][j];
            vp_prob->loads[i][j] = 0.0;
        }
        for (j = 0; j < flex_prob->num_fluid; j++) {
            vp_prob->bin_capacities[i][flex_prob->num_rigid + j] = 
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
    for (i = 0; i < vp_prob->num_vectors; i++) {
        free(vp_prob->vectors[i]);
    }
    free(vp_prob->vectors);
    for (i = 0; i < vp_prob->num_bins; i++) {
        free(vp_prob->bin_capacities[i]);
        free(vp_prob->loads[i]);
    }
    free(vp_prob->loads);
    free(vp_prob->misc);
    free(vp_prob->bin_capacities);
    free(vp_prob);
    return;
}

int vp_vector_can_fit_in_bin(vp_problem vp_prob, int v, int b)
{
    int i;
    for (i = 0; i < vp_prob->num_dims; i++) {
        if (vp_prob->loads[b][i] + vp_prob->vectors[v][i] > 
            vp_prob->bin_capacities[b][i])
            return 0;
    }
    return 1;
}

int vp_put_vector_in_bin_safe(vp_problem vp_prob, int v, int b)
{
    int i;
    if (!vp_vector_can_fit_in_bin(vp_prob, v, b)) return 1;
    vp_prob->mapping[v] = b;
    for (i = 0; i < vp_prob->num_dims; i++) {
        vp_prob->loads[b][i] += vp_prob->vectors[v][i];
    }
    return 0;
}

// load computed as the sum of all dims of all objects
float vp_compute_sum_load(vp_problem vp_prob, int b)
{
  return array_sum(vp_prob->loads[b], vp_prob->num_dims);
}

// NOT threadsafe, but apparently we don't have qsort_r or qsort_b on tomate...
vp_problem global_vp_prob;

// compare lexicographically
int rcmp_vp_vector_idxs_lex(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    int i;

    for (i = 0; i < global_vp_prob->num_dims; i++) {
        if (global_vp_prob->vectors[x][i] < global_vp_prob->vectors[y][i]) 
            return 1;
        if (global_vp_prob->vectors[x][i] > global_vp_prob->vectors[y][i]) 
            return -1;
    }
    return 0;
}

/* max comparison of vectors */
// for descinding order sort...
int rcmp_vp_vector_idxs_max(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(array_max(global_vp_prob->vectors[x], global_vp_prob->num_dims),
        array_max(global_vp_prob->vectors[y], global_vp_prob->num_dims));
}

/* Sum comparison of vectors */
// for descinding order sort...
int rcmp_vp_vector_idxs_sum(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(array_sum(global_vp_prob->vectors[x], global_vp_prob->num_dims),
        array_sum(global_vp_prob->vectors[y], global_vp_prob->num_dims));
}

// for descinding order sort...
int rcmp_vp_vector_idxs_maxratio(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(array_max(global_vp_prob->vectors[x], global_vp_prob->num_dims)
        / array_min(global_vp_prob->vectors[x], global_vp_prob->num_dims),
        array_max(global_vp_prob->vectors[y], global_vp_prob->num_dims) 
        / array_min(global_vp_prob->vectors[y], global_vp_prob->num_dims));
}

// for descinding order sort...
int rcmp_vp_vector_idxs_maxdiff(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(array_max(global_vp_prob->vectors[x], global_vp_prob->num_dims) 
        - array_min(global_vp_prob->vectors[x], global_vp_prob->num_dims),
        array_max(global_vp_prob->vectors[y], global_vp_prob->num_dims) 
        - array_min(global_vp_prob->vectors[y], global_vp_prob->num_dims));
}

/* Misc comparison of vectors  */
/* By decreasing order of MISC */
// FIXME: only really used by Maruyama, which is not currently implemented
int rcmp_vp_vector_idxs_misc(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(global_vp_prob->misc[x], global_vp_prob->misc[y]);
}

/* A helper function to compute the thingies from
 * the Maruyama article, and puts them into the misc 
 * field of the vp_prob->vectors
 */
void vp_compute_degrees_of_dominance(vp_problem vp_prob)
{
    int i, j;
    float degrees[vp_prob->num_dims];

    // For each dimension compute its degree of dominance
    for (i = 0; i < vp_prob->num_dims; i++) {
        degrees[i] = 0.0;
        // sets the sums to zero
        for (j = 0; j < vp_prob->num_bins; j++) 
            degrees[i] += ceilf(vp_prob->loads[j][i] - 1.0);

        degrees[i] /= vp_prob->num_bins;
    }

    for (i = 0; i < vp_prob->num_vectors; i++) {
        vp_prob->misc[i] = 0.0;
        for (j = 0; j < vp_prob->num_dims; j++)
            vp_prob->misc[i] += degrees[j] * vp_prob->vectors[i][j];
    }

    return;
}

/* Implementation of The standard "Fit" algorithms for
 * solving binpacking:
 *    "FIRST"   or "BEST":  fitting policy
 *    "LEX", "MAX", or "SUM": sorting policy
 *    "MARUYAMA" sorting: the Maruyama optimization // FIXME: not implemented
 *
 * Returns the number of used bins
 */
int solve_vp_problem_FITD(
    vp_problem vp_prob, const char *fit_type, const char *sort_type)
{
    int i, j, k;
    int (*compar)(void *, const void *, const void *);
    int sortmap[vp_prob->num_vectors];

    // set up vector sort map
    for (i = 0; i < vp_prob->num_vectors; i++) sortmap[i] = i;

    // Sort the vectors in the instance according to the sort type
    global_vp_prob = vp_prob;
    if (!strcmp(sort_type, "LEX")) {
        qsort(sortmap, vp_prob->num_vectors, sizeof(int), 
            rcmp_vp_vector_idxs_lex);
    } else if (!strcmp(sort_type, "MAX")) {
        qsort(sortmap, vp_prob->num_vectors, sizeof(int),
            rcmp_vp_vector_idxs_max);
    } else if (!strcmp(sort_type, "SUM")) {
        qsort(sortmap, vp_prob->num_vectors, sizeof(int),
            rcmp_vp_vector_idxs_sum);
    } else if (strcmp(sort_type, "NONE")) {
        fprintf(stderr,"Invalid VP sort type '%s'\n",sort_type);
        exit(1);
    }

    // Place vectors into bins
    if (!strcmp(fit_type, "FIRST")) { // First Fit
        for (i = 0; i < vp_prob->num_vectors; i++) {
            for (j = 0; j < vp_prob->num_bins; j++)
                if (!vp_put_vector_in_bin_safe(vp_prob, sortmap[i], j)) break;
            if (j >= vp_prob->num_bins) return 1;
        }
    } else if (!strcmp(fit_type, "BEST")) {
        float sumloads[vp_prob->num_bins];
        for (i = 0; i < vp_prob->num_vectors; i++) {
            for (j = 0; j < vp_prob->num_bins; j++) {
                if (vp_vector_can_fit_in_bin(vp_prob, sortmap[i], j)) {
                    sumloads[j] = vp_compute_sum_load(vp_prob, j);
                } else {
                    sumloads[j] = -1.0;
                }
            }
            j = array_argmax(sumloads, vp_prob->num_bins);
            if (sumloads[j] < 0.0 || 
                vp_put_vector_in_bin_safe(vp_prob, sortmap[i], j)) return 1;
        }
    } else {
        fprintf(stderr,"Invalid VP fit type '%s'\n",fit_type);
        exit(1);
    }

    return 0;
}

int cmp_ints(const void *x_ptr, const void *y_ptr) {
    return CMP(*((int *)x_ptr), *((int *)y_ptr));
}

float *global_float_values;
int *global_int_values;

int cmp_float_array_idxs(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return CMP(global_float_values[x], global_float_values[y]);
}

int rcmp_float_array_idxs(const void *x_ptr, const void *y_ptr)
{
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);

    return RCMP(global_float_values[x], global_float_values[y]);
}

// NOT a sort function
int cmp_int_arrays_lex(int a1[], int a2[], int length) {
    int i;
    for (i = 0; i < length; i++) {
        if (a1[i] < a2[i]) return -1;
        if (a1[i] > a2[i]) return 1;
    }
    return 0;
}

int **global_qsort_keys;
int global_length;
int (*global_compar)(const void *, const void *);

int cmp_perm_idxs(const void *x_ptr, const void *y_ptr) {
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);
    return cmp_int_arrays_lex(global_qsort_keys[x], global_qsort_keys[y], 
        global_length);
}

int cmp_vector_idxs_by_perm_then_compar(const void *x_ptr, const void *y_ptr) {
    int x = *((int *)x_ptr);
    int y = *((int *)y_ptr);
    int retval = 
        cmp_int_arrays_lex(global_qsort_keys[x], global_qsort_keys[y], 
        global_length);
    if (!retval) return global_compar(x_ptr, y_ptr);
    return retval;
}


int MCB_put_unmapped_vector_in_bin(vp_problem vp_prob, int **lists, int
        *list_keys[], int *list_sizes, int num_lists, int w, int isCP, int b)
{
    int i, j;
    int v;
    int bin_key[vp_prob->num_dims], bin_key_inverse[vp_prob->num_dims];
    int **list_keys_remapped;
    int lists_sortmap[num_lists];

    // compute bin permutation
    // FIXME: in next version this will be reverse order by avail. capacity...
    for (i = 0; i < vp_prob->num_dims; i++) bin_key[i] = i;
    global_float_values = vp_prob->loads[b];
    qsort(bin_key, vp_prob->num_dims, sizeof(int), cmp_float_array_idxs);

    // compute bin key inverse
    i  = 0;
    if (isCP) { // CP considers the first w elements to have equal weight
        while (i < w) {
            bin_key_inverse[bin_key[i]] = 0;
            i++;
        }
    }
    while (i < vp_prob->num_dims) {
        bin_key_inverse[bin_key[i]] = i;
        i++;
    }

    list_keys_remapped = (int **)calloc(num_lists, sizeof(int *));

    // apply bin key inverse to list keys to get how they permute the bin key in
    // the first w elements...
    for (i = 0; i < num_lists; i++) {
        lists_sortmap[i] = i;
        list_keys_remapped[i] = (int *)calloc(w, sizeof (int));
        for (j = 0; j < w; j++) {
            list_keys_remapped[i][j] = bin_key_inverse[list_keys[i][j]];
        }
        if (isCP) { // CP doesn't care about incomming order
            qsort(list_keys_remapped[i], w, sizeof(int), cmp_ints);
        }
    }

    // sort list keys using sortmap
    global_qsort_keys = list_keys_remapped;
    global_length = w;
    qsort(lists_sortmap, num_lists, sizeof(int), cmp_perm_idxs);

    if (isCP) { // merge lists that are now equal
        int *vectormap = (int *)calloc(vp_prob->num_vectors, sizeof(int));
        int **new_lists = (int **)calloc(vp_prob->num_vectors, sizeof(int *));
        int *new_list_sizes = (int *)calloc(vp_prob->num_vectors, sizeof(int));
        int **new_list_keys = (int **)calloc(vp_prob->num_vectors, sizeof(int));
        int new_num_lists = 0;
        int v = 0, k;
        i = 0;
        while (i < num_lists) {
            new_lists[new_num_lists] = &vectormap[v];
            new_list_keys[new_num_lists] = list_keys[i];
            new_list_sizes[new_num_lists] = list_sizes[i];
            for (k = 0; k < list_sizes[i]; k++) {
                vectormap[v++] = lists[i][k];
            }
            for (j = i+1; j < num_lists && 
                !cmp_int_arrays_lex(list_keys[i], list_keys[j], w); j++) {
                new_list_sizes[new_num_lists] += list_sizes[j];
                for (k = 0; k < list_sizes[j]; k++) {
                    vectormap[v++] = lists[j][k];
                }
            }
            i = j;
            // FIXME: I hate that we need to sort again, and use a global...
            qsort(new_lists[new_num_lists], list_sizes[new_num_lists],
                sizeof(int), global_compar);
            new_num_lists++;
        }
        lists = new_lists;
        num_lists = new_num_lists;
        list_keys = new_list_keys;
        list_sizes = new_list_sizes;
    }

    for (i = 0; i < num_lists; i++) {
        free(list_keys_remapped[i]);
    }
    free(list_keys_remapped);

    // FIXME: CP needs to free memory...

    for (i = 0; i < num_lists; i++) {
        for (j = 0; j < list_sizes[lists_sortmap[i]]; j++) {
            v = lists[lists_sortmap[i]][j];
            if (!vp_put_vector_in_bin_safe(vp_prob, v, b)) {
                // FIXME: does it make sense to remove size 0 lists? 
                list_sizes[lists_sortmap[i]]--;
                for(; j < list_sizes[lists_sortmap[i]]; j++) {
                    lists[lists_sortmap[i]][j] = lists[lists_sortmap[i]][j+1];
                }
                return 0;
            }
        }
    }

    return 1;
}


int solve_vp_problem_MCB(vp_problem vp_prob, int w, char *pack, char *sort_type)
{
    int i, j;

    int vector_sortmap[vp_prob->num_vectors];
    int **lists, num_lists;

    // there's an open question here as to whether it makes more sense to
    // allocate num_vectors*num_dims memory, or to just allocate
    // num_vectors*w+num_dims memory and repeatedly copy the first w elements
    // of the array into the keys...
    int **vector_keys;

    // really the number of vectors is only an upper bound on the number
    // of lists, but it isn't much space and this seems more elegant than
    // repeated calls to malloc.
    int list_sizes[vp_prob->num_vectors], *list_keys[vp_prob->num_vectors];

    int isCP = 0;
    int (*compar)(const void *, const void *);

    if (!strcmp(pack, "PP")) {
        isCP = 0;
    } else if (!strcmp(pack, "CP")) {
        isCP = 1;
    } else {
        fprintf(stderr, "MCB: unknown pack type '%s'\n", pack);
        exit(1);
    }

    // Sort the lists according to the criteria
    // Select the sorting function
    if (!strcmp(sort_type,"MAX")) {
        compar = rcmp_vp_vector_idxs_max;
    } else if (!strcmp(sort_type,"SUM")) {
        compar = rcmp_vp_vector_idxs_sum;
    } else if (!strcmp(sort_type,"MAXDIFF")) {
        compar = rcmp_vp_vector_idxs_maxdiff;
    } else if (!strcmp(sort_type,"MAXRATIO")) {
        compar = rcmp_vp_vector_idxs_maxratio;
    } else {
        fprintf(stderr, "Invalid MCB compare type '%s'\n", sort_type);
        exit(1);
    }

    // initialize vector permutations
    vector_keys = (int **)calloc(vp_prob->num_vectors, sizeof(int *));
    for (i = 0; i < vp_prob->num_vectors; i++) {
        vector_sortmap[i] = i;
        vector_keys[i] = (int *)calloc(vp_prob->num_dims, sizeof(int));
        for (j = 0; j < vp_prob->num_dims; j++) vector_keys[i][j] = j;
        global_float_values = vp_prob->vectors[i];
        qsort(vector_keys[i], vp_prob->num_dims, sizeof(int), 
            rcmp_float_array_idxs);
    }


    // FIXME: theoretically at least there should be a better way than qsort...
    global_qsort_keys = vector_keys;
    global_length = w;
    global_compar = compar;
    global_vp_prob = vp_prob;
    qsort(vector_sortmap, vp_prob->num_vectors, sizeof(int),
            cmp_vector_idxs_by_perm_then_compar);


    // assign vectors to lists based on permutation...
    // FIXME: do realloc here?
    lists = (int **)calloc(vp_prob->num_vectors, sizeof(int *));
    num_lists = 1;
    lists[0] = &vector_sortmap[0];
    list_sizes[0] = 1;
    list_keys[0] = vector_keys[vector_sortmap[0]];
    for (i = 1; i < vp_prob->num_vectors; i++) {
        if (cmp_int_arrays_lex(list_keys[num_lists - 1], 
            vector_keys[vector_sortmap[i]], w)) {
            lists[num_lists] = &vector_sortmap[i];
            list_sizes[num_lists] = 1;
            list_keys[num_lists] = vector_keys[vector_sortmap[i]];
            num_lists++;
        } else {
            list_sizes[num_lists - 1]++;
        }
    }

    // FIXME: in next version we want to sort the bins by some criteria...
    i = 0; // mapped vectors
    j = 0; // current bin
    while (i < vp_prob->num_vectors && j < vp_prob->num_bins) {

        if(MCB_put_unmapped_vector_in_bin(vp_prob, lists, list_keys, 
            list_sizes, num_lists, w, isCP, j)) {
            j++;
        } else {
            i++;
        }

    }

    for (i = 0; i < vp_prob->num_vectors; i++) {
        free(vector_keys[i]);
    }

    free(vector_keys);
    free(lists);

    return (j >= vp_prob->num_bins) ? 1 : 0;
}
    
/* solve_vp_instance() 
 *  Returns the number of bins used
 */
int solve_vp_problem(vp_problem vp_prob, char *vp_algorithm)
{
    int retval;
    if (!strcmp(vp_algorithm, "FFDLEX")) {
        retval = solve_vp_problem_FITD(vp_prob, "FIRST", "LEX");
    } else if (!strcmp(vp_algorithm, "FFDMAX")) {
        retval = solve_vp_problem_FITD(vp_prob, "FIRST", "MAX");
    } else if (!strcmp(vp_algorithm, "FFDSUM")) {
        retval = solve_vp_problem_FITD(vp_prob, "FIRST", "SUM");
    } else if (!strcmp(vp_algorithm, "BFDLEX")) {
        retval = solve_vp_problem_FITD(vp_prob, "BEST", "LEX");
    } else if (!strcmp(vp_algorithm, "BFDMAX")) {
        retval = solve_vp_problem_FITD(vp_prob, "BEST", "MAX");
    } else if (!strcmp(vp_algorithm, "BFDSUM")) {
        retval = solve_vp_problem_FITD(vp_prob, "BEST", "SUM");
    } else if (!strcmp(vp_algorithm, "CPMAX")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "CP", "MAX");
    } else if (!strcmp(vp_algorithm, "CPSUM")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "CP", "SUM");
    } else if (!strcmp(vp_algorithm, "CPMAXRATIO")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "CP", "MAXRATIO");
    } else if (!strcmp(vp_algorithm, "CPMAXDIFF")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "CP", "MAXDIFF");
    } else if (!strcmp(vp_algorithm, "PPMAX")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "PP", "MAX");
    } else if (!strcmp(vp_algorithm, "PPSUM")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "PP", "SUM");
    } else if (!strcmp(vp_algorithm, "PPMAXRATIO")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "PP", "MAXRATIO");
    } else if (!strcmp(vp_algorithm, "PPMAXDIFF")) {
        retval = solve_vp_problem_MCB(vp_prob, 2, "PP", "MAXDIFF");
    } else {
        fprintf(stderr,"Unknown vp_algorithm '%s'\n",vp_algorithm);
        exit(1);
    }
    return retval;
}

/* 
 * VP_scheduler() 
 */
flexsched_solution VP_scheduler(
    char *vp_algorithm, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution(vp_algorithm);
    double yield, yieldlb, yieldub, best_yield;
    vp_problem vp_prob = NULL;
    int i, status;

    yieldlb = 0.0;
    yieldub = compute_LP_bound();

    best_yield = -1.0;
    yield = 0.0;


    while(yieldub - yieldlb > 0.001) {
        yield = (yieldub + yieldlb) / 2.0;

        // Generate the VP instance
        vp_prob = new_vp_problem(yield);

        // Solve the VP instance
        if (solve_vp_problem(vp_prob, vp_algorithm)) {
            yieldub = yield;
        } else {
            yieldlb = yield;
            // Save the computed mapping
            flex_soln->success = 1;
            for (i = 0; i < flex_prob->num_services; i++) {
                flex_soln->mapping[i] = vp_prob->mapping[i];
            }
        }

        free_vp_problem(vp_prob);

    }

    if (flex_soln->success) {
        maximize_minimum_then_average_yield(flex_soln);
    }
  
    return flex_soln;

}
