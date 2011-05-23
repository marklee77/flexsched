#include "flexsched.h"
#include "vector.h"

#define FIRST_FIT 0
#define BEST_FIT 1

int solve_vp_problem_HETERO_FITD(vp_problem vp_prob, int fit_type, 
    int reshuffle_bins, qsort_cmp_func *cmp_item_idxs, 
    qsort_cmp_func *cmp_bin_idxs)
{
    int i, j;
    int v, b;
    int item_sortmap[vp_prob->num_items];
    int bin_sortmap[vp_prob->num_bins];
    float sumloads[vp_prob->num_bins];
    struct vector_array_t va;
    va.num_dims = vp_prob->num_dims;

    // set up vector sort map
    for (i = 0; i < vp_prob->num_items; i++) item_sortmap[i] = i;
    va.vectors = vp_prob->items;
    qsort_r(item_sortmap, vp_prob->num_items, sizeof(int), &va, cmp_item_idxs);

    // Place vectors into bins
    switch(fit_type) {
        case FIRST_FIT:

        // set up bin sort map
        for (j = 0; j < vp_prob->num_bins; j++) bin_sortmap[j] = j;
        va.vectors = vp_prob->capacities;
        qsort_r(bin_sortmap, vp_prob->num_bins, sizeof(int), &va, cmp_bin_idxs);

        for (i = 0; i < vp_prob->num_items; i++) {
            v = item_sortmap[i];
            for (j = 0; j < vp_prob->num_bins; j++) {
                b = bin_sortmap[j];
                if (!vp_put_vector_in_bin_safe(vp_prob, v, b)) break;
            }
            if (j >= vp_prob->num_bins) return 1;
            if (reshuffle_bins) {
                qsort_r(bin_sortmap, vp_prob->num_bins, sizeof(int), &va, 
                    cmp_bin_idxs);
            }
        }
        break;
        case BEST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            v = item_sortmap[i];
            for (j = 0; j < vp_prob->num_bins; j++) {
                if (vp_vector_can_fit_in_bin(vp_prob, v, j)) {
                    sumloads[j] = vp_compute_sum_load(vp_prob, j);
                } else {
                    sumloads[j] = -1.0;
                }
            }
            b = array_argmax(sumloads, vp_prob->num_bins);
            if (sumloads[j] < 0.0 || 
                vp_put_vector_in_bin_safe(vp_prob, v, b)) return 1;
        }
        break;
        default:
        fprintf(stderr,"Invalid VP fit type '%d'\n", fit_type);
        exit(1);
    }

    return 0;
}

// solve vp problem using Permutation Pack or Choose Pack
int solve_vp_problem_HETERO_MCB(vp_problem vp_prob, int w, int isCP, 
    int rescale_items, int reshuffle_bins, qsort_cmp_func *cmp_item_idxs, 
    qsort_cmp_func *cmp_bin_idxs)
{
    int i, j;

    int dims[vp_prob->num_dims];
    int vector_dims[vp_prob->num_items][w];

    int unmapped_vectors[vp_prob->num_items];
    int num_unmapped_vectors;

    int b, v, best_v, best_v_idx, cmp_val;

    float rescaled_item[vp_prob->num_dims];

    int bin_dim_positions[vp_prob->num_dims];

    int *v_perm, *best_perm, *tmp_perm;

    int bin_idxs[vp_prob->num_bins];
    int *open_bins = bin_idxs;
    int num_open_bins;
    float *bin_capacities[vp_prob->num_bins];

    struct vector_array_t va;
    va.num_dims = vp_prob->num_dims;

    // initialize dims
    for (j = 0; j < vp_prob->num_dims; j++) dims[j] = j;

    // initialize unmapped vectors and vector dims
    if (rescale_items) {
        for (i = 0; i < vp_prob->num_items; i++) unmapped_vectors[i] = i;
    } else {
        for (i = 0; i < vp_prob->num_items; i++) {
            unmapped_vectors[i] = i;
            qsort_r(dims, vp_prob->num_dims, sizeof(int), 
                vp_prob->items[i], rcmp_float_array_idxs);
            for (j = 0; j < w; j++) vector_dims[i][j] = dims[j];
        }
    }

    // initialize open bins
    for (i = 0; i < vp_prob->num_bins; i++) {
        open_bins[i] = i;
        bin_capacities[i] = (float *)calloc(vp_prob->num_dims, sizeof(float));
        for (j = 0; j < vp_prob->num_dims; j++) {
            bin_capacities[i][j] = vp_prob->bins[i][j];
        }
    }

    if (cmp_bin_idxs) {
        va.vectors = bin_capacities; 
        qsort_r(open_bins, vp_prob->num_bins, sizeof(int), &va, cmp_bin_idxs);
    }

    // initialize pointers to v_perm and best_perm
    v_perm = (int *)calloc(w, sizeof(int));
    best_perm = (int *)calloc(w, sizeof(int));

    num_open_bins = vp_prob->num_bins;
    num_unmapped_vectors = vp_prob->num_items;

    while (num_open_bins > 0 && num_unmapped_vectors > 0) {

        b = *open_bins;

        best_v = -1;
        best_v_idx = -1;

        // Find the first vector that can be put in the bin
        for (i = 0; i < num_unmapped_vectors; i++) {
            v = unmapped_vectors[i];
            if (vp_vector_can_fit_in_bin(vp_prob, v, b)) {
                best_v = v;
                best_v_idx = i;
                break;
            }
        }

        // This probably over-complicates things, but check and make sure that 
        // there are at least 2 vectors before doing all the PP/CP stuff
        for (i++; i < num_unmapped_vectors; i++) {
            v = unmapped_vectors[i];
            if (vp_vector_can_fit_in_bin(vp_prob, v, b)) {
                break;
            }
        }

        if (i < num_unmapped_vectors) {

            // compute bin permutation
            qsort_r(dims, vp_prob->num_dims, sizeof(int), 
                bin_capacities[b], rcmp_float_array_idxs);

            // compute bin dim positions
            j = 0;
            if (isCP) { // CP treats first w positions as the same...
                while (j < w) {
                    bin_dim_positions[dims[j]] = 0;
                    j++;
                }
            }
            while (j < vp_prob->num_dims) {
                bin_dim_positions[dims[j]] = j;
                j++;
            }

            // apply bin key inverse to best vector key to get how it permutes
            // the bin key in the first w elements...
            // FIXME: should this be done to current or absolute bin capacity?
            if (rescale_items) {
                for (j = 0; j < vp_prob->num_dims; j++) {
                    rescaled_item[j] = 
                        vp_prob->items[best_v][j] / bin_capacities[b][j];
                }
                qsort_r(dims, vp_prob->num_dims, sizeof(int), 
                    rescaled_item, rcmp_float_array_idxs);
                for (j = 0; j < w; j++) 
                    best_perm[j] = bin_dim_positions[dims[j]];
            } else {
                for (j = 0; j < w; j++) 
                    best_perm[j] = bin_dim_positions[vector_dims[best_v][j]];
            }

            if (isCP) { // CP ignores position
                qsort(best_perm, w, sizeof(int), cmp_ints);
            }

            // start with the current vector and look for the "best" one by the
            // CP/PP and secondary criteria...
            for (; i < num_unmapped_vectors; i++) {

                v = unmapped_vectors[i];

                if (!vp_vector_can_fit_in_bin(vp_prob, v, b)) continue;

                // apply bin key inverse to vector keys to get how they permute 
                // the bin key in the first w elements...
                if (rescale_items) {
                    for (j = 0; j < vp_prob->num_dims; j++) {
                        rescaled_item[j] = 
                            vp_prob->items[v][j] / bin_capacities[b][j];
                    }
                    qsort_r(dims, vp_prob->num_dims, sizeof(int), 
                        rescaled_item, rcmp_float_array_idxs);
                    for (j = 0; j < w; j++) 
                        v_perm[j] = bin_dim_positions[dims[j]];
                } else {
                    for (j = 0; j < w; j++) 
                        v_perm[j] = bin_dim_positions[vector_dims[v][j]];
                }

                if (isCP) { // CP ignores position of the first w
                    qsort(v_perm, w, sizeof(int), cmp_ints);
                }

                cmp_val = cmp_int_arrays_lex(w, v_perm, best_perm);
                va.vectors = vp_prob->items;

#ifdef NO_QSORT_R
                global_qsort_vptr = &va;
                if (cmp_val < 0 || (0 == cmp_val && 
                    cmp_item_idxs(&v, &best_v) < 0))
#else
                if (cmp_val < 0 || (0 == cmp_val && 
                    cmp_item_idxs(&va, &v, &best_v) < 0))
#endif
                {
#if 0
                    printf("bin state is: (");
                    for (j = 0; j < vp_prob->num_dims; j++) {
                        printf("%.3f ", vp_prob->loads[b][j]);
                    }   
                    printf ("), (");
                    for (j = 0; j < vp_prob->num_dims; j++) {
                        printf("%d ", bin_dim_positions[j]);
                    }   
                    printf(")\n");
                    printf("cmp is (%d, %d)\n", cmp_val, cmp_item_idxs(&va, &v, &best_v));
                    printf("selecting item %d (", v); 
                    for (j = 0; j < vp_prob->num_dims; j++) {
                        printf("%.3f ", vp_prob->items[v][j]);
                    }   
                    printf("; ");
                    for (j = 0; j < w; j++) {
                        printf("%d ", v_perm[j]);
                    }   
                    printf(") over item %d (", best_v);
                    for (j = 0; j < vp_prob->num_dims; j++) {
                        printf("%.3f ", vp_prob->items[best_v][j]);
                    }   
                    printf("; ");
                    for (j = 0; j < w; j++) {
                        printf("%d ", best_perm[j]);
                    }   
                    printf(")\n");
#endif

                    best_v = v;
                    best_v_idx = i;
                    tmp_perm = best_perm;
                    best_perm = v_perm;
                    v_perm = tmp_perm;
                }

            }

        }

        // if we found a vector put it in the bin and delete from the list of
        // unmapped vectors -- otherwise advance to the next bin
        if (best_v > -1) {
            vp_put_vector_in_bin(vp_prob, best_v, b);
            num_unmapped_vectors--;
            unmapped_vectors[best_v_idx] = 
                unmapped_vectors[num_unmapped_vectors];
            for (j = 0; j < vp_prob->num_dims; j++) {
                bin_capacities[b][j] = 
                    vp_prob->bins[b][j] - vp_prob->loads[b][j];
            }
            if (reshuffle_bins && cmp_bin_idxs) {
                va.vectors = bin_capacities;
                qsort_r(open_bins, num_open_bins, sizeof(int), &va, 
                    cmp_bin_idxs);
            }
        } else {
            open_bins++;
            num_open_bins--;
        }

    }

    for (i = 0; i < vp_prob->num_bins; i++) {
        free(bin_capacities[i]);
    }
    free(v_perm);
    free(best_perm);

    return num_unmapped_vectors;
}
    
int hvp_solver_FF(vp_problem vp_prob, 
    qsort_cmp_func cmp_item_idxs, qsort_cmp_func cmp_bin_idxs)
{
    return solve_vp_problem_HETERO_FITD(vp_prob, FIRST_FIT, 0, cmp_item_idxs,
            cmp_bin_idxs);
}

int hvp_solver_FFR(vp_problem vp_prob, 
    qsort_cmp_func cmp_item_idxs, qsort_cmp_func cmp_bin_idxs)
{
    return solve_vp_problem_HETERO_FITD(vp_prob, FIRST_FIT, 1, cmp_item_idxs,
            cmp_bin_idxs);
}

int hvp_solver_BF(vp_problem vp_prob, 
    qsort_cmp_func cmp_item_idxs, qsort_cmp_func cmp_bin_idxs)
{
    return solve_vp_problem_HETERO_FITD(vp_prob, BEST_FIT, 0, cmp_item_idxs,
            cmp_bin_idxs);
}

// FIXME: work on this...
/* solve_vp_instance() 
 *  Returns the number of bins used
 */
int (*get_hvp_solver(char *vp_algorithm))(vp_problem, qsort_cmp_func,
        qsort_cmp_func)
{
    if (!strcmp(vp_algorithm, "FF")) {
        return hvp_solver_FF;
    } else if (!strcmp(vp_algorithm, "FFR")) {
        return hvp_solver_FFR;
    } else if (!strcmp(vp_algorithm, "BF")) {
        return hvp_solver_BF;
    } else if (!strcmp(vp_algorithm, "CPMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_max, NULL);
    } else if (!strcmp(vp_algorithm, "CPSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_sum, NULL);
    } else if (!strcmp(vp_algorithm, "PPMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_max, NULL);
    } else if (!strcmp(vp_algorithm, "PPSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_sum, NULL);
    } else if (!strcmp(vp_algorithm, "CPMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_max, NULL);
    } else if (!strcmp(vp_algorithm, "CPSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_sum, NULL);
    } else if (!strcmp(vp_algorithm, "PPMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_max, NULL);
    } else if (!strcmp(vp_algorithm, "PPSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_sum, NULL);
    } else if (!strcmp(vp_algorithm, "CPMAXMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPSUMMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPMAXMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPSUMMAX")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPMAXMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPSUMMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPMAXMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPSUMMAXR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPMAXSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPSUMSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPMAXSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPSUMSUM")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPMAXSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPSUMSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPMAXSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPSUMSUMR")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 0,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPMAXMAXS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPSUMMAXS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPMAXMAXS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPSUMMAXS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPMAXMAXRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPSUMMAXRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPMAXMAXRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPSUMMAXRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPMAXSUMS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPSUMSUMS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 0, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPMAXSUMS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPSUMSUMS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 0, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPMAXSUMRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPSUMSUMRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 1, 1, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPMAXSUMRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 1,
            rcmp_vector_array_idxs_max, rcmp_vector_array_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPSUMSUMRS")) {
        retval = solve_vp_problem_HETERO_MCB(vp_prob, 2, 0, 1, 1,
            rcmp_vector_array_idxs_sum, rcmp_vector_array_idxs_sum);
    } else {
        fprintf(stderr, "Unknown vp_algorithm '%s'\n", vp_algorithm);
        exit(1);
    }
    return retval;
}

qsort_cmp_func get_cmp_function(char *cmp_name) {
    if (NULL == cmp_name) {
        return NULL;
    if (!strcmp(cmp_name, "ALEX")) {
        return cmp_vector_array_idxs_lex;
    } else if (!strcmp(cmp_name, "AMAX")) {
        return cmp_vector_array_idxs_max;
    } else if (!strcmp(cmp_name, "ASUM")) {
        return cmp_vector_array_idxs_sum;
    } else if (!strcmp(cmp_name, "AMAXRATIO")) {
        return cmp_vector_array_idxs_maxratio;
    } else if (!strcmp(cmp_name, "AMAXDIFF")) {
        return cmp_vector_array_idxs_maxdiff;
    if (!strcmp(cmp_name, "DLEX")) {
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

/* 
 * VP_scheduler() 
 */
flexsched_solution HVP_scheduler(
    char *vp_algorithm, char *item_cmp_name, char *bin_cmp_name)
{
    flexsched_solution flex_soln = new_flexsched_solution(vp_algorithm);
    double yield, yieldlb, yieldub, best_yield;
    vp_problem vp_prob = NULL;
    int i, status;
    int (*hvp_solver)(vp_prob, qsort_cmp_func, qsort_cmp_func) 
        = get_solver_func(vp_algorithm);
    qsort_cmp_func cmp_item_idxs = get_cmp_func(item_cmp_name);
    qsort_cmp_func cmp_bin_idxs = get_cmp_func(bin_cmp_name);

    yieldlb = 0.0;
    yieldub = compute_LP_bound();

    best_yield = -1.0;
    yield = 0.0;


    while(yieldub - yieldlb > 0.001) {
        yield = (yieldub + yieldlb) / 2.0;

        // Generate the VP instance
        vp_prob = new_vp_problem(yield);

        // Solve the VP instance
        if (hvp_solver(vp_prob, cmp_item_idxs, cmp_bin_idxs)) {
            yieldub = yield;
        } else {
            yieldlb = yield;
            // Save the computed mapping
            flex_soln->success = 1;
            for (i = 0; i < flex_prob->num_services; i++) {
                flex_soln->mapping[i] = vp_prob->mapping[i];
                flex_soln->scaled_yields[i] = yield;
            }
        }

        free_vp_problem(vp_prob);

    }

    return flex_soln;
}
