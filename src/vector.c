#include "flexsched.h"
#include "vector.h"

#define FIRST_FIT 0
#define BEST_FIT 1

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
    vp_problem vp_prob, int fit_type, 
#ifdef NO_QSORT_R
    int (*rcmp_vp_vector_idxs)(const void *, const void *),
#else
    int (*rcmp_vp_vector_idxs)(void *, const void *, const void *)
#endif
    )
{
    int i, j, k;
    int (*compar)(void *, const void *, const void *);
    int sortmap[vp_prob->num_vectors];
    float sumloads[vp_prob->num_bins];

    // set up vector sort map
    for (i = 0; i < vp_prob->num_vectors; i++) sortmap[i] = i;
    qsort_r(sortmap, vp_prob->num_vectors, sizeof(int), vp_prob, 
        rcmp_vp_vector_idxs);

    // Place vectors into bins
    switch(fit_type) {
        case FIRST_FIT:
        for (i = 0; i < vp_prob->num_vectors; i++) {
            for (j = 0; j < vp_prob->num_bins; j++)
                if (!vp_put_vector_in_bin_safe(vp_prob, sortmap[i], j)) break;
            if (j >= vp_prob->num_bins) return 1;
        }
        break;
        case BEST_FIT:
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
        break;
        default:
        fprintf(stderr,"Invalid VP fit type '%s'\n",fit_type);
        exit(1);
    }

    return 0;
}

// solve vp problem using Permutation Pack or Choose Pack
// FIXME: the comparitor doesn't really need to do all this playing around
// with pointers since we don't sort the vectors with qsort anymore...
int solve_vp_problem_MCB(vp_problem vp_prob, int w, int isCP, 
#ifdef NO_QSORT_R
    int (*rcmp_vp_vector_idxs)(const void *, const void *)
#else
    int (*rcmp_vp_vector_idxs)(void *, const void *, const void *)
#endif
    )
{
    int i, j;

    int dim_perm[vp_prob->num_dims];
    int vector_dims[vp_prob->num_vectors][w];

    int unmapped_vectors[vp_prob->num_vectors];
    int num_unmapped_vectors;

    int b, v, best_v, best_v_idx, cmp_val;

    int bin_dim_positions[vp_prob->num_dims];

    int *v_perm, *best_perm, *tmp_perm;

    // initialize dim_perm
    for (j = 0; j < vp_prob->num_dims; j++) {
        dim_perm[j] = j;
    }

    // initialize unmapped vectors and vector dims
    for (i = 0; i < vp_prob->num_vectors; i++) {
        unmapped_vectors[i] = i;
        qsort_r(dim_perm, vp_prob->num_dims, sizeof(int), 
            vp_prob->vectors[i], rcmp_float_array_idxs);
        for (j = 0; j < w; j++) vector_dims[i][j] = dim_perm[j];
    }

    // initialize pointers to v_perm and best_perm
    v_perm = (int *)calloc(w, sizeof(int));
    best_perm = (int *)calloc(w, sizeof(int));

    b = 0;
    num_unmapped_vectors = vp_prob->num_vectors;
    while (b < vp_prob->num_bins && num_unmapped_vectors > 0) {

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
            // FIXME: in next version this will be reverse order by capacity...
            qsort_r(dim_perm, vp_prob->num_dims, sizeof(int), 
                vp_prob->loads[j], cmp_float_array_idxs);

            // compute bin dim positions
            j = 0;
            if (isCP) { // CP treats first w positions as the same...
                while (j < w) {
                    bin_dim_positions[dim_perm[j]] = 0;
                    j++;
                }
            }
            while (j < vp_prob->num_dims) {
                bin_dim_positions[dim_perm[j]] = j;
                j++;
            }

            // apply bin key inverse to best vector key to get how it permutes
            // the bin key in the first w elements...
            for (j = 0; j < w; j++) {
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
                for (j = 0; j < w; j++) 
                    v_perm[j] = bin_dim_positions[vector_dims[v][j]];

                if (isCP) { // CP ignores position of the first w
                    qsort(v_perm, w, sizeof(int), cmp_ints);
                }

                cmp_val = cmp_int_arrays_lex(w, v_perm, best_perm);

#ifdef NO_QSORT_R
                global_qsort_vptr = vp_prob;
                if (cmp_val > 0 || (0 == cmp_val && 
                    rcmp_vp_vector_idxs(&v, &best_v) < 0))
#else
                if (cmp_val > 0 || (0 == cmp_val && 
                    rcmp_vp_vector_idxs(vp_prob, &v, &best_v) < 0))
#endif
                {
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
        } else {
            b++;
        }

    }

    free(v_perm);
    free(best_perm);

    return num_unmapped_vectors;
}
    
/* solve_vp_instance() 
 *  Returns the number of bins used
 */
int solve_vp_problem(vp_problem vp_prob, char *vp_algorithm)
{
    int retval;
    if (!strcmp(vp_algorithm, "FFDLEX")) {
        retval = solve_vp_problem_FITD(vp_prob, FIRST_FIT, 
            rcmp_vp_vector_idxs_lex);
    } else if (!strcmp(vp_algorithm, "FFDMAX")) {
        retval = solve_vp_problem_FITD(vp_prob, FIRST_FIT, 
            rcmp_vp_vector_idxs_max);
    } else if (!strcmp(vp_algorithm, "FFDSUM")) {
        retval = solve_vp_problem_FITD(vp_prob, FIRST_FIT, 
            rcmp_vp_vector_idxs_sum);
    } else if (!strcmp(vp_algorithm, "BFDLEX")) {
        retval = solve_vp_problem_FITD(vp_prob, BEST_FIT, 
            rcmp_vp_vector_idxs_lex);
    } else if (!strcmp(vp_algorithm, "BFDMAX")) {
        retval = solve_vp_problem_FITD(vp_prob, BEST_FIT, 
            rcmp_vp_vector_idxs_max);
    } else if (!strcmp(vp_algorithm, "BFDSUM")) {
        retval = solve_vp_problem_FITD(vp_prob, BEST_FIT, 
            rcmp_vp_vector_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPMAX")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 1, rcmp_vp_vector_idxs_max);
    } else if (!strcmp(vp_algorithm, "CPSUM")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 1, rcmp_vp_vector_idxs_sum);
    } else if (!strcmp(vp_algorithm, "CPMAXRATIO")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 1, 
            rcmp_vp_vector_idxs_maxratio);
    } else if (!strcmp(vp_algorithm, "CPMAXDIFF")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 1, 
            rcmp_vp_vector_idxs_maxdiff);
    } else if (!strcmp(vp_algorithm, "PPMAX")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 0, rcmp_vp_vector_idxs_max);
    } else if (!strcmp(vp_algorithm, "PPSUM")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 0, rcmp_vp_vector_idxs_sum);
    } else if (!strcmp(vp_algorithm, "PPMAXRATIO")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 0,
            rcmp_vp_vector_idxs_maxratio);
    } else if (!strcmp(vp_algorithm, "PPMAXDIFF")) {
        retval = solve_vp_problem_MCB(vp_prob, 1, 0, 
            rcmp_vp_vector_idxs_maxdiff);
    } else {
        fprintf(stderr, "Unknown vp_algorithm '%s'\n", vp_algorithm);
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
                flex_soln->scaled_yields[i] = yield;
            }
        }

        free_vp_problem(vp_prob);

    }

    return flex_soln;
}
