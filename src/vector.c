#include "flexsched.h"
#include "vector.h"

#define FIRST_FIT 0
#define BEST_FIT 1

/* A helper function to compute the thingies from
 * the Maruyama article, and puts them into the misc 
 * field of the vp_prob->items
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

    for (i = 0; i < vp_prob->num_items; i++) {
        vp_prob->misc[i] = 0.0;
        for (j = 0; j < vp_prob->num_dims; j++)
            vp_prob->misc[i] += degrees[j] * vp_prob->items[i][j];
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
int solve_vp_problem_FITD(vp_problem vp_prob, int args[], 
    qsort_cmp_func *cmp_item_idxs)
{
    int i, j;
    int fit_type = args[0];
    int sortmap[vp_prob->num_items];
    float sumloads[vp_prob->num_bins];
    struct vector_array_t va;
    va.num_dims = vp_prob->num_dims;

    // set up vector sort map
    for (i = 0; i < vp_prob->num_items; i++) sortmap[i] = i;
    va.vectors = vp_prob->items;
    qsort_r(sortmap, vp_prob->num_items, sizeof(int), &va, 
        cmp_item_idxs);

    // Place vectors into bins
    switch(fit_type) {
        case FIRST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            for (j = 0; j < vp_prob->num_bins; j++)
                if (!vp_put_vector_in_bin_safe(vp_prob, sortmap[i], j)) break;
            if (j >= vp_prob->num_bins) return 1;
        }
        break;
        case BEST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
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
        fprintf(stderr, "Invalid VP fit type '%d'\n", fit_type);
        exit(1);
    }

    return 0;
}

// solve vp problem using Permutation Pack or Choose Pack
// FIXME: the comparitor doesn't really need to do all this playing around
// with pointers since we don't sort the vectors with qsort anymore...
int solve_vp_problem_MCB(vp_problem vp_prob, int args[], 
    qsort_cmp_func *cmp_item_idxs)
{
    int i, j;
    int isCP = args[0];
    int w = MIN(args[1], vp_prob->num_dims);

    int dim_perm[vp_prob->num_dims];
    int vector_dims[vp_prob->num_items][w];

    int unmapped_vectors[vp_prob->num_items];
    int num_unmapped_vectors;

    int b, v, best_v, best_v_idx, cmp_val;

    int bin_dim_positions[vp_prob->num_dims];

    int *v_perm, *best_perm, *tmp_perm;

    struct vector_array_t va;

    va.num_dims = vp_prob->num_dims;

    // initialize dim_perm
    for (j = 0; j < vp_prob->num_dims; j++) {
        dim_perm[j] = j;
    }

    // initialize unmapped vectors and vector dims
    for (i = 0; i < vp_prob->num_items; i++) {
        unmapped_vectors[i] = i;
        qsort_r(dim_perm, vp_prob->num_dims, sizeof(int), 
            vp_prob->items[i], rcmp_float_array_idxs);
        for (j = 0; j < w; j++) vector_dims[i][j] = dim_perm[j];
    }

    // initialize pointers to v_perm and best_perm
    v_perm = (int *)calloc(w, sizeof(int));
    best_perm = (int *)calloc(w, sizeof(int));

    b = 0;
    num_unmapped_vectors = vp_prob->num_items;
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
            qsort_r(dim_perm, vp_prob->num_dims, sizeof(int), 
                vp_prob->loads[b], cmp_float_array_idxs);

#if 0
            printf("sort of dims: (");
            for (j = 0; j < vp_prob->num_dims; j++) {
                printf("%.3f ", vp_prob->loads[b][j]);
            }
            printf("), (");
            for (j = 0; j < vp_prob->num_dims; j++) {
                printf("%d ", dim_perm[j]);
            }
            printf(")\n");
#endif

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
        } else {
            b++;
        }

    }

    free(v_perm);
    free(best_perm);

    return num_unmapped_vectors;
}
    
void VP_solver(flexsched_solution flex_soln, 
    int (*solve_vp_problem)(vp_problem, int[], qsort_cmp_func), 
    int args[], qsort_cmp_func cmp_item_idxs) 
{
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
        if (solve_vp_problem(vp_prob, args, cmp_item_idxs)) {
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
}

/* 
 * VP_scheduler() 
 */
flexsched_solution VP_scheduler(char *name, char **options)
{
    int (*solve_vp_problem)(vp_problem, int [], qsort_cmp_func) = NULL;
    int args[2];
    qsort_cmp_func *cmp_item_idxs = NULL;
    flexsched_solution flex_soln = new_flexsched_solution();
    char **opt;

    // could be a little more rigorous here...
    for (opt = options; *opt; opt++) {
        if (!strcmp(*opt, "FF")) {
            solve_vp_problem = solve_vp_problem_FITD;
            args[0] = FIRST_FIT;
        } else if (!strcmp(*opt, "BF")) {
            solve_vp_problem = solve_vp_problem_FITD;
            args[0] = BEST_FIT;
        } else if (!strcmp(*opt, "PP")) {
            solve_vp_problem = solve_vp_problem_MCB;
            args[0] = 0;
            args[1] = 1; // default W
        } else if (!strcmp(*opt, "CP")) {
            solve_vp_problem = solve_vp_problem_MCB;
            args[0] = 1;
            args[1] = 1;
        } else if ('W' == **opt) {
            args[1] = atoi(*opt + 1);
        } else {
            cmp_item_idxs = get_vp_cmp_func(*opt);
        }
    }

    VP_solver(flex_soln, solve_vp_problem, args, cmp_item_idxs);
    
    return flex_soln;
}
