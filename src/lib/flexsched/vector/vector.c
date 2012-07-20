#include "flexsched.h"
#include "vector.h"

#define FIRST_FIT 0
#define BEST_FIT 1

/* Implementation of The standard "Fit" algorithms for
 * solving binpacking:
 *    "FIRST"   or "BEST":  fitting policy
 *    "LEX", "MAX", or "SUM": sorting policy
 *
 * Returns the number of used bins
 */
vp_solution_t solve_vp_problem_FITD(vp_problem_t vp_prob, int args[], 
    qsort_cmp_func *cmp_item_idxs)
{
    int i, j;
    int fit_type = args[0];
    int sortmap[vp_prob->num_items];
    double sumloads[vp_prob->num_bins];
    vp_solution_t vp_soln = new_vp_solution(vp_prob);

    // set up vector sort map
    for (i = 0; i < vp_prob->num_items; i++) sortmap[i] = i;
    if (cmp_item_idxs) {
        QSORT_R(sortmap, vp_prob->num_items, sizeof(int), vp_prob->items, 
            cmp_item_idxs);
    }

    // Place vectors into bins
    switch(fit_type) {
        case FIRST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            for (j = 0; j < vp_prob->num_bins; j++)
                if (!vp_put_item_in_bin_safe(vp_soln, sortmap[i], j)) break;
            if (j >= vp_prob->num_bins) return vp_soln;
        }
        break;
        case BEST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            for (j = 0; j < vp_prob->num_bins; j++) {
                if (vp_item_can_fit_in_bin(vp_soln, sortmap[i], j)) {
                    sumloads[j] = vp_compute_sum_load(vp_soln, j);
                } else {
                    sumloads[j] = -1.0;
                }
            }
            j = double_array_argmax(sumloads, vp_prob->num_bins);
            if (sumloads[j] < 0.0 || 
                vp_put_item_in_bin_safe(vp_soln, sortmap[i], j)) return vp_soln;
        }
        break;
        default:
        fprintf(stderr, "Invalid VP fit type '%d'\n", fit_type);
        exit(1);
    }

    vp_soln->success = 1;
    return vp_soln;
}

// solve vp problem using Permutation Pack or Choose Pack
// FIXME: don't like all the 0 indexes...
vp_solution_t solve_vp_problem_MCB(vp_problem_t vp_prob, int args[], 
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

    vp_solution_t vp_soln = new_vp_solution(vp_prob);

    // initialize dim_perm
    for (j = 0; j < vp_prob->num_dims; j++) {
        dim_perm[j] = j;
    }

    // initialize unmapped vectors and vector dims
    for (i = 0; i < vp_prob->num_items; i++) {
        unmapped_vectors[i] = i;
        QSORT_R(dim_perm, vp_prob->items[i]->num_dims, sizeof(int), 
            vp_prob->items[i]->totals, rcmp_double_array_idxs);
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
            if (vp_item_can_fit_in_bin(vp_soln, v, b)) {
                best_v = v;
                best_v_idx = i;
                break;
            }
        }

        // This probably over-complicates things, but check and make sure that 
        // there are at least 2 vectors before doing all the PP/CP stuff
        for (i++; i < num_unmapped_vectors; i++) {
            v = unmapped_vectors[i];
            if (vp_item_can_fit_in_bin(vp_soln, v, b)) {
                break;
            }
        }

        if (i < num_unmapped_vectors) {

            // compute bin permutation
            QSORT_R(dim_perm, vp_prob->bins[b]->num_dims, sizeof(int), 
                vp_soln->loads[b], cmp_double_array_idxs);

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
            while (j < vp_prob->bins[b]->num_dims) {
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

                if (!vp_item_can_fit_in_bin(vp_soln, v, b)) continue;

                // apply bin key inverse to vector keys to get how they permute 
                // the bin key in the first w elements...
                for (j = 0; j < w; j++) 
                    v_perm[j] = bin_dim_positions[vector_dims[v][j]];

                if (isCP) { // CP ignores position of the first w
                    qsort(v_perm, w, sizeof(int), cmp_ints);
                }

                cmp_val = cmp_int_arrays_lex(w, v_perm, best_perm);

                if (cmp_val < 0 || (0 == cmp_val && 
                    QSORT_CMP_CALL(cmp_item_idxs, vp_prob->items, &v, &best_v) 
                    < 0)) {
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
            vp_put_item_in_bin(vp_soln, best_v, b);
            num_unmapped_vectors--;
            unmapped_vectors[best_v_idx] = 
                unmapped_vectors[num_unmapped_vectors];
        } else {
            b++;
        }

    }

    free(v_perm);
    free(best_perm);

    if (!num_unmapped_vectors) vp_soln->success = 1;

    return vp_soln;
}

vp_solution_t solve_vp_problem_META(vp_problem_t vp_prob, int notargs[],
    qsort_cmp_func cmp_item_idxs)
{
    int args[2] = {0, 0};
    char *sortnames[] = {"ALEX", "AMAX", "ASUM", "AMAXRATIO", "AMAXDIFF", 
                         "DLEX", "DMAX", "DSUM", "DMAXRATIO", "DMAXDIFF",
                         "NONE", NULL };
    char **item_sort_name;
    int isCP, w;
    int i;
    vp_solution_t vp_soln = NULL;

    for (isCP = 0; isCP <= 1; isCP++) {
        args[0] = isCP;
        for (item_sort_name = sortnames; *item_sort_name; item_sort_name++) {
            vp_soln = solve_vp_problem_FITD(vp_prob, args, 
                get_vp_cmp_func(*item_sort_name));
            if (vp_soln && vp_soln->success) {
                sprintf(vp_soln->misc_output, "%s %s", isCP ? "BF" : "FF", 
                    *item_sort_name);
                return vp_soln;
            }
            free_vp_solution(vp_soln);
        }
    }

    for (isCP = 0; isCP <= 1; isCP++) {
        args[0] = isCP;
        for (item_sort_name = sortnames; *item_sort_name; item_sort_name++) {
            for (w = 0; w <= vp_prob->num_dims; w++) {
                args[1] = w;
                vp_soln = solve_vp_problem_MCB(vp_prob, args,
                    get_vp_cmp_func(*item_sort_name));
                if (vp_soln && vp_soln->success) {
                    sprintf(vp_soln->misc_output, "%s %s W%d", 
                        isCP ? "CP" : "PP", *item_sort_name, w);
                    return vp_soln;
                }
                free_vp_solution(vp_soln);
            }
        }
    }

    return new_vp_solution(vp_prob);
}
    
flexsched_solution_t VP_solver(flexsched_problem_t flex_prob, 
    vp_solution_t (*solve_vp_problem)(vp_problem_t, int[], qsort_cmp_func), 
    int args[], qsort_cmp_func cmp_item_idxs) 
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    double yield, yieldlb, yieldub, best_yield;
    vp_problem_t vp_prob = NULL;
    vp_solution_t vp_soln = NULL;
    int i, status;

    yieldlb = 0.0;
    yieldub = 1.0;

    best_yield = -1.0;
    yield = 0.0;

    while(yieldub - yieldlb > 0.001) {
        yield = (yieldub + yieldlb) / 2.0;

        // Generate the VP instance
        vp_prob = new_vp_problem(flex_prob, yield);

        vp_soln = solve_vp_problem(vp_prob, args, cmp_item_idxs);

        // Solve the VP instance
        if (vp_soln->success) {
            yieldlb = yield;
            // Save the computed mapping
            flex_soln->success = 1;
            for (i = 0; i < flex_prob->num_jobs; i++) {
                flex_soln->mapping[i] = vp_soln->mapping[i];
                flex_soln->yields[i] = yield;
            }
            strcpy(flex_soln->misc_output, vp_soln->misc_output);
        } else {
            yieldub = yield;
        }

        free_vp_solution(vp_soln);
        free_vp_problem(vp_prob);

    }

    return flex_soln;
}

/* 
 * VP_scheduler() 
 */
flexsched_solution_t VP_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    vp_solution_t (*solve_vp_problem)(vp_problem_t, int [], qsort_cmp_func) ;
    int args[2] = {0, 0};
    qsort_cmp_func *cmp_item_idxs = NULL;
    char **opt;

    solve_vp_problem = NULL;

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

    if (!strcmp(name, "METAVP")) {
        solve_vp_problem = solve_vp_problem_META;
    }

    return VP_solver(flex_prob, solve_vp_problem, args, cmp_item_idxs);
}

