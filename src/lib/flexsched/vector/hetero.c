#include "flexsched.h"
#include "vector.h"

#define FIRST_FIT 0
#define BEST_FIT 1

vp_solution_t solve_hvp_problem_FITD(vp_problem_t vp_prob, int args[], 
    qsort_cmp_func *cmp_item_idxs, qsort_cmp_func *cmp_bin_idxs)
{
    int fit_type = args[0];
    int i, j, k;
    int v, b;
    int item_sortmap[vp_prob->num_items];
    int bin_sortmap[vp_prob->num_bins];
    double sumcapacities[vp_prob->num_bins];
    vp_solution_t vp_soln = new_vp_solution(vp_prob);

    // set up vector sort map
    for (i = 0; i < vp_prob->num_items; i++) item_sortmap[i] = i;
    if (cmp_item_idxs)
        QSORT_R(item_sortmap, vp_prob->num_items, sizeof(int), vp_prob->items, 
            cmp_item_idxs);

    // set up bin sort map
    for (j = 0; j < vp_prob->num_bins; j++) bin_sortmap[j] = j;
    if (cmp_bin_idxs)
        QSORT_R(bin_sortmap, vp_prob->num_bins, sizeof(int), vp_prob->bins, 
            cmp_bin_idxs);

    // Place vectors into bins
    switch(fit_type) {
        case FIRST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            v = item_sortmap[i];
            for (j = 0; j < vp_prob->num_bins; j++) {
                b = bin_sortmap[j];
                if (!vp_put_item_in_bin_safe(vp_soln, v, b)) break;
            }
            if (j >= vp_prob->num_bins) return vp_soln; 
        }
        break;
        case BEST_FIT:
        for (i = 0; i < vp_prob->num_items; i++) {
            v = item_sortmap[i];
            for (j = 0; j < vp_prob->num_bins; j++) {
                b = bin_sortmap[j];
                if (vp_item_can_fit_in_bin(vp_soln, v, b)) {
                    sumcapacities[b] = 0.0;
                    for (k = 0; k < vp_prob->num_dims; k++) {
                        sumcapacities[b] += 
                            vp_prob->bins[b]->totals[k] - vp_soln->loads[b][k];
                    }
                } else {
                    sumcapacities[b] = 2.0 * vp_prob->num_dims + 1.0;
                }
            }
            b = double_array_argmin(sumcapacities, vp_prob->num_bins);
            if (sumcapacities[j] > 2.0 * vp_prob->num_dims || 
                vp_put_item_in_bin_safe(vp_soln, v, b)) return vp_soln;
        }
        break;
        default:
        fprintf(stderr,"Invalid VP fit type '%d'\n", fit_type);
        exit(1);
    }
    vp_soln->success = 1;
    return vp_soln;
}

// solve vp problem using Permutation Pack or Choose Pack
vp_solution_t solve_hvp_problem_MCB(vp_problem_t vp_prob, int args[],
    qsort_cmp_func *cmp_item_idxs, qsort_cmp_func *cmp_bin_idxs)
{
    int i, j;
    int isCP = args[0];
    int w = MIN(args[1], vp_prob->num_dims);

    int dims[vp_prob->num_dims];
    int vector_dims[vp_prob->num_items][w];

    int unmapped_vectors[vp_prob->num_items];
    int num_unmapped_vectors = vp_prob->num_items;

    int b, v, best_v, best_v_idx, cmp_val;

    int bin_dim_ranks[vp_prob->num_dims];

    int *v_perm, *best_perm, *tmp_perm;

    int bin_idxs[vp_prob->num_bins];
    int *open_bins = bin_idxs;
    int num_open_bins = vp_prob->num_bins;

    vp_solution_t vp_soln = new_vp_solution(vp_prob);

    // initialize dims
    for (j = 0; j < vp_prob->num_dims; j++) dims[j] = j;

    // initialize unmapped vectors and vector dims
    for (i = 0; i < vp_prob->num_items; i++) {
        unmapped_vectors[i] = i;
        QSORT_R(dims, vp_prob->num_dims, sizeof(int), 
            vp_prob->items[i]->totals, rcmp_double_array_idxs);
        for (j = 0; j < w; j++) vector_dims[i][j] = dims[j];
    }

    // initialize open bins
    for (i = 0; i < vp_prob->num_bins; i++) open_bins[i] = i;

    if (cmp_bin_idxs)
        QSORT_R(open_bins, num_open_bins, sizeof(int), vp_prob->bins, 
            cmp_bin_idxs);

    // initialize pointers to v_perm and best_perm
    v_perm = (int *)calloc(w, sizeof(int));
    best_perm = (int *)calloc(w, sizeof(int));

    while (num_open_bins > 0 && num_unmapped_vectors > 0) {

        b = *open_bins;

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
            QSORT_R(dims, vp_prob->num_dims, sizeof(int),
                vp_soln->capacities[b], rcmp_double_array_idxs);

            // compute bin dim positions
            bin_dim_ranks[dims[0]] = 0; // make sure j > 0
            j = 1;
            // CP treats first w positions as the same...
            if (isCP) for (; j < w; j++) bin_dim_ranks[dims[j]] = 0;

            // we assume j > 0 from above...
            for (; j < vp_prob->num_dims; j++) {
                if (vp_soln->capacities[b][dims[j-1]] == 
                    vp_soln->capacities[b][dims[j]]) {
                    bin_dim_ranks[dims[j]] = bin_dim_ranks[dims[j-1]];
                } else {
                    bin_dim_ranks[dims[j]] = bin_dim_ranks[dims[j-1]]+1;
                }
            }

            // apply bin key inverse to best vector key to get how it permutes
            // the bin key in the first w elements...
            for (j = 0; j < w; j++) 
                best_perm[j] = bin_dim_ranks[vector_dims[best_v][j]];

            // CP ignores position
            if (isCP) qsort(best_perm, w, sizeof(int), cmp_ints);

            // start with the current vector and look for the "best" one by the
            // CP/PP and secondary criteria...
            for (; i < num_unmapped_vectors; i++) {

                v = unmapped_vectors[i];

                if (!vp_item_can_fit_in_bin(vp_soln, v, b)) continue;

                // apply bin key inverse to vector keys to get how they permute 
                // the bin key in the first w elements...
                for (j = 0; j < w; j++) 
                    v_perm[j] = bin_dim_ranks[vector_dims[v][j]];

                // CP ignores position of the first w
                if (isCP) qsort(v_perm, w, sizeof(int), cmp_ints);

                cmp_val = cmp_int_arrays_lex(w, v_perm, best_perm);

                if (cmp_val < 0 || (0 == cmp_val && cmp_item_idxs &&
                    QSORT_CMP_CALL(cmp_item_idxs, vp_prob->items, &v, &best_v) 
                    < 0)) 
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
            vp_put_item_in_bin(vp_soln, best_v, b);
            num_unmapped_vectors--;
            unmapped_vectors[best_v_idx] = 
                unmapped_vectors[num_unmapped_vectors];
        } else {
            open_bins++;
            num_open_bins--;
        }

    }

    free(v_perm);
    free(best_perm);

    if (!num_unmapped_vectors) vp_soln->success = 1;

    return vp_soln;
}

double pairwise_distance(vp_solution_t vp_soln, int b, int v) {
    int i, j;
    double x;
    double retval = 0;
    for (i = 0; i < vp_soln->prob->num_dims - 1; i++) {
        for (j = i+1; j < vp_soln->prob->num_dims; j++) {
            x = (vp_soln->capacities[b][i] + 
                    vp_soln->prob->items[v]->totals[i]) - 
                    (vp_soln->capacities[b][j] + 
                        vp_soln->prob->items[v]->totals[j]);
            retval += x*x;
        }
    }
    //printf("%d %d distance: %g\n", i, j, retval);
    return retval;
}

// solve vp problem using Permutation Pack or Choose Pack
vp_solution_t solve_hvp_problem_MinCD(vp_problem_t vp_prob, int args[],
    qsort_cmp_func *cmp_item_idxs, qsort_cmp_func *cmp_bin_idxs)
{
    int i, j;
    int isCP = args[0];

    int unmapped_vectors[vp_prob->num_items];
    int num_unmapped_vectors = vp_prob->num_items;

    int b, v, best_v, best_v_idx;
    double best_val, val;

    int bin_idxs[vp_prob->num_bins];
    int *open_bins = bin_idxs;
    int num_open_bins = vp_prob->num_bins;

    vp_solution_t vp_soln = new_vp_solution(vp_prob);

    // initialize unmapped vectors and vector dims
    for (i = 0; i < vp_prob->num_items; i++) unmapped_vectors[i] = i;

    // initialize open bins
    for (i = 0; i < vp_prob->num_bins; i++) open_bins[i] = i;

    if (cmp_bin_idxs)
        QSORT_R(open_bins, num_open_bins, sizeof(int), vp_prob->bins, 
            cmp_bin_idxs);

    while (num_open_bins > 0 && num_unmapped_vectors > 0) {

        b = *open_bins;

        best_v = -1;
        best_v_idx = -1;

        // Find the first vector that can be put in the bin
        for (i = 0; i < num_unmapped_vectors; i++) {
            v = unmapped_vectors[i];
            if (vp_item_can_fit_in_bin(vp_soln, v, b)) {
                best_v = v;
                best_v_idx = i;
                best_val = pairwise_distance(vp_soln, b, v);
                break;
            }
        }

        for (i++; i < num_unmapped_vectors; i++) {
            v = unmapped_vectors[i];
            if (!vp_item_can_fit_in_bin(vp_soln, v, b)) continue;

            val = pairwise_distance(vp_soln, b, v);

            // FIXME: pick the one that minimizes pairwise distance...
            if (val < best_val || (val == best_val && cmp_item_idxs &&
                QSORT_CMP_CALL(cmp_item_idxs, vp_prob->items, &v, &best_v) < 0))
            {
                best_v = v;
                best_v_idx = i;
                best_val = val;
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
            open_bins++;
            num_open_bins--;
        }

    }

    if (!num_unmapped_vectors) vp_soln->success = 1;

    return vp_soln;
}
    
vp_solution_t solve_hvp_problem_META(vp_problem_t vp_prob, int notargs[],
    qsort_cmp_func cmp_item_idxs, qsort_cmp_func cmp_bin_idxs)
{

    int args[2] = {0, 0};
    char *sortnames[] = { "ALEX", "AMAX", "ASUM", "AMAXRATIO", "AMAXDIFF",
        "DLEX", "DMAX", "DSUM", "DMAXRATIO", "DMAXDIFF", "NONE", NULL };
    char **item_sort_name, **bin_sort_name;
    int isCP, w;
    int i;
    vp_solution_t vp_soln = NULL;

    for (isCP = 0; isCP <= 1; isCP++) {
        args[0] = isCP;
        for (item_sort_name = sortnames; *item_sort_name; item_sort_name++) {
            for (bin_sort_name = sortnames; *bin_sort_name; bin_sort_name++) {
                vp_soln = solve_hvp_problem_FITD(vp_prob, args, 
                    get_vp_cmp_func(*item_sort_name), 
                    get_vp_cmp_func(*bin_sort_name));
                if (vp_soln && vp_soln->success) {
                    sprintf(vp_soln->misc_output, "%s %s %s", 
                        isCP ? "BF" : "FF", *item_sort_name, *bin_sort_name);
                    return vp_soln;
                }
                free_vp_solution(vp_soln);
            }
        }
    }

    for (isCP = 0; isCP <= 1; isCP++) {
        args[0] = isCP;
        for (item_sort_name = sortnames; *item_sort_name; item_sort_name++) {
            for (bin_sort_name = sortnames; *bin_sort_name; bin_sort_name++) {
                for (w = 0; w <= vp_prob->num_dims; w++) {
                    args[1] = w;
                    vp_soln = solve_hvp_problem_MCB(vp_prob, args, 
                        get_vp_cmp_func(*item_sort_name), 
                        get_vp_cmp_func(*bin_sort_name)); 
                    if (vp_soln && vp_soln->success) {
                        sprintf(vp_soln->misc_output, 
                            "%s %s %s W%d", 
                            isCP ? "CP" : "PP", 
                            *item_sort_name, *bin_sort_name, w);
                        return vp_soln;
                    }
                    free_vp_solution(vp_soln);
                }
            }
        }
    }

    return new_vp_solution(vp_prob);
}

vp_solution_t solve_hvp_problem_METALIGHT(vp_problem_t vp_prob, int notargs[],
    qsort_cmp_func cmp_item_idxs, qsort_cmp_func cmp_bin_idxs)
{

    int args[2] = {0, 0};
    char *itemsorts[] = { "DMAX", "DSUM", NULL };
    char *binsorts[] = { "AMAX", "ASUM", NULL };
    char **item_sort_name, **bin_sort_name;
    int isCP, w;
    int i;
    vp_solution_t vp_soln = NULL;

    args[0] = FIRST_FIT;
    for (item_sort_name = itemsorts; *item_sort_name; item_sort_name++) {
        vp_soln = solve_hvp_problem_FITD(vp_prob, args, 
            get_vp_cmp_func(*item_sort_name), NULL);
        if (vp_soln && vp_soln->success) { 
            sprintf(vp_soln->misc_output, "FF %s NONE", *item_sort_name);
            return vp_soln;
        }
        free_vp_solution(vp_soln);
    }

    // don't really know about CP vs PP for only 2 dims...
    args[0] = 0;
    for (item_sort_name = itemsorts; *item_sort_name; item_sort_name++) {
        for (bin_sort_name = binsorts; *bin_sort_name; bin_sort_name++) {
            for (w = 1; w < vp_prob->num_dims; w++) {
                args[1] = w;
                vp_soln = solve_hvp_problem_MCB(vp_prob, args, 
                    get_vp_cmp_func(*item_sort_name), 
                    get_vp_cmp_func(*bin_sort_name)); 
                if (vp_soln && vp_soln->success) {
                    sprintf(vp_soln->misc_output, "PP %s %s W%d", 
                        *item_sort_name, *bin_sort_name, w);
                    return vp_soln;
                }
                free_vp_solution(vp_soln);
            }
        }
    }

    return new_vp_solution(vp_prob);
}

flexsched_solution_t HVP_solver(flexsched_problem_t flex_prob,
    vp_solution_t (*solve_hvp_problem)(vp_problem_t, int[], qsort_cmp_func, 
        qsort_cmp_func), int args[], qsort_cmp_func cmp_item_idxs, 
    qsort_cmp_func cmp_bin_idxs) 
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    double yield, yieldlb, yieldub;
    vp_problem_t vp_prob = NULL;
    vp_solution_t vp_soln = NULL;
    int i, status;

    yieldlb = 0.0;
    yieldub = 1.0;

    yield = 0.0;

    while(yieldub - yieldlb > 0.001) {
        yield = (yieldub + yieldlb) / 2.0;

        // Generate the VP instance
        vp_prob = new_vp_problem(flex_prob, yield);

        // Solve the VP instance
        vp_soln = solve_hvp_problem(vp_prob, args, cmp_item_idxs, cmp_bin_idxs);
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
 * HVP_scheduler() 
 */
flexsched_solution_t HVP_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    vp_solution_t (*solve_hvp_problem)(vp_problem_t, int [], 
        qsort_cmp_func, qsort_cmp_func)
        = NULL;
    int args[2] = {0, 0};
    qsort_cmp_func *cmp_tmp = NULL;
    qsort_cmp_func *cmp_item_idxs = NULL;
    qsort_cmp_func *cmp_bin_idxs = NULL;
    char **opt;
    int item_sort_specified = 0;

    // could be a little more rigorous here...
    for (opt = options; *opt; opt++) {
        if (!strcmp(*opt, "FF")) {
            solve_hvp_problem = solve_hvp_problem_FITD;
            args[0] = FIRST_FIT;
        } else if (!strcmp(*opt, "BF")) {
            solve_hvp_problem = solve_hvp_problem_FITD;
            args[0] = BEST_FIT;
        } else if (!strcmp(*opt, "PP")) {
            solve_hvp_problem = solve_hvp_problem_MCB;
            args[0] = 0;
            args[1] = 1; // default W
        } else if (!strcmp(*opt, "MinCD")) {
            solve_hvp_problem = solve_hvp_problem_MinCD;
            args[0] = 0;
            args[1] = 1; // default W
        } else if (!strcmp(*opt, "CP")) {
            solve_hvp_problem = solve_hvp_problem_MCB;
            args[0] = 1;
            args[1] = 1;
        } else if ('W' == **opt) {
            args[1] = atoi(*opt + 1);
        } else {
            cmp_tmp = get_vp_cmp_func(*opt);
            if (!item_sort_specified) {
                cmp_item_idxs = cmp_tmp;
                item_sort_specified = 1;
            } else {
                cmp_bin_idxs = cmp_tmp;
            }
        }
    }

    if (!strcmp(name, "METAHVP")) {
        solve_hvp_problem = solve_hvp_problem_META;
    }

    if (!strcmp(name, "METAHVPLIGHT")) {
        solve_hvp_problem = solve_hvp_problem_METALIGHT;
    }

    return HVP_solver(flex_prob, solve_hvp_problem, args, cmp_item_idxs, 
        cmp_bin_idxs);

}
