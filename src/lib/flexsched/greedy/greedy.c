#include "flexsched.h"

// FIXME: check on ordering...
// FIXME: this is actually kind of hard to do in a way..., do we need two sets
// of greedy algs the same way we have 2 sets of vector packing algorithms?

/* Keeping the jobs in their original order */
QSORT_CMP_FUNC_DECL(cmp_greedy_S1)
{
    return 0;
}

/* Sorting the jobs by decreasing of max fluid needs */
QSORT_CMP_FUNC_DECL(cmp_greedy_S2)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    // y is first to make it DECREASING
    return CMP(double_array_max(flex_prob->jobs[y]->total_fluid_needs, 
        flex_prob->num_resources),
        double_array_max(flex_prob->jobs[x]->total_fluid_needs,
            flex_prob->num_resources));
}

/* Sorting the jobs by decreasing sum of fluid needs */
QSORT_CMP_FUNC_DECL(cmp_greedy_S3)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    // y is first to make it DECREASING
    return CMP(double_array_sum(flex_prob->jobs[y]->total_fluid_needs,
        flex_prob->num_resources),
        double_array_sum(flex_prob->jobs[x]->total_fluid_needs,
            flex_prob->num_resources));
}

/* Sorting the jobs by decreasing max of rigid requirements */
QSORT_CMP_FUNC_DECL(cmp_greedy_S4)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    // y is first to make it DECREASING
    return CMP(
        double_array_max(flex_prob->jobs[y]->total_rigid_requirements,
            flex_prob->num_resources),
        double_array_max(flex_prob->jobs[x]->total_rigid_requirements,
            flex_prob->num_resources));
}

/* Sorting the jobs by decreasing sum of rigid requirements */
QSORT_CMP_FUNC_DECL(cmp_greedy_S5)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    // y is first to make it DECREASING
    return CMP(
        double_array_sum(flex_prob->jobs[y]->total_rigid_requirements,
            flex_prob->num_resources),
        double_array_sum(flex_prob->jobs[x]->total_rigid_requirements,
            flex_prob->num_resources));
}

/* Sorting by decreasing max of sum rigid requirements and fluid needs*/
QSORT_CMP_FUNC_DECL(cmp_greedy_S6)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;
    double xload[flex_prob->num_resources];
    double yload[flex_prob->num_resources];
    int i;

    for (i = 0; i < flex_prob->num_resources; i++) {
        xload[i] = flex_prob->jobs[x]->total_rigid_requirements[i] +
            flex_prob->jobs[x]->total_fluid_needs[i];
        yload[i] = flex_prob->jobs[y]->total_rigid_requirements[i] +
            flex_prob->jobs[y]->total_fluid_needs[i];
    }

    return CMP(double_array_max(yload, flex_prob->num_resources), 
        double_array_max(xload, flex_prob->num_resources));
}

/* Sorting the jobs by decreasing sum of all requirements and needs */
QSORT_CMP_FUNC_DECL(cmp_greedy_S7)
{
    flexsched_problem_t flex_prob = (flexsched_problem_t)qsort_thunk_vp;
    int x = *(int *)xvp, y = *(int *)yvp;

    return CMP(double_array_sum(flex_prob->jobs[y]->total_rigid_requirements,
        flex_prob->num_resources) + 
        double_array_sum(flex_prob->jobs[y]->total_fluid_needs, 
            flex_prob->num_resources),
        double_array_sum(flex_prob->jobs[x]->total_rigid_requirements, 
            flex_prob->num_resources) + 
        double_array_sum(flex_prob->jobs[x]->total_fluid_needs, 
            flex_prob->num_resources));
}

void GREEDY_sort_jobs(int *sortmap, flexsched_problem_t flex_prob, 
    qsort_cmp_func cmp_items)
{
    int i;

    // Initialize index array with the original order
    for (i = 0; i < flex_prob->num_jobs; i++) sortmap[i] = i;

    // Sort the jobs
    QSORT_R(sortmap, flex_prob->num_jobs, sizeof(int), flex_prob, 
        cmp_items);

    return;
}

// most available resource in max fluid need dimension
int GREEDY_pick_node_P1(flexsched_solution_t flex_soln, int job)
{
    int i;
    int dim = double_array_argmax(
        flex_soln->prob->jobs[job]->total_fluid_needs, 
        flex_soln->prob->num_resources);
    double available_resource, max_available_resource = 0.0;
    int picked = -1;
    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        if (!job_can_fit_on_node_fast(flex_soln, job, i)) continue;
        available_resource = compute_available_resource_fast(flex_soln, i, dim);
        if(available_resource > max_available_resource) {
            max_available_resource = available_resource;
            picked = i;
        }
    }

    return picked;
}

// new: minimum (sumload / sumresource) *after* placement
int GREEDY_pick_node_P2(flexsched_solution_t flex_soln, int job)
{
    int i, j;
    double sumload, sumresource;
    double ratio, minratio = -1.0;
    int picked;
  
    picked = -1;

    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        if (!job_can_fit_on_node_fast(flex_soln, job, i)) continue;
        sumload = 0.0;
        sumresource = 0.0;
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            sumload += compute_fluid_load_fast(flex_soln, i, j) +
                flex_soln->prob->jobs[job]->total_fluid_needs[j];
            sumresource += compute_available_resource_fast(flex_soln, i, j) -
                flex_soln->prob->jobs[job]->total_rigid_requirements[j];
        }
        if (sumresource != 0.0) {
            ratio = sumload / sumresource;
            if (picked < 0 || ratio < minratio) {
                minratio = ratio;
                picked = i;
            }
        }
    }

    if (picked > -1 && minratio < 0.0) {
        fprintf(stderr, "Error: something went wrong in P2.\n");
        //exit(1);
    }

    return picked;
}

// consider largest rigid resource requirement dimension and pick a node based
// on either best fit or worst fit.
int GREEDY_pick_node_P3_P5(
    flexsched_solution_t flex_soln, int job, int mode)
{
    int i, picked;
    int dim = double_array_argmax(
        flex_soln->prob->jobs[job]->total_rigid_requirements, 
            flex_soln->prob->num_resources);
    double available_resource, obj_available_resource = -1.0;

    picked = -1;
    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        if (!job_can_fit_on_node_fast(flex_soln, job, i)) continue;

        available_resource = compute_available_resource_fast(flex_soln, i, dim);

        if (picked < 0 || 
            (mode == BEST_FIT && available_resource < obj_available_resource) ||
            (mode == WORST_FIT && available_resource > obj_available_resource)) 
        {
            picked = i;
            obj_available_resource = available_resource;
        }

    }

    return picked;
} 

int GREEDY_pick_node_P3(flexsched_solution_t flex_soln, int job)
{
  return GREEDY_pick_node_P3_P5(flex_soln, job, BEST_FIT);
}

int GREEDY_pick_node_P5(flexsched_solution_t flex_soln, int job)
{
  return GREEDY_pick_node_P3_P5(flex_soln, job, WORST_FIT);
}

// smallest or biggest total available resource
int GREEDY_pick_node_P4_P6(
    flexsched_solution_t flex_soln, int job, int mode)
{
    int i, j, picked;
    double sumresource, obj_sumresource = -1.0;

    picked = -1;

    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        if (!job_can_fit_on_node_fast(flex_soln, job, i)) continue;
        sumresource = 0.0;
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            sumresource += compute_available_resource_fast(flex_soln, i, j);
        }

        if (picked < 0 || 
            (mode == BEST_FIT && 
                sumresource < obj_sumresource) ||
            (mode == WORST_FIT && 
                sumresource > obj_sumresource)) 
        {
            picked = i;
            obj_sumresource = sumresource;
        }
    }

    return picked;
}

int GREEDY_pick_node_P4(flexsched_solution_t flex_soln, int job)
{
  return GREEDY_pick_node_P4_P6(flex_soln, job, BEST_FIT);
}

int GREEDY_pick_node_P6(flexsched_solution_t flex_soln, int job)
{
  return GREEDY_pick_node_P4_P6(flex_soln, job, WORST_FIT);
}

// first one it can fit on
int GREEDY_pick_node_P7(flexsched_solution_t flex_soln, int job)
{
    int i, picked;

    picked = -1;
    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        if (!job_can_fit_on_node_fast(flex_soln, job, i)) continue;
        picked = i;
        break;
    }

    return picked;
}

void GREEDY_compute_mapping(flexsched_solution_t flex_soln, 
    int (*pick_node)(flexsched_solution_t, int), int *sortmap)
{
    int i, j; 
    int job, node;

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        job = sortmap[i];
        node = pick_node(flex_soln, job);

        if (-1 == node) return;
        
        put_job_on_node_fast(flex_soln, job, node);

    }

    flex_soln->success = 1;

    return;
}

flexsched_solution_t GREEDY_solver(flexsched_problem_t flex_prob, 
    qsort_cmp_func *cmp_items, int (*pick_node)(flexsched_solution_t, int))
{
    int sortmap[flex_prob->num_jobs];
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);

    if (NULL == cmp_items || NULL == pick_node) {
        fprintf(stderr, "Invalid GREEDY sorting or picking function.\n");
        return flex_soln;
    }

    // needed for _fast
    initialize_global_resource_availabilities_and_loads(flex_prob);

    // Sort the jobs appropriately
    GREEDY_sort_jobs(sortmap, flex_prob, cmp_items);

    // Map each job to a particular host
    GREEDY_compute_mapping(flex_soln, pick_node, sortmap);

    // greedy has to do this since it doesn't have a concept of allocation
    if (flex_soln->success) {
        maximize_minimum_yield(flex_soln);
    }

    // cleanup
    free_global_resource_availabilities_and_loads(flex_prob);

    return flex_soln;
}

flexsched_solution_t GREEDY_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options) 
{
    qsort_cmp_func *cmp_items = NULL;
    int (*pick_node)(flexsched_solution_t, int) = NULL;
    char **opt;

    for (opt = options; *opt; opt++) {
        if (!strcmp(*opt, "S1")) {
            cmp_items = cmp_greedy_S1;
        } else if (!strcmp(*opt, "S2")) {
            cmp_items = cmp_greedy_S2;
        } else if (!strcmp(*opt, "S3")) {
            cmp_items = cmp_greedy_S3;
        } else if (!strcmp(*opt, "S4")) {
            cmp_items = cmp_greedy_S4;
        } else if (!strcmp(*opt, "S5")) {
            cmp_items = cmp_greedy_S5;
        } else if (!strcmp(*opt, "S6")) {
            cmp_items = cmp_greedy_S6;
        } else if (!strcmp(*opt, "S7")) {
            cmp_items = cmp_greedy_S7;
        } else if (!strcmp(*opt, "P1")) {
            pick_node = GREEDY_pick_node_P1;
        } else if (!strcmp(*opt, "P2")) {
            pick_node = GREEDY_pick_node_P2;
        } else if (!strcmp(*opt, "P3")) {
            pick_node = GREEDY_pick_node_P3;
        } else if (!strcmp(*opt, "P4")) {
            pick_node = GREEDY_pick_node_P4;
        } else if (!strcmp(*opt, "P5")) {
            pick_node = GREEDY_pick_node_P5;
        } else if (!strcmp(*opt, "P6")) {
            pick_node = GREEDY_pick_node_P6;
        } else if (!strcmp(*opt, "P7")) {
            pick_node = GREEDY_pick_node_P7;
        }
    }

    return GREEDY_solver(flex_prob, cmp_items, pick_node);
}

flexsched_solution_t METAGREEDY_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);

    qsort_cmp_func *sorting[] = {cmp_greedy_S1, cmp_greedy_S2, 
        cmp_greedy_S3, cmp_greedy_S4, cmp_greedy_S5, cmp_greedy_S6, 
        cmp_greedy_S7, NULL};
    int (*picking[])(flexsched_solution_t, int) = {GREEDY_pick_node_P1, 
        GREEDY_pick_node_P2, GREEDY_pick_node_P3, GREEDY_pick_node_P4, 
        GREEDY_pick_node_P5, GREEDY_pick_node_P6, GREEDY_pick_node_P7, 
        NULL};
    int i, is, ip;

    flexsched_solution_t curr_soln = NULL;
    double minyield, maxminyield = -1;
    double avgyield, maxavgyield = -1;

    for (is = 0; sorting[is]; is++) {
        for (ip = 0; picking[ip]; ip++) {

            // call the GREEDY scheduler
            curr_soln = GREEDY_solver(flex_prob, sorting[is], picking[ip]);

            // If worked and improved, then great
            if (curr_soln->success) {
                flex_soln->success = 1;

                // compute minimum yield
                minyield = compute_minimum_yield(curr_soln);

                // compute maximum yield
                avgyield = compute_average_yield(curr_soln);

                if (minyield > maxminyield || (minyield > maxminyield - EPSILON
                    && avgyield > maxavgyield)) {
                    maxminyield = minyield;
                    maxavgyield = avgyield;
                    // save the returned answer
                    for (i=0; i < flex_prob->num_jobs; i++) {
                        flex_soln->mapping[i] = curr_soln->mapping[i];
                        flex_soln->yields[i] = curr_soln->yields[i];
                    }
                    sprintf(flex_soln->misc_output, "S%d P%d", is+1, ip+1);
                }
            }
            free_flexsched_solution(curr_soln);
        }
    }
    return flex_soln;
}

flexsched_solution_t METAGREEDYLIGHT_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);

    qsort_cmp_func *sorting[] = {cmp_greedy_S3, cmp_greedy_S2, cmp_greedy_S5, 
        cmp_greedy_S7, cmp_greedy_S6, cmp_greedy_S2, cmp_greedy_S3, 
        cmp_greedy_S6, cmp_greedy_S5, NULL};
    int (*picking[])(flexsched_solution_t, int) = {GREEDY_pick_node_P1, 
        GREEDY_pick_node_P1, GREEDY_pick_node_P1, GREEDY_pick_node_P1, 
        GREEDY_pick_node_P1, GREEDY_pick_node_P2, GREEDY_pick_node_P2, 
        GREEDY_pick_node_P2, GREEDY_pick_node_P4, NULL};
    int i, isp;

    flexsched_solution_t curr_soln = NULL;
    double minyield, maxminyield = -1;
    double avgyield, maxavgyield = -1;

    for (isp = 0; sorting[isp] && picking[isp]; isp++) {

        // call the GREEDY scheduler
        curr_soln = GREEDY_solver(flex_prob, sorting[isp], picking[isp]);

        // If worked and improved, then great
        if (curr_soln->success) {
            flex_soln->success = 1;

            // compute minimum yield
            minyield = compute_minimum_yield(curr_soln);

            // compute maximum yield
            avgyield = compute_average_yield(curr_soln);

            if (minyield > maxminyield || 
                (minyield > maxminyield - EPSILON && avgyield > maxavgyield)) {
                maxminyield = minyield;
                maxavgyield = avgyield;
                // save the returned answer
                for (i=0; i < flex_prob->num_jobs; i++) {
                    flex_soln->mapping[i] = curr_soln->mapping[i];
                    flex_soln->yields[i] = curr_soln->yields[i];
                }
            }
        }

        free_flexsched_solution(curr_soln);
    }

    return flex_soln;
}
