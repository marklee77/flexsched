#include "flexsched.h"
#ifdef NO_QSORT_R
void *qsort_thunk_vp = NULL;
#endif

double double_array_sum(double *array, int size)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < size; i++) sum += array[i];
    return sum;
}

/* array_max(): Utility function */
double double_array_max(double *array, int size)
{
    int i;
    double max = array[0];
    for (i = 1; i < size; i++) if (max < array[i]) max = array[i];
    return max;
}

/* array_min(): Utility function */
double double_array_min(double *array, int size)
{
    int i;
    double min = array[0];
    for (i = 1; i < size; i++) if (min > array[i]) min = array[i];
    return min;
}

/* array_argmax(): Utility function */
int double_array_argmax(double *array, int size)
{
    int i;
    int argmax = 0;
    double max = array[0];
    for (i = 1; i < size; i++) 
        if (max < array[i]) {
            argmax = i;
            max = array[i];
        }
    return argmax;
}

/* array_argmax(): Utility function */
int double_array_argmin(double *array, int size)
{
    int i;
    int argmin = 0;
    double min = array[0];
    for (i = 1; i < size; i++) 
        if (array[i] < min) {
            argmin = i;
            min = array[i];
        }
    return argmin;
}

void delete_int_array_index(int *array, int idx, int size) {
    int i;
    for (i = idx; i < size - 1; i++) {
        array[i] = array[i+1];
    }
    return;
}

flexsched_solution_t new_flexsched_solution(flexsched_problem_t flex_prob) 
{
    int i;
    flexsched_solution_t flex_soln = 
        (flexsched_solution_t)malloc(sizeof(struct flexsched_solution_s));

    if (!flex_soln) {
        fprintf(stderr, 
                "couldn't allocate sufficient memory for new solution!\n");
        exit(1);
    }

    flex_soln->prob = flex_prob;
    flex_soln->success = 0;
    flex_soln->misc_output[0] = '\0';

    if (!(flex_soln->mapping = 
        (int *)calloc(flex_prob->num_jobs, sizeof(int)))) {
        fprintf(stderr, 
                "couldn't allocate sufficient memory for new mapping!\n");
        exit(1);
    }
    if (!(flex_soln->yields =
        (double *)calloc(flex_prob->num_jobs, sizeof(double)))) {
        fprintf(stderr, "couldn't allocate sufficient memory for yields!\n");
        exit(1);
    }
    for (i = 0; i < flex_prob->num_jobs; i++) {
        flex_soln->mapping[i] = -1;
        flex_soln->yields[i] = 0.0;
    }
    return flex_soln;
}

void free_flexsched_solution(flexsched_solution_t flex_soln)
{
    if (!flex_soln) return;
    free(flex_soln->mapping);
    free(flex_soln->yields);
    free(flex_soln);
    return;
}

inline int unit_requirements_within_capacity_in_dim(
    flexsched_problem_t flex_prob, int job, int node, int dim)
{
    return flex_prob->jobs[job]->unit_rigid_requirements[dim]
        <= flex_prob->nodes[node]->unit_capacities[dim] - EPSILON;
}

int unit_requirements_within_capacity(
    flexsched_problem_t flex_prob, int job, int node) 
{
    int i;
    for (i = 0; i < flex_prob->num_resources; i++) {
        if (!unit_requirements_within_capacity_in_dim(flex_prob, job, 
            node, i)) return 0;
    }
    return 1;
}

inline int unit_allocation_at_yield_within_capacity_in_dim(flexsched_problem_t 
    flex_prob, int job, double yield, int node, int dim) 
{
    return flex_prob->jobs[job]->unit_rigid_requirements[dim] +
        yield * flex_prob->jobs[job]->unit_fluid_needs[dim]
        < flex_prob->nodes[node]->unit_capacities[dim] + EPSILON;
}

// FIXME: painfully hacky
inline int unit_allocation_at_yield_within_capacity_in_dim2(flexsched_problem_t 
    flex_prob, int job, double yield, int node, int dim) 
{
    return flex_prob->jobs[job]->unit_rigid_requirements[dim] +
        yield * flex_prob->jobs[job]->actual_unit_fluid_needs[dim]
        < flex_prob->nodes[node]->unit_capacities[dim] + EPSILON;
}

int unit_allocation_at_yield_within_capacity(
    flexsched_problem_t flex_prob, int job, double yield, int node) 
{
    int i;
    for (i = 0; i < flex_prob->num_resources; i++) {
        if (!unit_allocation_at_yield_within_capacity_in_dim(flex_prob, job,
            yield, node, i)) return 0;
    }
    return 1;
}

// not really necessary...
void put_job_on_node(
    flexsched_solution_t flex_soln, int job, int node) 
{
    flex_soln->mapping[job] = node;
    return;
}

double compute_allocated_resource(
    flexsched_solution_t flex_soln, int node, int dim)
{
    double allocated_resource = EPSILON;
    int i;
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        if (flex_soln->mapping[i] != node) continue;
        allocated_resource +=
            flex_soln->prob->jobs[i]->total_rigid_requirements[dim] +
            flex_soln->yields[i]*
                flex_soln->prob->jobs[i]->total_fluid_needs[dim];
    }
    return allocated_resource;
}

double compute_available_resource(flexsched_solution_t flex_soln, int node, 
    int dim)
{
    double allocated_resource = EPSILON; 
    int i;

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        if (flex_soln->mapping[i] != node) continue;
        allocated_resource += 
            flex_soln->prob->jobs[i]->total_rigid_requirements[dim];
    }

    return MAX(0.0, flex_soln->prob->nodes[node]->total_capacities[dim] - 
        allocated_resource); 
}

double compute_fluid_load(flexsched_solution_t flex_soln, int node, int dim)
{
    double fluid_load = 0.0;
    int i;

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        if (flex_soln->mapping[i] != node) continue;
        fluid_load += flex_soln->prob->jobs[i]->total_fluid_needs[dim];
    }

    return fluid_load; 
}

int total_requirements_within_available_capacity_in_dim(
    flexsched_solution_t flex_soln, int job, int node, int dim)
{
    return flex_soln->prob->jobs[job]->total_rigid_requirements[dim]
        <= compute_available_resource(flex_soln, node, dim) - EPSILON;
}

int job_can_fit_on_node(
    flexsched_solution_t flex_soln, int job, int node)
{
    int i;

    for (i = 0; i < flex_soln->prob->num_resources; i++)
        if (!unit_requirements_within_capacity_in_dim(flex_soln->prob, job, 
            node, i) || 
            !total_requirements_within_available_capacity_in_dim(flex_soln, 
                job, node, i)) return 0;

    return 1;
}

double **global_available_resources;
double **global_fluid_loads;

void initialize_global_resource_availabilities_and_loads(
    flexsched_problem_t flex_prob) 
{
    int i, j;

    global_available_resources = 
        (double **)calloc(flex_prob->num_nodes, sizeof(double *));
    global_fluid_loads = 
        (double **)calloc(flex_prob->num_nodes, sizeof(double *));

    for (i = 0; i < flex_prob->num_nodes; i++) {
        global_available_resources[i] = 
            (double *)calloc(flex_prob->num_resources, sizeof(double));
        global_fluid_loads[i] = 
            (double *)calloc(flex_prob->num_resources, sizeof(double));
        for (j = 0; j < flex_prob->num_resources; j++) {
            global_available_resources[i][j] =
                flex_prob->nodes[i]->total_capacities[j];
            global_fluid_loads[i][j] = 0.0;
        }
    }
}

void free_global_resource_availabilities_and_loads(
    flexsched_problem_t flex_prob) 
{
    int i;

    for (i = 0; i < flex_prob->num_nodes; i++) {
        free(global_available_resources[i]);
        free(global_fluid_loads[i]);
    }
    free(global_available_resources);
    free(global_fluid_loads);
    global_available_resources = NULL;
    global_fluid_loads = NULL;

    return;
}

void put_job_on_node_fast(
    flexsched_solution_t flex_soln, int job, int node) 
{
    int i;
    
    flex_soln->mapping[job] = node;

    for (i = 0; i < flex_soln->prob->num_resources; i++) {
        global_available_resources[node][i] -=
            flex_soln->prob->jobs[job]->total_rigid_requirements[i];
        global_fluid_loads[node][i] +=
            flex_soln->prob->jobs[job]->total_fluid_needs[i];
    }

    return;
}

inline double compute_available_resource_fast(
    flexsched_solution_t flex_soln, int node, int dim)
{
    return global_available_resources[node][dim];
}

inline double compute_fluid_load_fast(
    flexsched_solution_t flex_soln, int node, int dim)
{
    return global_fluid_loads[node][dim];
}

inline int total_requirements_within_available_capacity_in_dim_fast(
    flexsched_solution_t flex_soln, int job, int node, int dim)
{
    return flex_soln->prob->jobs[job]->total_rigid_requirements[dim]
        <= global_available_resources[node][dim] - EPSILON;
}

int job_can_fit_on_node_fast(
    flexsched_solution_t flex_soln, int job, int node)
{
    int i;

    for (i = 0; i < flex_soln->prob->num_resources; i++)
        if (!unit_requirements_within_capacity_in_dim(flex_soln->prob, job, 
            node, i) || 
            !total_requirements_within_available_capacity_in_dim_fast(flex_soln,
                job, node, i)) return 0;

    return 1;
}

void maximize_minimum_yield_on_node(
    flexsched_solution_t flex_soln, int node)
{
    int i, j;
    double load;
    double minyield = 1.0;

    // yield limited by total needs
    for (i = 0; i < flex_soln->prob->num_resources; i++) {
        load = compute_fluid_load(flex_soln, node, i);
        if (load > 0.0) { 
            minyield = MIN(minyield, 
                compute_available_resource(flex_soln, node, i) / load);
        }
    }

    // yield limited by unit needs
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        if (flex_soln->mapping[i] != node) continue;
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            if (flex_soln->prob->jobs[i]->unit_fluid_needs[j] > 0.0) {
                minyield = MIN(minyield, 
                    (flex_soln->prob->nodes[node]->unit_capacities[j] -
                      flex_soln->prob->jobs[i]->unit_rigid_requirements[j])
                    / flex_soln->prob->jobs[i]->unit_fluid_needs[j]);
            }
        }
    }

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        if (flex_soln->mapping[i] != node) continue;
        flex_soln->yields[i] = MAX(EPSILON, minyield - EPSILON);
    }

    return;
}

double compute_minimum_yield(flexsched_solution_t flex_soln)
{
    return double_array_min(flex_soln->yields, flex_soln->prob->num_jobs);
}

void maximize_minimum_yield(flexsched_solution_t flex_soln)
{
    int i;
    double minyield;

    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        maximize_minimum_yield_on_node(flex_soln, i);
    }
    minyield = compute_minimum_yield(flex_soln);
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        flex_soln->yields[i] = minyield;
    }

    return;
}

// FIXME: organizationally I don't know if it makes more sense to put this in
// or linearprog.c
/*
void maximize_average_yield_given_minimum(flexsched_solution_t, double);

void maximize_minimum_then_average_yield(flexsched_solution_t flex_soln)
{
    double minyield;

    maximize_minimum_yield(flex_soln);
    minyield = compute_minimum_yield(flex_soln);
    maximize_average_yield_given_minimum(flex_soln, minyield);
    return;
}
*/

double compute_average_yield(flexsched_solution_t flex_soln)
{
    int i;
    double sumyield = 0.0;
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        sumyield += flex_soln->yields[i];
    }
    return (sumyield / flex_soln->prob->num_jobs);
}

double compute_utilization(flexsched_solution_t flex_soln)
{
    int i;
    double total_capacity = 0.0;
    double total_alloc = 0.0;

    for (i = 0; i < flex_soln->prob->num_nodes; i++) {
        total_capacity +=
            double_array_sum(flex_soln->prob->nodes[i]->total_capacities,
                flex_soln->prob->num_resources);
    }

    if (total_capacity <= 0.0) return 1.0;

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        total_alloc +=
            double_array_sum(
                flex_soln->prob->jobs[i]->total_rigid_requirements, 
                flex_soln->prob->num_resources) + 
            flex_soln->yields[i] * 
              double_array_sum(flex_soln->prob->jobs[i]->total_fluid_needs,
                    flex_soln->prob->num_resources);
    }

    return (total_alloc / total_capacity);
}


/* sanity_check() 
 *
 * Checks that resource capacities are not overcome and that
 * the mapping is valid
 */
int sanity_check(flexsched_problem_t flex_prob, int mapping[], double yields[])
{
    int i, j, k;
    int retval = 0;
    double allocated_resources[flex_prob->num_nodes][flex_prob->num_resources];

    // check that each job is mapped to a node and yields are <= 1.0
    for (i = 0; i < flex_prob->num_jobs; i++) {
        if (mapping[i] < 0 || 
            mapping[i] >= flex_prob->num_nodes) {
            fprintf(stderr, 
                "Error: Service %d is mapped to invalid node %d\n.", i, 
                mapping[i]);
            retval = 1;
        }
        if (yields[i] < EPSILON || yields[i] > 1.0) {
            fprintf(stderr, "Error: Allocation of job %d is %.2f.\n",
                i, yields[i]);
            retval = 1;
        }
    }

    if (retval) return retval;

    // initialize allocated resources
    for (i = 0; i < flex_prob->num_nodes; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            allocated_resources[i][j] = 0.0;
        }
    }

    // check unit allocations and sum up allocations
    for (i = 0; i < flex_prob->num_jobs; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            if(!unit_allocation_at_yield_within_capacity_in_dim(flex_prob,
                i, yields[i], mapping[i], j))
            {
                fprintf(stderr, 
                    "Error: Allocation of job %d to node %d with yield %.3f exceeds unit capacity in resource dimension %d (%.3f + %.3f = %.3f/%.3f).\n", 
                    i, mapping[i], yields[i], j,
                    flex_prob->jobs[i]->unit_rigid_requirements[j],
                    yields[i] * 
                        flex_prob->jobs[i]->unit_fluid_needs[j], 
                    flex_prob->jobs[i]->unit_rigid_requirements[j] + 
                        yields[i] * 
                        flex_prob->jobs[i]->unit_fluid_needs[j], 
                    flex_prob->nodes[mapping[i]]->unit_capacities[j]);
                retval = 1;
            }
            allocated_resources[mapping[i]][j] += 
                flex_prob->jobs[i]->total_rigid_requirements[j] + 
                yields[i] * 
                    flex_prob->jobs[i]->total_fluid_needs[j];
        }
    }

    for (i = 0; i < flex_prob->num_nodes; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            if (flex_prob->nodes[i]->total_capacities[j] + EPSILON <= 
                allocated_resources[i][j]) {
                fprintf(stderr, 
                    "Error: Total resource allocation on node %d exceeds capacity in dimension %d (%.6f/%.6f).\n",
                    i, j, allocated_resources[i][j], 
                    flex_prob->nodes[i]->total_capacities[j]);
                retval = 1;
            }
        }
    }

    return retval;
}

int sanity_check2(flexsched_problem_t flex_prob, int mapping[], double yields[])
{
    int i, j, k;
    int retval = 0;
    double allocated_resources[flex_prob->num_nodes][flex_prob->num_resources];

    // check that each job is mapped to a node and yields are <= 1.0
    for (i = 0; i < flex_prob->num_jobs; i++) {
        if (mapping[i] < 0 || 
            mapping[i] >= flex_prob->num_nodes) {
            fprintf(stderr, 
                "Error: Service %d is mapped to invalid node %d\n.", i, 
                mapping[i]);
            retval = 1;
        }
        if (yields[i] < EPSILON || yields[i] > 1.0) {
            fprintf(stderr, "Error: Allocation of job %d is %.2f.\n",
                i, yields[i]);
            retval = 1;
        }
    }

    if (retval) return retval;

    // initialize allocated resources
    for (i = 0; i < flex_prob->num_nodes; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            allocated_resources[i][j] = 0.0;
        }
    }

    // check unit allocations and sum up allocations
    for (i = 0; i < flex_prob->num_jobs; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            if(!unit_allocation_at_yield_within_capacity_in_dim2(flex_prob,
                i, yields[i], mapping[i], j))
            {
                fprintf(stderr, 
                    "Error: Allocation of job %d to node %d with yield %.3f exceeds unit capacity in resource dimension %d (%.3f + %.3f = %.3f/%.3f).\n", 
                    i, mapping[i], yields[i], j,
                    flex_prob->jobs[i]->unit_rigid_requirements[j],
                    yields[i] * 
                        flex_prob->jobs[i]->actual_unit_fluid_needs[j], 
                    flex_prob->jobs[i]->unit_rigid_requirements[j] + 
                        yields[i] * 
                        flex_prob->jobs[i]->actual_unit_fluid_needs[j], 
                    flex_prob->nodes[mapping[i]]->unit_capacities[j]);
                retval = 1;
            }
            allocated_resources[mapping[i]][j] += 
                flex_prob->jobs[i]->total_rigid_requirements[j] + 
                yields[i] * 
                    flex_prob->jobs[i]->actual_total_fluid_needs[j];
        }
    }

    for (i = 0; i < flex_prob->num_nodes; i++) {
        for (j = 0; j < flex_prob->num_resources; j++) {
            if (flex_prob->nodes[i]->total_capacities[j] + EPSILON <= 
                allocated_resources[i][j]) {
                fprintf(stderr, 
                    "Error: Total resource allocation on node %d exceeds capacity in dimension %d (%.6f/%.6f).\n",
                    i, j, allocated_resources[i][j], 
                    flex_prob->nodes[i]->total_capacities[j]);
                retval = 1;
            }
        }
    }

    return retval;
}
