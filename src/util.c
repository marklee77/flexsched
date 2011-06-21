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
        (int *)calloc(flex_prob->num_services, sizeof(int)))) {
        fprintf(stderr, 
                "couldn't allocate sufficient memory for new mapping!\n");
        exit(1);
    }
    if (!(flex_soln->yields =
        (double *)calloc(flex_prob->num_services, sizeof(double)))) {
        fprintf(stderr, "couldn't allocate sufficient memory for yields!\n");
        exit(1);
    }
    for (i = 0; i < flex_prob->num_services; i++) {
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
    flexsched_problem_t flex_prob, int service, int server, int dim)
{
    return flex_prob->services[service]->unit_rigid_requirements[dim]
        <= flex_prob->servers[server]->unit_capacities[dim] - EPSILON;
}

int unit_requirements_within_capacity(
    flexsched_problem_t flex_prob, int service, int server) 
{
    int i;
    for (i = 0; i < flex_prob->num_resources; i++) {
        if (!unit_requirements_within_capacity_in_dim(flex_prob, service, 
            server, i)) return 0;
    }
    return 1;
}

inline int unit_allocation_at_yield_within_capacity_in_dim(flexsched_problem_t 
    flex_prob, int service, double yield, int server, int dim) 
{
    return flex_prob->services[service]->unit_rigid_requirements[dim] +
        yield * flex_prob->services[service]->unit_fluid_needs[dim]
        < flex_prob->servers[server]->unit_capacities[dim] + EPSILON;
}

int unit_allocation_at_yield_within_capacity(
    flexsched_problem_t flex_prob, int service, double yield, int server) 
{
    int i;
    for (i = 0; i < flex_prob->num_resources; i++) {
        if (!unit_allocation_at_yield_within_capacity_in_dim(flex_prob, service,
            yield, server, i)) return 0;
    }
    return 1;
}

// not really necessary...
void put_service_on_server(
    flexsched_solution_t flex_soln, int service, int server) 
{
    flex_soln->mapping[service] = server;
    return;
}

double compute_allocated_resource(
    flexsched_solution_t flex_soln, int server, int dim)
{
    double allocated_resource = EPSILON;
    int i;
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        allocated_resource +=
            flex_soln->prob->services[i]->total_rigid_requirements[dim] +
            flex_soln->yields[i]*
                flex_soln->prob->services[i]->total_fluid_needs[dim];
    }
    return allocated_resource;
}

double compute_available_resource(flexsched_solution_t flex_soln, int server, 
    int dim)
{
    double allocated_resource = EPSILON; 
    int i;

    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        allocated_resource += 
            flex_soln->prob->services[i]->total_rigid_requirements[dim];
    }

    return MAX(0.0, flex_soln->prob->servers[server]->total_capacities[dim] - 
        allocated_resource); 
}

double compute_fluid_load(flexsched_solution_t flex_soln, int server, int dim)
{
    double fluid_load = 0.0;
    int i;

    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        fluid_load += flex_soln->prob->services[i]->total_fluid_needs[dim];
    }

    return fluid_load; 
}

int total_requirements_within_available_capacity_in_dim(
    flexsched_solution_t flex_soln, int service, int server, int dim)
{
    return flex_soln->prob->services[service]->total_rigid_requirements[dim]
        <= compute_available_resource(flex_soln, server, dim) - EPSILON;
}

int service_can_fit_on_server(
    flexsched_solution_t flex_soln, int service, int server)
{
    int i;

    for (i = 0; i < flex_soln->prob->num_resources; i++)
        if (!unit_requirements_within_capacity_in_dim(flex_soln->prob, service, 
            server, i) || 
            !total_requirements_within_available_capacity_in_dim(flex_soln, 
                service, server, i)) return 0;

    return 1;
}

double **global_available_resources;
double **global_fluid_loads;

void initialize_global_resource_availabilities_and_loads(
    flexsched_problem_t flex_prob) 
{
    int i, j;

    global_available_resources = 
        (double **)calloc(flex_prob->num_servers, sizeof(double *));
    global_fluid_loads = 
        (double **)calloc(flex_prob->num_servers, sizeof(double *));

    for (i = 0; i < flex_prob->num_servers; i++) {
        global_available_resources[i] = 
            (double *)calloc(flex_prob->num_resources, sizeof(double));
        global_fluid_loads[i] = 
            (double *)calloc(flex_prob->num_resources, sizeof(double));
        for (j = 0; j < flex_prob->num_resources; j++) {
            global_available_resources[i][j] =
                flex_prob->servers[i]->total_capacities[j];
            global_fluid_loads[i][j] = 0.0;
        }
    }
}

void free_global_resource_availabilities_and_loads(
    flexsched_problem_t flex_prob) 
{
    int i;

    for (i = 0; i < flex_prob->num_servers; i++) {
        free(global_available_resources[i]);
        free(global_fluid_loads[i]);
    }
    free(global_available_resources);
    free(global_fluid_loads);
    global_available_resources = NULL;
    global_fluid_loads = NULL;

    return;
}

void put_service_on_server_fast(
    flexsched_solution_t flex_soln, int service, int server) 
{
    int i;
    
    flex_soln->mapping[service] = server;

    for (i = 0; i < flex_soln->prob->num_resources; i++) {
        global_available_resources[server][i] -=
            flex_soln->prob->services[service]->total_rigid_requirements[i];
        global_fluid_loads[server][i] +=
            flex_soln->prob->services[service]->total_fluid_needs[i];
    }

    return;
}

inline double compute_available_resource_fast(
    flexsched_solution_t flex_soln, int server, int dim)
{
    return global_available_resources[server][dim];
}

inline double compute_fluid_load_fast(
    flexsched_solution_t flex_soln, int server, int dim)
{
    return global_fluid_loads[server][dim];
}

inline int total_requirements_within_available_capacity_in_dim_fast(
    flexsched_solution_t flex_soln, int service, int server, int dim)
{
    return flex_soln->prob->services[service]->total_rigid_requirements[dim]
        <= global_available_resources[server][dim] - EPSILON;
}

int service_can_fit_on_server_fast(
    flexsched_solution_t flex_soln, int service, int server)
{
    int i;

    for (i = 0; i < flex_soln->prob->num_resources; i++)
        if (!unit_requirements_within_capacity_in_dim(flex_soln->prob, service, 
            server, i) || 
            !total_requirements_within_available_capacity_in_dim_fast(flex_soln,
                service, server, i)) return 0;

    return 1;
}

void maximize_minimum_yield_on_server(
    flexsched_solution_t flex_soln, int server)
{
    int i, j;
    double load;
    double minyield = 1.0;

    // yield limited by total needs
    for (i = 0; i < flex_soln->prob->num_resources; i++) {
        load = compute_fluid_load(flex_soln, server, i);
        if (load > 0.0) { 
            minyield = MIN(minyield, 
                compute_available_resource(flex_soln, server, i) / load);
        }
    }

    // yield limited by unit needs
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            if (flex_soln->prob->services[i]->unit_fluid_needs[j] > 0.0) {
                minyield = MIN(minyield, 
                    (flex_soln->prob->servers[server]->unit_capacities[j] -
                      flex_soln->prob->services[i]->unit_rigid_requirements[j])
                    / flex_soln->prob->services[i]->unit_fluid_needs[j]);
            }
        }
    }

    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        flex_soln->yields[i] = minyield - EPSILON;
    }

    return;
}

double compute_minimum_yield(flexsched_solution_t flex_soln)
{
    return double_array_min(flex_soln->yields, flex_soln->prob->num_services);
}

void maximize_minimum_yield(flexsched_solution_t flex_soln)
{
    int i;
    double minyield;

    for (i = 0; i < flex_soln->prob->num_servers; i++) {
        maximize_minimum_yield_on_server(flex_soln, i);
    }
    minyield = compute_minimum_yield(flex_soln);
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        flex_soln->yields[i] = minyield;
    }

    return;
}

// FIXME: organizationally I don't know if it makes more sense to put this in
// or linearprog.c
void maximize_average_yield_given_minimum(flexsched_solution_t, double);

void maximize_minimum_then_average_yield(flexsched_solution_t flex_soln)
{
    double minyield;

    maximize_minimum_yield(flex_soln);
    minyield = compute_minimum_yield(flex_soln);
    maximize_average_yield_given_minimum(flex_soln, minyield);
    return;
}

double compute_average_yield(flexsched_solution_t flex_soln)
{
    int i;
    double sumyield = 0.0;
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        sumyield += flex_soln->yields[i];
    }
    return (sumyield / flex_soln->prob->num_services);
}

double compute_utilization(flexsched_solution_t flex_soln)
{
    int i;
    double total_capacity = 0.0;
    double total_alloc = 0.0;

    for (i = 0; i < flex_soln->prob->num_servers; i++) {
        total_capacity +=
            double_array_sum(flex_soln->prob->servers[i]->total_capacities,
                flex_soln->prob->num_resources);
    }

    if (total_capacity <= 0.0) return 1.0;

    for (i = 0; i < flex_soln->prob->num_services; i++) {
        total_alloc +=
            double_array_sum(
                flex_soln->prob->services[i]->total_rigid_requirements, 
                flex_soln->prob->num_resources) + 
            flex_soln->yields[i] * 
              double_array_sum(flex_soln->prob->services[i]->total_fluid_needs,
                    flex_soln->prob->num_resources);
    }

    return (total_alloc / total_capacity);
}


/* sanity_check() 
 *
 * Checks that resource capacities are not overcome and that
 * the mapping is valid
 */
int sanity_check(flexsched_solution_t flex_soln)
{
    int i, j, k;
    int retval = 0;
    double allocated_resources[flex_soln->prob->num_servers][flex_soln->prob->num_resources];

    // check that each service is mapped to a server and yields are <= 1.0
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        if (flex_soln->mapping[i] < 0 || 
            flex_soln->mapping[i] >= flex_soln->prob->num_servers) {
            fprintf(stderr, 
                "Error: Service %d is mapped to invalid server %d\n.", i, 
                flex_soln->mapping[i]);
            retval = 1;
        }
        if (flex_soln->yields[i] < EPSILON || flex_soln->yields[i] > 1.0) {
            fprintf(stderr, "Error: Allocation of service %d is %.2f.\n",
                i, flex_soln->yields[i]);
            retval = 1;
        }
    }

    if (retval) return retval;

    // initialize allocated resources
    for (i = 0; i < flex_soln->prob->num_servers; i++) {
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            allocated_resources[i][j] = 0.0;
        }
    }

    // check unit allocations and sum up allocations
    for (i = 0; i < flex_soln->prob->num_services; i++) {
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            if(!unit_allocation_at_yield_within_capacity_in_dim(flex_soln->prob,
                i, flex_soln->yields[i], flex_soln->mapping[i], j))
            {
                fprintf(stderr, 
                    "Error: Allocation of service %d to server %d with yield %.3f exceeds unit capacity in resource dimension %d (%.3f + %.3f = %.3f/%.3f).\n", 
                    i, flex_soln->mapping[i], flex_soln->yields[i], j,
                    flex_soln->prob->services[i]->unit_rigid_requirements[j],
                    flex_soln->yields[i] * 
                        flex_soln->prob->services[i]->unit_fluid_needs[j], 
                    flex_soln->prob->services[i]->unit_rigid_requirements[j] + 
                        flex_soln->yields[i] * 
                        flex_soln->prob->services[i]->unit_fluid_needs[j], 
                    flex_soln->prob->servers[flex_soln->mapping[i]]->unit_capacities[j]);
                retval = 1;
            }
            allocated_resources[flex_soln->mapping[i]][j] += 
                flex_soln->prob->services[i]->total_rigid_requirements[j] + 
                flex_soln->yields[i] * 
                    flex_soln->prob->services[i]->total_fluid_needs[j];
        }
    }

    for (i = 0; i < flex_soln->prob->num_servers; i++) {
        for (j = 0; j < flex_soln->prob->num_resources; j++) {
            if (flex_soln->prob->servers[i]->total_capacities[j] + EPSILON <= 
                allocated_resources[i][j]) {
                fprintf(stderr, 
                    "Error: Total resource allocation on server %d exceeds capacity in dimension %d (%.6f/%.6f).\n",
                    i, j, allocated_resources[i][j], 
                    flex_soln->prob->servers[i]->total_capacities[j]);
                retval = 1;
            }
        }
    }

    return retval;
}
