#include "flexsched.h"

/* This is a global variable that keeps track of server load
*  *     and that is used to speed up a bunch of the naively implemented
*  *        algorithms (especially greedy ones) */
float **global_server_rigid_loads = NULL;
float **global_server_fluid_loads = NULL;
float **global_server_fluidmin_loads = NULL;

void initialize_global_server_loads()
{
    int i, j;

    // Allocated memory for global_server_loads
    if (!global_server_rigid_loads)
        global_server_rigid_loads = 
            (float **)calloc(flex_prob->num_servers, sizeof(float*));
    if (!global_server_fluid_loads)
        global_server_fluid_loads = 
            (float **)calloc(flex_prob->num_servers, sizeof(float*));
    if (!global_server_fluidmin_loads)
        global_server_fluidmin_loads = 
            (float **)calloc(flex_prob->num_servers, sizeof(float*));

    for (i=0; i < flex_prob->num_servers; i++) {
        if(!global_server_rigid_loads[i])
            global_server_rigid_loads[i] = 
                (float *)calloc(flex_prob->num_rigid, sizeof(float));
        for (j=0; j < flex_prob->num_rigid; j++)  {
            global_server_rigid_loads[i][j] = 0.0;
        }
        if(!global_server_fluid_loads[i])
            global_server_fluid_loads[i] = 
                (float *)calloc(flex_prob->num_fluid, sizeof(float));
        if(!global_server_fluidmin_loads[i])
            global_server_fluidmin_loads[i] = 
                (float *)calloc(flex_prob->num_fluid, sizeof(float));
        for (j=0; j < flex_prob->num_fluid; j++)  {
            global_server_fluid_loads[i][j] = 0.0;
            global_server_fluidmin_loads[i][j] = 0.0;
        }
    }
  
    return;
}

void free_global_server_loads() 
{
    int i;
    for (i = 0; i < flex_prob->num_servers; i++) {
        free(global_server_rigid_loads[i]);
        free(global_server_fluid_loads[i]);
        free(global_server_fluidmin_loads[i]);
    }
    free(global_server_rigid_loads);
    global_server_rigid_loads = NULL;
    free(global_server_fluid_loads);
    global_server_fluid_loads = NULL;
    free(global_server_fluidmin_loads);
    global_server_fluidmin_loads = NULL;
    return;
}

void add_service_load_to_server(int service, int server) {
    int j;

    for (j = 0; j < flex_prob->num_rigid; j++) {
        global_server_rigid_loads[server][j] +=  
            flex_prob->rigid_needs[service][j];
    }   
    for (j = 0; j < flex_prob->num_fluid; j++) {
        global_server_fluid_loads[server][j] +=  
            flex_prob->fluid_needs[service][j];
        global_server_fluidmin_loads[server][j] +=
            flex_prob->slas[service] * flex_prob->fluid_needs[service][j];
    }   

    return;
}

double compute_LP_bound()
{
    int i, j;

    float total_sla_need, total_non_sla_need, total_capacity;
    float yield, minyield;

    minyield = 1.0 + EPSILON;
    for (j = 0; j < flex_prob->num_fluid; j++) {
        total_capacity = 0.0;
        for (i = 0; i < flex_prob->num_servers; i++) {
            total_capacity += flex_prob->fluid_capacities[i][j];
        }
        total_sla_need = 0.0;
        total_non_sla_need = 0.0;
        for (i=0; i < flex_prob->num_services; i++) {
            total_sla_need += flex_prob->slas[i] * flex_prob->fluid_needs[i][j];
            total_non_sla_need += 
                (1 - flex_prob->slas[i]) * flex_prob->fluid_needs[i][j];
        }
        yield = MIN((total_capacity - total_sla_need) / total_non_sla_need, 
            1.0);
        if (yield < minyield) minyield = yield;
    }

    return minyield;
}

/* array_sum(): Utility function */
float array_sum(float *array, int size)
{
    int i;
    float sum = 0.0;
    for (i = 0; i < size; i++) sum += array[i];
    return sum;
}

/* array_max(): Utility function */
float array_max(float *array, int size)
{
    int i;
    float max = array[0];
    for (i = 1; i < size; i++) if (max < array[i]) max = array[i];
    return max;
}

/* array_min(): Utility function */
float array_min(float *array, int size)
{
    int i;
    float min = array[0];
    for (i = 1; i < size; i++) if (min > array[i]) min = array[i];
    return min;
}

/* array_argmax(): Utility function */
int array_argmax(float *array, int size)
{
    int i;
    int argmax = 0;
    float max = array[0];
    for (i = 1; i < size; i++) 
        if (max < array[i]) {
            argmax = i;
            max = array[i];
        }
    return argmax;
}

/*
 * compute_sum_server_load:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
// FIXME: strcmp!?
float compute_sum_server_load(flexsched_solution flex_soln, int server, 
    const char *type)
{
    int i;
    float sumload = 0.0;

    if (!strcmp(type,"rigid")) {
        for (i = 0; i < flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            sumload += 
                array_sum(flex_prob->rigid_needs[i], flex_prob->num_rigid);
        }
    } else if (!strcmp(type,"fluid")) {
        for (i = 0; i < flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            sumload += 
                array_sum(flex_prob->fluid_needs[i], flex_prob->num_fluid);
        }
    } else if (!strcmp(type,"fluidmin")) {
        for (i = 0; i < flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            sumload += flex_prob->slas[i] * 
                array_sum(flex_prob->fluid_needs[i], flex_prob->num_fluid);
        }
    } else {
        fprintf(stderr, "compute_sum_server_load(): unknown type '%s'\n",
            type);
        exit(1);
    }
    return sumload;
}

/*
 * compute_sum_server_load_fast:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
// FIXME: strcmp
float compute_sum_server_load_fast(int server, const char *type)
{
    if (!strcmp(type,"rigid")) {
        if (!global_server_rigid_loads || !global_server_rigid_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return 
            array_sum(global_server_rigid_loads[server], flex_prob->num_rigid);
    } else if (!strcmp(type,"fluid")) {
        if (!global_server_fluid_loads || !global_server_fluid_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return 
            array_sum(global_server_fluid_loads[server], flex_prob->num_fluid);
    } else if (!strcmp(type,"fluidmin")) {
        if (!global_server_fluidmin_loads || 
            !global_server_fluidmin_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return array_sum(global_server_fluidmin_loads[server], 
            flex_prob->num_fluid);
    }

    fprintf(stderr,"compute_sum_server_load(): unknown type '%s'\n",type);
    exit(1);
    return -1.0;
}

/*
 * compute_server_load_in_dimension:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_server_load_in_dimension(flexsched_solution flex_soln, int server,
    const char *type, int dim)
{
    int i;
    float load = 0.0;

    if (!strcmp(type,"rigid")) {
        for (i=0; i<flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            load += flex_prob->rigid_needs[i][dim];
        }
    } else if (!strcmp(type,"fluid")) {
        for (i=0; i<flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            load += flex_prob->fluid_needs[i][dim];
        }
    } else if (!strcmp(type,"fluidmin")) {
        for (i=0; i<flex_prob->num_services; i++) {
            if (flex_soln->mapping[i] != server) continue;
            load += flex_prob->slas[i] * flex_prob->fluid_needs[i][dim];
        }
    } else {
        fprintf(stderr, "compute_server_load_in_dimension: unknown type '%s'\n",
            type);
        exit(1);
    }

    return load;
}

/*
 * compute_server_load_in_dimension_fast:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_server_load_in_dimension_fast(int server, const char *type, 
    int dim)
{
    if (!strcmp(type,"rigid")) {
        if (!global_server_rigid_loads || !global_server_rigid_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return global_server_rigid_loads[server][dim];
    } else if (!strcmp(type,"fluid")) {
        if (!global_server_fluid_loads || !global_server_fluid_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return global_server_fluid_loads[server][dim];
    } else if (!strcmp(type,"fluidmin")) {
        if (!global_server_fluidmin_loads || 
            !global_server_fluidmin_loads[server]) {
            fprintf(stderr, 
                "need to run global initializer before _fast'%s'\n", type);
            exit(1);
        }
        return global_server_fluidmin_loads[server][dim];
    }

    fprintf(stderr, "compute_server_load_in_dimension: unknown type '%s'\n",
        type);
    exit(1);
    return -1.0;
}

int service_can_fit_on_server(
    flexsched_solution flex_soln, int service, int server)
{
    int i, j;
    float load;

    // Rigid needs
    for (j = 0; j < flex_prob->num_rigid; j++) {
        load = compute_server_load_in_dimension(flex_soln, server, "rigid", j);
        // BEWARE OF THE EPSILON
        if (flex_prob->rigid_capacities[server][j] - load + EPSILON < 
            flex_prob->rigid_needs[service][j]) return 0;
    }

    // Fluid needs
    for (j = 0; j < flex_prob->num_fluid; j++) {
        load = 
            compute_server_load_in_dimension(flex_soln, server, "fluidmin", j);
        // BEWARE OF THE EPSILON
        if (flex_prob->fluid_capacities[server][j] - load + EPSILON < 
            flex_prob->slas[service] * flex_prob->fluid_needs[service][j]) 
            return 0;
    }

    return 1;

}

int service_can_fit_on_server_fast(int service, int server)
{
    int j;

    if (!global_server_rigid_loads || !global_server_rigid_loads[server] ||
        !global_server_fluidmin_loads || 
        !global_server_fluidmin_loads[server]) {
        fprintf(stderr, 
            "need to run global initializer before _fast\n");
        exit(1);
    }

    // Rigid needs
    for (j = 0; j < flex_prob->num_rigid; j++) {
        // BEWARE OF THE EPSILON
        if (flex_prob->rigid_capacities[server][j] - 
            global_server_rigid_loads[server][j] + EPSILON < 
            flex_prob->rigid_needs[service][j]) return 0;
    }

    // Fluid needs
    for (j = 0; j < flex_prob->num_fluid; j++) {
        // BEWARE OF THE EPSILON
        if (flex_prob->fluid_capacities[server][j] -
            global_server_fluidmin_loads[server][j] + EPSILON <
            flex_prob->slas[service] * flex_prob->fluid_needs[service][j]) 
            return 0;
    }

    return 1;
}

void maximize_minimum_yield_on_server(flexsched_solution flex_soln, int server)
{
    int i, j;
    float total_sla_need, total_non_sla_need, available_resource;
    float yield, minyield;

    minyield = 1.0 + EPSILON;
    for (j = 0; j < flex_prob->num_fluid; j++) {
        total_sla_need = compute_server_load_in_dimension(flex_soln, server,
            "fluidmin", j);
        total_non_sla_need = compute_server_load_in_dimension(flex_soln, server,
            "fluid", j) - total_sla_need;
        available_resource = 
            flex_prob->fluid_capacities[server][j] - total_sla_need;
        yield = MIN(available_resource / total_non_sla_need, 1.0);
        if (yield < minyield) minyield = yield;
    }


    for (i = 0; i < flex_prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        flex_soln->scaled_yields[i] = minyield;
    }

    return;
}

/* compute_minimym_yield() */
float compute_minimum_yield(flexsched_solution flex_soln)
{
    return array_min(flex_soln->scaled_yields, flex_prob->num_services);
}

void maximize_minimum_yield(flexsched_solution flex_soln)
{
    int i;
    float minyield;

    for (i = 0; i < flex_prob->num_servers; i++) {
        maximize_minimum_yield_on_server(flex_soln, i);
    }
    minyield = compute_minimum_yield(flex_soln);
    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->scaled_yields[i] = minyield;
    }

    return;
}

// FIXME: go over this from paper again and make sure you agree.
void maximize_average_yield_on_server_given_minimum(
    flexsched_solution flex_soln, int server, float minyield)
{
    int i, j;
    float available_resources[flex_prob->num_fluid];

    for (i = 0; i < flex_prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;
        flex_soln->scaled_yields[i] = minyield;
    }

    // FIXME: doesn't do anything yet...
    /*
    for (j = 0; j < flex_prob->num_fluid; j++) {
        available_resources[j] = flex_prob->fluid_capacities[j] -
            compute_load_on_server_in_dimension(flex_soln, server, "fluidmin", 
                j);
    */

    return;
}

void maximize_minimum_then_average_yield(flexsched_solution flex_soln)
{
    int i;
    float minyield;

    for (i = 0; i < flex_prob->num_servers; i++) {
        maximize_minimum_yield_on_server(flex_soln, i);
    }
    minyield = compute_minimum_yield(flex_soln);
    for (i = 0; i < flex_prob->num_servers; i++) {
        maximize_average_yield_on_server_given_minimum(flex_soln, i, minyield);
    }

    return;
}

float compute_average_yield(flexsched_solution flex_soln)
{
    int i;
    float sumyield = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        sumyield += flex_soln->scaled_yields[i];
    }
    return (sumyield / flex_prob->num_services);
}

float compute_server_sum_alloc(flexsched_solution flex_soln, 
    int server)
{
    float alloc = 0.0;
    int i, j;

    for (i = 0; i < flex_prob->num_services; i++) {
        if (flex_soln->mapping[i] != server) continue;

        // rigid needs
        for (j = 0; j < flex_prob->num_rigid; j++) 
            alloc += flex_prob->rigid_needs[i][j];
        
        // fluid needs
        for (j=0; j < flex_prob->num_fluid; j++) {
            alloc += flex_prob->fluid_needs[i][j] * (flex_prob->slas[i] + 
                flex_soln->scaled_yields[i] * (1.0 - flex_prob->slas[i]));
        }
    }
    return alloc;
}

float compute_utilization(flexsched_solution flex_soln)
{
    int i;
    float total_alloc = 0.0;
    float total_capacity = 0.0;

    for (i = 0; i < flex_prob->num_servers; i++) {
        total_alloc += compute_server_sum_alloc(flex_soln, i);
        total_capacity += 
            array_sum(flex_prob->rigid_capacities[i], flex_prob->num_rigid) +
            array_sum(flex_prob->fluid_capacities[i], flex_prob->num_fluid);
    }

    return (total_alloc / total_capacity);
}


/* sanity_check() 
 *
 * Checks that resource capacities are not overcome and that
 * the mapping is valid
 */
int sanity_check(flexsched_solution flex_soln)
{
    int i, j, k;
    int return_value = 0;
    float alloc;

    // check that each service is mapped to a server and yields are <= 1.0
    for (i = 0; i < flex_prob->num_services; i++) {
        if (flex_soln->mapping[i] < 0 || 
            flex_soln->mapping[i] >= flex_prob->num_servers) {
            fprintf(stderr, 
                "Error: Service %d is mapped to invalid server %d\n", i, 
                flex_soln->mapping[i]);
            return_value = 1;
        }
        if (flex_soln->scaled_yields[i] > 1.0) {
            fprintf(stderr, "Error: Allocation of service %d is > 1.0 (%.2f)\n",
                i, flex_soln->scaled_yields[i]);
            return_value = 1;
        }
    }

    // check that server capacities are respected
    for (j = 0; j < flex_prob->num_servers; j++) {
        // checking rigid needs
        for (k = 0; k < flex_prob->num_rigid; k++) {
            alloc = compute_server_load_in_dimension(flex_soln, j, "rigid", k);
            if(flex_prob->rigid_capacities[j][k] + EPSILON < alloc) {
                fprintf(stderr,
                    "Error: Rigid Capacity %d of server %d exceeded (%f/%f)\n", 
                    k, j, alloc, flex_prob->rigid_capacities[j][k]);
                return_value = 1;
            }
        }

        // checking fluid needs
        // FIXME: I really don't like how this is done inconsistently
        // FIXME: get a compute_sum_alloc_in_dimension?
        for (k = 0; k < flex_prob->num_fluid; k++) {
            alloc = 0.0;
            for (i = 0; i < flex_prob->num_services; i++) {
                if (flex_soln->mapping[i] != j) continue;
                alloc += flex_prob->fluid_needs[i][k] * (flex_prob->slas[i] + 
                    flex_soln->scaled_yields[i] * (1.0 - flex_prob->slas[i]));
            }
            if (flex_prob->fluid_capacities[j][k] + EPSILON < alloc) {
                fprintf(stderr,
                    "Error: Fluid Capacity %d of server %d exceeded (%f/%f)\n",
                    k, j, alloc, flex_prob->fluid_capacities[j][k]);
                return_value = 1;
            }
        }
    }

    return return_value;

}

flexsched_solution new_flexsched_solution(const char *algorithm) 
{
    int i;
    flexsched_solution flex_soln =
        (flexsched_solution)calloc(1, sizeof(struct flexsched_solution_struct));

    if (!flex_soln) {
        fprintf(stderr, 
                "couldn't allocate sufficient memory for new solution!\n");
        exit(1);
    }

    flex_soln->algorithm = strdup(algorithm);
    if (!(flex_soln->mapping = 
        (int *)calloc(flex_prob->num_services, sizeof(int)))) {
        fprintf(stderr, 
                "couldn't allocate sufficient memory for new mapping!\n");
        exit(1);
    }
    if (!(flex_soln->scaled_yields =
        (float *)calloc(flex_prob->num_services, sizeof(float)))) {
        fprintf(stderr, "couldn't allocate sufficient memory for yields!\n");
        exit(1);
    }
    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->mapping[i] = -1;
        flex_soln->scaled_yields[i] = 0.0;
    }
    return flex_soln;
}

void free_flexsched_solution(flexsched_solution flex_soln)
{
    if (!flex_soln) return;
    free(flex_soln->algorithm);
    free(flex_soln->mapping);
    free(flex_soln->scaled_yields);
    free(flex_soln);
    return;
}

#if 0
/* increment_binary_counter()
 *   returns the sum of the bits
 */
void increment_binary_counter(struct binary_counter_t *counter)
{
  int i;

  for (i = counter->size - 1; i >= 0; i--) {
    counter->bits[i] = 1 - counter->bits[i];
    if (counter->bits[i] == 1) {
      counter->sum++;
      break;
    } else {
      counter->sum--;
    }
  }
  return;
}

void print_binary_counter(struct binary_counter_t *counter)
{
  int i;
  for (i = 0; i < counter->size; i++)
    fprintf(stderr,"%d",counter->bits[i]);
  return;
}

/* find_subset_of_size() */
int *find_subset_of_size(struct vp_instance *vp,
                         int *remaining_unmapped_vectors,
                         int num_remaining_unmapped_vectors,
                         int subset_size)
{
  struct binary_counter_t binary_counter;
  int i, j, count;
  float load[vp->num_dims];
  int *found_subset;

  binary_counter.size = num_remaining_unmapped_vectors;
  binary_counter.bits = (char *)calloc(binary_counter.size,sizeof(char));
  binary_counter.sum = 0;

  found_subset = NULL;

  // initialize the binary counter to 0..01..1
  for (i=binary_counter.size - subset_size; i < binary_counter.size; i++) {
    binary_counter.bits[i] = 1;
  }
  binary_counter.sum = subset_size;
//  fprintf(stderr,"  initial value: ");
//  print_binary_counter(&binary_counter);
//  fprintf(stderr,"    (sum = %d)\n",binary_counter.sum);


  while(binary_counter.sum) {

    // Do we have the right number of bits?
    if (binary_counter.sum == subset_size) {

//      fprintf(stderr,"  current value: ");
//      print_binary_counter(&binary_counter);
//      fprintf(stderr,"    (sum = %d)\n",binary_counter.sum);
//      fprintf(stderr,"  computing load...");
      // Compute the load due to all vectors in the subset
      for (j=0; j < vp->num_dims; j++)
        load[j] = 0.0;
      for (i=0; i < binary_counter.size; i++) {
        if (!binary_counter.bits[i])
          continue;
        for (j=0; j < vp->num_dims; j++) {
          load[j] += vp->vectors[remaining_unmapped_vectors[i]].x[j];
        }
      }
//      fprintf(stderr,"  load = ");
//      for (j=0; j < vp->num_dims; j++)
//        fprintf(stderr,"%.3f ",load[j]);
//      fprintf(stderr,"\n");

      // See if the whole subset can fit in a single bin
      int fit = 1;
      for (j=0; j < vp->num_dims; j++) {
        if (load[j] > 1.0) {
          fit = 0;
          break;
        }
      }

      if (fit) {
//        fprintf(stderr,"  it fits!!\n");
        found_subset = (int *)calloc(subset_size,sizeof(int));
        count = 0;
        for (i = 0; i < binary_counter.size; i++) {
          if (binary_counter.bits[i]) {
            found_subset[count++] = remaining_unmapped_vectors[i];
          }
        }
        break;
     } else {
//        fprintf(stderr,"  it doesn't fit\n");
      }

    }
    increment_binary_counter(&binary_counter);
  }

  free(binary_counter.bits);
  return found_subset;
}

/* find_maximum_subset() 
 *
 *    Maximum size of the subset: vp->num_dims
 */
int find_maximum_subset(struct vp_instance *vp,
                        int *remaining_unmapped_vectors,
                        int num_remaining_unmapped_vectors,
                        int **subset,
                        int max_possible_subset)
{
  int subset_size;
  int *found_subset=NULL;
  int max_subset_size;

//  fprintf(stderr,"*############# %d %d\n",vp->num_dims, num_remaining_unmapped_vectors);
  if (max_possible_subset == -1)
    max_subset_size = MIN(vp->num_dims, num_remaining_unmapped_vectors);
  else
    max_subset_size = MIN(max_possible_subset, MIN(vp->num_dims, num_remaining_unmapped_vectors));

  for (subset_size = max_subset_size;
       subset_size >= 1; subset_size--) {
#ifdef CHEKURI_DEBUG
    fprintf(stderr,"  Looking for a subset of size %d\n",subset_size);
#endif
    if (found_subset = find_subset_of_size(vp,remaining_unmapped_vectors,
                                num_remaining_unmapped_vectors,subset_size)) {
#ifdef CHEKURI_DEBUG
      fprintf(stderr,"  Found one!\n");
#endif
      *subset = found_subset;
      return subset_size;
    }
  }
  fprintf(stderr,"CHEKURI: find_maximum_subset(): Can't find a subset\n");
 return 0;
}
#endif

