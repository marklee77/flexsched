#include "flexsched.h"

void free_2D_int_array(int **a, int n)
{
   int i;
   for (i=0; i<n; i++) free(a[i]);
   free(a);

   return;
}

void initialize_global_server_loads()
{
    int i, j;
    for (i=0; i < INS.numservers; i++) {
        for (j=0; j < INS.numrigid; j++)  {
            global_server_rigid_loads[i][j] = 0.0;
        }
        for (j=0; j < INS.numfluid; j++)  {
            global_server_fluid_loads[i][j] = 0.0;
            global_server_fluidmin_loads[i][j] = 0.0;
        }
    }
  
    return;
}

double compute_LP_bound()
{
    int i, j;

    float total_qos_minimum, total_scaled_needs, total_capacity;
    float min_scaled_yield, scaled_yield;

    // Compute the bound directly
    min_scaled_yield = -1.0;
    for (j = 0; j < INS.numfluid; j++) {
        total_capacity = 0.0;
        for (i = 0; i < INS.numservers; i++) {
            total_capacity += INS.fluidcapacities[i][j];
        }
        total_qos_minimum = 0.0;
        total_scaled_needs = 0.0;
        for (i=0; i < INS.numservices; i++) {
            total_qos_minimum += INS.slas[i] * INS.fluidneeds[i][j];
            total_scaled_needs += (1-INS.slas[i]) * INS.fluidneeds[i][j];
        }
        scaled_yield = (total_capacity - total_qos_minimum) / 
            total_scaled_needs;
        if ((min_scaled_yield == -1.0) || (scaled_yield < min_scaled_yield)) {
            min_scaled_yield = scaled_yield;
        }
    }

    return MIN(1.0, min_scaled_yield);
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
float compute_sum_server_load(int server, const char *type)
{
    int i;
    float sumload = 0.0;

    if (!strcmp(type,"rigid")) {
        for (i = 0; i < INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            sumload += array_sum(INS.rigidneeds[i], INS.numrigid);
        }
    } else if (!strcmp(type,"fluid")) {
        for (i = 0; i < INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            sumload += array_sum(INS.fluidneeds[i], INS.numfluid);
        }
    } else if (!strcmp(type,"fluidmin")) {
        for (i = 0; i < INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            sumload += INS.slas[i] * array_sum(INS.fluidneeds[i], INS.numfluid);
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
float compute_sum_server_load_fast(int server, const char *type)
{
    if (!strcmp(type,"rigid")) {
        return array_sum(global_server_rigid_loads[server], INS.numrigid);
    } else if (!strcmp(type,"fluid")) {
        return array_sum(global_server_fluid_loads[server], INS.numfluid);
    } else if (!strcmp(type,"fluidmin")) {
        return array_sum(global_server_fluidmin_loads[server], INS.numfluid);
    }

    fprintf(stderr,"compute_sum_server_load(): unknown type '%s'\n",type);
    exit(1);
    return -1.0;
}

/*
 * compute_server_load_in_dimension:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_server_load_in_dimension(int server, const char *type, int dim)
{
    int i;
    float load = 0.0;

    if (!strcmp(type,"rigid")) {
        for (i=0; i<INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            load += INS.rigidneeds[i][dim];
        }
    } else if (!strcmp(type,"fluid")) {
        for (i=0; i<INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            load += INS.fluidneeds[i][dim];
        }
    } else if (!strcmp(type,"fluidmin")) {
        for (i=0; i<INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            load += INS.slas[i] * INS.fluidneeds[i][dim];
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
        return global_server_rigid_loads[server][dim];
    } else if (!strcmp(type,"fluid")) {
        return global_server_fluid_loads[server][dim];
    } else if (!strcmp(type,"fluidmin")) {
        return global_server_fluidmin_loads[server][dim];
    }

    fprintf(stderr, "compute_server_load_in_dimension: unknown type '%s'\n",
        type);
    exit(1);
    return -1.0;
}

int service_can_fit_on_server(int service, int server)
{
    int i, j;
    float load;

    // Rigid needs
    for (j = 0; j < INS.numrigid; j++) {
        load = 0.0;
        for (i = 0; i < INS.numservices; i++) {
            if (i == service || INS.mapping[i] != server) continue;
            load += INS.rigidneeds[i][j];
        }
        // BEWARE OF THE EPSILON
        if (INS.rigidcapacities[server][j] - load + EPSILON < 
            INS.rigidneeds[service][j]) return 0;
    }

    // Fluid needs
    for (j = 0; j < INS.numfluid; j++) {
        load = 0.0;
        for (i = 0; i < INS.numservices; i++) {
            if (i == service || INS.mapping[i] != server) continue;
            load += INS.slas[i] * INS.fluidneeds[i][j];
        }
        // BEWARE OF THE EPSILON
        if (INS.fluidcapacities[server][j] - load + EPSILON < 
            INS.slas[service] * INS.fluidneeds[service][j]) return 0;
    }

    return 1;

}

int service_can_fit_on_server_fast(int service, int server)
{
    int j;

    // Rigid needs
    for (j = 0; j < INS.numrigid; j++) {
        // BEWARE OF THE EPSILON
        if (INS.rigidcapacities[server][j] - 
            global_server_rigid_loads[server][j] + EPSILON < 
            INS.rigidneeds[service][j]) return 0;
    }

    // Fluid needs
    for (j = 0; j < INS.numfluid; j++) {
        // BEWARE OF THE EPSILON
        if (INS.fluidcapacities[server][j] -
            global_server_fluidmin_loads[server][j] + EPSILON <
            INS.slas[service] * INS.fluidneeds[service][j]) return 0;
    }

    return 1;
}

/*
 Used for GREEDY algorithms
 Right now, gives every services on the same server
 the same yield, which is the highest feasible such yield.
*/
void compute_allocations_given_mapping(int server)
{
    int i, j;
    float yield, yield_min;
    float sum1 = 0.0, sum2 = 0.0;

    yield_min = -1.0;
    for (j = 0; j < INS.numfluid; j++) {
        for (i = 0; i < INS.numservices; i++) {
            if (INS.mapping[i] != server) continue;
            sum1 += INS.slas[i] * INS.fluidneeds[i][j];
            sum2 += (1.0 - INS.slas[i]) * INS.fluidneeds[i][j];
        }
        yield = MIN(1.0, (INS.fluidcapacities[server][j] - sum1) / sum2);
        if ((yield_min == -1.0) || (yield < yield_min)) yield_min = yield;
    }

    for (i = 0; i < INS.numservices; i++) {
        if (INS.mapping[i] != server) continue;
        INS.allocation[i] = yield_min;
    }

    return;
}

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

/* compute_minimym_yield() */
float compute_minimum_yield()
{
    int i;
    float min_yield = 1.0;
    for (i = 0; i < INS.numservices; i++) {
      if (INS.allocation[i] < min_yield)
        min_yield = INS.allocation[i];
    }
    return min_yield;
}

float compute_average_yield()
{
    int i;
    float ave_yield = 0.0;
    for (i = 0; i < INS.numservices; i++) {
      ave_yield += INS.allocation[i];
    }
    return (ave_yield / INS.numservices);
}

float compute_overall_server_load(int server)
{
  float load = 0.0;
  int i,j;

  for (i=0; i < INS.numservices; i++) {
    if (INS.mapping[i] != server)
      continue;
    // rigid needs
    for (j=0; j < INS.numrigid; j++) {
      load += INS.rigidneeds[i][j];
    }
    // fluid needs
    for (j=0; j < INS.numfluid; j++) {
      load += INS.fluidneeds[i][j] * (INS.allocation[i] + INS.slas[i]*(1.0 - INS.allocation[i]));
    }
  }
  return load;

}

float compute_utilization()
{
  float overall_load = 0.0;
  float overall_capacity;
  int i;

  for (i=0; i < INS.numservers; i++) {
    overall_load += compute_overall_server_load(i);
  }

  overall_capacity = INS.numservers * (INS.numrigid * 1.0 + INS.numfluid * 1.0);

  return (overall_load / overall_capacity);
}

float compute_potential_yield_increase(int service, int server)
{
  float current_yield = INS.allocation[service];
  float min_potential_yield_increase = -1.0;
  int i, j;

  for (j=0; j<INS.numfluid; j++) {
    float free_resource = INS.fluidcapacities[server][j];
    float potential_yield_increase;
    float potential_new_yield;
    // compute free resource
    for (i=0; i<INS.numservices; i++) {
      if (INS.mapping[i] != server)
        continue;
      free_resource -= INS.fluidneeds[i][j] *
               (INS.allocation[i] + INS.slas[i]*(1.0 - INS.allocation[i]));
    }
    // compute potential new yield
    potential_new_yield = MIN(1.0, (1 / (1 - INS.slas[service])) *
            (INS.allocation[service] * (1 - INS.slas[service]) +
                  free_resource / INS.fluidneeds[service][j]));
    // compute potential yield increase
    potential_yield_increase = potential_new_yield - INS.allocation[service];
    if ((min_potential_yield_increase == -1.0) ||
        (min_potential_yield_increase > potential_yield_increase)) {
      min_potential_yield_increase = potential_yield_increase;
    }
  }

  //fprintf(stderr,"Potential yield increase for service %d on server %d: %f\n",service,server,min_potential_yield_increase);
  return min_potential_yield_increase;
}

void maximize_average_yield()
{
  int service, server;
  float min_yield = compute_minimum_yield();

  // set everybody's allocation to the minimum
  for (service=0; service < INS.numservices; service++) {
    INS.allocation[service] = min_yield;
  }

  // For each serer, increase yields
  for (server=0; server < INS.numservers; server++) {
    while(1) {
      float max_increase;
      int service_picked_for_increase;

      // pick a service for increase
      service_picked_for_increase = -1;
      for (service=0; service < INS.numservices; service++) {
        float increase;
        if (INS.mapping[service] != server)
          continue;
        increase = compute_potential_yield_increase(service,server);
        if ((service_picked_for_increase == -1) ||
            (increase > max_increase)) {
          service_picked_for_increase = service;
          max_increase = increase;
        }
      }

      // are we done?
      if (max_increase < EPSILON)
        break;

      // do the increase
      INS.allocation[service_picked_for_increase] += max_increase;
    }
  }

  return;
}

float calc_stddev(int hosts, int tasks, int *taskhost, float *taskattr) {
    int i;
    float hostload[hosts];
    float avghostload = 0.0;
    float stddev = 0.0;
    for (i = 0; i < hosts; i++) {
        hostload[i] = 0;
    }
    for (i = 0; i < tasks; i++) {
        hostload[taskhost[i]] += taskattr[i];
        avghostload += taskattr[i];
    }
    avghostload /= (float)hosts;
    for (i = 0; i < hosts; i++) {
        stddev += (hostload[i] - avghostload) * (hostload[i] - avghostload);
    }
    stddev = sqrtf(stddev / (hosts - 1));
    return stddev;
}

float find_min_memload(int hosts, int tasks, int *taskhost, float *taskmem) {
    int i;
    float memload[hosts];
    float minload = 1.0;

    for (i = 0; i < hosts; i++) {
        memload[i] = 0.0;
    }
    for (i = 0; i < tasks; i++) {
        memload[taskhost[i]] += taskmem[i];
    }
    for (i = 0; i < hosts; i++) {
        if (memload[i] < minload) {
            minload = memload[i];
        }
    }
    return minload;
}

float find_max_memload(int hosts, int tasks, int *taskhost, float *taskmem) {
    int i;
    float memload[hosts];
    float maxload = 0.0;

    for (i = 0; i < hosts; i++) {
        memload[i] = 0.0;
    }
    for (i = 0; i < tasks; i++) {
        memload[taskhost[i]] += taskmem[i];
    }
    for (i = 0; i < hosts; i++) {
        if (memload[i] > maxload) {
            maxload = memload[i];
        }
    }
    return maxload;
}

/* sanity_check() 
 *
 * Checks that resource capacities are not overcome and that
 * the mapping is valid
 */
int sanity_check()
{
  int j, server, service;
  int return_value = RESOURCE_ALLOCATION_SUCCESS;

  // check that each service is mapped to a serer
  for (service=0; service < INS.numservices; service++) {
    if ((INS.mapping[service] < 0) || (INS.mapping[service] > INS.numservers -1)) {
      fprintf(stderr,"Error: Service %d is mapped to invalid server %d\n",
                 service, INS.mapping[service]);
      return_value = RESOURCE_ALLOCATION_FAILURE;
    }
  }

  // check that server capacities are respected
  for (server=0; server < INS.numservers; server++) {
    // checking rigid needs
    for (j=0; j < INS.numrigid; j++) {
      float load=0.0;
      for (service=0; service < INS.numservices; service++) {
        if (INS.mapping[service] != server)
          continue;
        load += INS.rigidneeds[service][j];
      }
      if (load > INS.rigidcapacities[server][j] + EPSILON) {
        fprintf(stderr,
          "Error: Rigid Capacity %d of server %d exceeded (%f/%f)\n", j, server,
          load, INS.rigidcapacities[server][j]);
        return_value = RESOURCE_ALLOCATION_FAILURE;
      } else {
    //    fprintf(stderr,"Rigid load #%d = %f\n",j,load);
      }
    }
    // checking fluid needs
    for (j=0; j < INS.numfluid; j++) {
      float load=0.0;
      for (service=0; service < INS.numservices; service++) {
        if (INS.mapping[service] != server)
          continue;
        load += INS.fluidneeds[service][j] * (INS.allocation[service] +
          INS.slas[service]*(1 - INS.allocation[service]));
      }
      if (load > INS.fluidcapacities[server][j] + EPSILON) {
        fprintf(stderr,
          "Error: Fluid Capacity %d of server %d exceeded (%f/%f)\n",
          j, server, load, INS.fluidcapacities[server][j]);
        return_value = RESOURCE_ALLOCATION_FAILURE;
      } else {
     //   fprintf(stderr,"Fluid load #%d = %f\n",j,load);
      }
    }
  }
  // Check that all allocations are <= 1.0
  for (service=0; service < INS.numservices; service++) {
    if (INS.allocation[service] > 1.0) {
      fprintf(stderr,"Error: Allocation of service %d is > 1.0 (%.2f)\n",
                        service, INS.allocation[service]);
      return_value = RESOURCE_ALLOCATION_FAILURE;
    }
  }
  return return_value;
}



