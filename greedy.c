#include "flexsched.h"

/* Keeping the services in their original order */
int GREEDY_compare_S1(const void *x, const void *y)
{
    return 0;
}

/* Sorting the services by decreasing of max fluid needs */
int GREEDY_compare_S2(const void *x, const void *y)
{
    return CMP(array_max(INS.fluidneeds[(int *)x], INS.numfluid),
        array_max(INS.fluidneeds[(int *)y], INS.numfluid));
}

/* Sorting the services by decreasing sum of fluid needs */
int GREEDY_compare_S3(const void *x, const void *y)
{
    return CMP(array_sum(INS.fluidneeds[(int *)x], INS.numfluid),
        array_sum(INS.fluidneeds[(int *)y], INS.numfluid));
}

/* Sorting the services by decreasing max of rigid and constrained fluid need */
int GREEDY_compare_S4(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float max_rigid_x = array_max(INS.rigidneeds[ix], INS.numrigid);
    float max_rigid_y = array_max(INS.rigidneeds[iy],INS.numrigid);
    float max_constrainedfluid_x = 
        INS.slas[ix] * array_max(INS.fluidneeds[ix], INS.numfluid);
    float max_constrainedfluid_y = 
        INS.slas[iy] * array_max(INS.fluidneeds[iy], INS.numfluid);

    return CMP(MAX(max_rigid_x, max_constrainedfluid_x), 
        MAX(max_rigid_y, max_constrainedfluid_y));
}

/* Sorting the services by decreasing sum of rigid and constrained fluid need */
int GREEDY_compare_S5(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float sum_rigid_x = array_sum(INS.rigidneeds[ix], INS.numrigid);
    float sum_rigid_y = array_sum(INS.rigidneeds[iy], INS.numrigid);
    float sum_constrainedfluid_x = 
        INS.slas[ix] * array_sum(INS.fluidneeds[ix], INS.numfluid);
    float sum_constrainedfluid_y = 
        INS.slas[iy] * array_sum(INS.fluidneeds[iy], INS.numfluid);

    return CMP(sum_rigid_x + sum_constrainedfluid_x, 
        sum_rigid_y + sum_constrainedfluid_y);
}

/* Sorting the services by decreasing max of all needs */
int GREEDY_compare_S6(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float max_rigid_x = array_max(INS.rigidneeds[ix], INS.numrigid);
    float max_rigid_y = array_max(INS.rigidneeds[iy], INS.numrigid);
    float max_fluid_x = array_max(INS.fluidneeds[ix], INS.numfluid);
    float max_fluid_y = array_max(INS.fluidneeds[iy], INS.numfluid);

    return CMP(MAX(max_rigid_x, max_fluid_x), MAX(max_rigid_y, max_fluid_y));
}

/* Sorting the services by decreasing sum of all needs */
int GREEDY_compare_S7(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float sum_x = array_sum(INS.rigidneeds[ix], INS.numrigid) + 
        array_sum(INS.fluidneeds[ix],INS.numfluid);
    float sum_y = array_sum(INS.rigidneeds[iy],INS.numrigid) +
        array_sum(INS.fluidneeds[iy],INS.numfluid);

    return CMP(sum_x, sum_y);
}

void GREEDY_sort_services(const char *S, int *sorted)
{
    int i;
    int (*compar)(const void *, const void *);

    // Select the appropriate sorting routine
    if (!strcmp(S, "S1")) {
        compar = GREEDY_compare_S1;
    } else if (!strcmp(S, "S2")) {
        compar = GREEDY_compare_S2;
    } else if (!strcmp(S, "S3")) {
        compar = GREEDY_compare_S3;
    } else if (!strcmp(S, "S4")) {
        compar = GREEDY_compare_S4;
    } else if (!strcmp(S, "S5")) {
        compar = GREEDY_compare_S5;
    } else if (!strcmp(S, "S6")) {
        compar = GREEDY_compare_S6;
    } else if (!strcmp(S, "S7")) {
        compar = GREEDY_compare_S7;
    } else {
        fprintf(stderr,"Greedy algorithm: unknown sorting procedure '%s'\n",S);
        exit(1);
    }

    // Initialize index array with the original order
    for (i = 0; i < INS.numservices; i++) sorted[i]=i;

    // Sort the jobs
    qsort(sorted, INS.numservices, sizeof(int), compar);

    return;
}

int GREEDY_pick_server_P1(int service)
{
    int i;
    float load;
    float minload;
    int picked;
  
    minload = -1.0;
    picked = -1;
    for (i=0; i < INS.numservers; i++) {
        if (!service_can_fit_on_server_fast(service, i)) continue;
        load = compute_server_load_in_dimension_fast(i, "fluid", 
            array_argmax(INS.fluidneeds[service], INS.numfluid));
        if ((minload == -1.0) || (load < minload)) {
            minload = load;
            picked = i; 
        }    
    }

    return picked;
}

int GREEDY_pick_server_P2(int service)
{
    int i;
    float load;
    float minload;
    int picked;
  
    minload = -1.0;
    picked = -1;
    for (i = 0; i < INS.numservers; i++) {
        if (!service_can_fit_on_server_fast(service,i)) continue;
        load = compute_sum_server_load_fast(i, "fluid");
        if ((minload == -1.0) || (load < minload)) {
            minload = load;
            picked = i;
        }
    }

    return picked;
}

int GREEDY_pick_server_P3_P5(int service, const char *mode)
{
    int i, picked;
    float objload;
    float load;

    picked = -1;
    for (i = 0; i < INS.numservers; i++) {
        if (!service_can_fit_on_server_fast(service, i)) continue;

        if (array_max(INS.rigidneeds[service], INS.numrigid) > 
            INS.slas[service] * array_max(INS.fluidneeds[service], 
            INS.numfluid)) {
            load = compute_server_load_in_dimension_fast(i, "rigid", 
                array_argmax(INS.rigidneeds[service],INS.numrigid));
        } else {
            load = compute_server_load_in_dimension_fast(i, "fluidmin", 
                array_argmax(INS.fluidneeds[service],INS.numfluid));
        }

        if (!strcmp(mode, "bestfit")) {
            if ((picked == -1) || (load > objload)) {
                objload = load;
                picked = i;
            }
        } else if (!strcmp(mode, "worstfit")) {
            if ((picked == -1) || (load < objload)) {
                objload = load;
                picked = i;
            }
        }
    }

    return picked;
} 

int GREEDY_pick_server_P3(int service)
{
  return GREEDY_pick_server_P3_P5(service, "bestfit");
}

int GREEDY_pick_server_P5(int service)
{
  return GREEDY_pick_server_P3_P5(service, "worstfit");
}

int GREEDY_pick_server_P4_P6(int service, const char *mode)
{
    int i, picked;
    float objload = 0.0;
    float load;

    picked = -1;

    for (i=0; i < INS.numservers; i++) {

        if (!service_can_fit_on_server_fast(service,i)) continue;

        load = compute_sum_server_load_fast(i, "rigid") + 
            compute_sum_server_load_fast(i, "fluidmin");

        if (!strcmp(mode,"bestfit")) {
            if ((picked == -1) || (load > objload)) {
                objload = load;
                picked = i;
            }
        } else if (!strcmp(mode,"worstfit")) {
            if ((picked == -1) || (load < objload)) {
                objload = load;
                picked = i;
            }
        }
    }

    return picked;
}

int GREEDY_pick_server_P4(int service)
{
  return GREEDY_pick_server_P4_P6(service, "bestfit");
}

int GREEDY_pick_server_P6(int service)
{
  return GREEDY_pick_server_P4_P6(service, "worstfit");
}

int GREEDY_pick_server_P7(int service)
{
    int i, picked;

    picked = -1;

    for (i = 0; i < INS.numservers; i++) {
        if (!service_can_fit_on_server_fast(service, i)) continue;
        picked = i;
        break;
    }

    return picked;
}

int GREEDY_pick_server(const char *P, int service)
{
    if (!strcmp(P,"P1")) {
        return GREEDY_pick_server_P1(service);
    } else if (!strcmp(P,"P2")) {
        return GREEDY_pick_server_P2(service);
    } else if (!strcmp(P,"P3")) {
        return GREEDY_pick_server_P3(service);
    } else if (!strcmp(P,"P4")) {
        return GREEDY_pick_server_P4(service);
    } else if (!strcmp(P,"P5")) {
        return GREEDY_pick_server_P5(service);
    } else if (!strcmp(P,"P6")) {
        return GREEDY_pick_server_P6(service);
    } else if (!strcmp(P,"P7")) {
        return GREEDY_pick_server_P7(service);
    }

    fprintf(stderr, "Greedy algorithm: unknown picking procedure '%s'\n", P);
    exit(1);
    return -1;
}

int GREEDY_compute_mapping(const char *P, int *sorted)
{
    int i, j; 
    int service, server;


    /* Initialize the server load (for speed of computation) */
    // FIXME: do we really need to call this now? 
    initialize_global_server_loads();

    for (i = 0; i < INS.numservices; i++) {
        service = sorted[i];
        server = GREEDY_pick_server(P, service)

        if (-1 == server) return RESOURCE_ALLOCATION_FAILURE;

        /* Update the server loads */
        for (j = 0; j < INS.numrigid; j++) {
            global_server_rigid_loads[server][j] += INS.rigidneeds[service][j];
        }
        for (j = 0; j < INS.numfluid; j++) {
            global_server_fluid_loads[server][j] += INS.fluidneeds[service][j];
            global_server_fluidmin_loads[server][j] +=
                INS.slas[service] * INS.fluidneeds[service][j];
        }
    }

    return RESOURCE_ALLOCATION_SUCCESS;
}

int GREEDY_scheduler(char *S, char *P, char *ignore)
{
    int i;
    int sorted[INS.numservices];

    // Sort the services appropriately
    GREEDY_sort_services(S,sorted);

    // Map each service to a particular host
    if (GREEDY_compute_mapping(P,sorted) == RESOURCE_ALLOCATION_FAILURE) {
        return RESOURCE_ALLOCATION_FAILURE; 
    }

    // For each host, compute its services' allocations
    for (i = 0; i < INS.numservers; i++) {
        compute_allocations_given_mapping(i);
    } 

    return RESOURCE_ALLOCATION_SUCCESS;
}

int METAGREEDY_scheduler(char *S, char *P, char *ignore)
{
    char *sorting[] = {"S1","S2","S3","S4","S5","S6","S7",NULL};
    char *picking[] = {"P1","P2","P3","P4","P5","P6","P7",NULL};
    int i, is, ip;
    int status;
    int metastatus = RESOURCE_ALLOCATION_FAILURE;
    float minyield, maxminyield = -1;
    float aveyield, maxaveyield = -1;
    int *mapping = NULL;
    float *allocation = NULL;

    for (is=0; sorting[is]; is++) {
        for (ip=0; picking[ip]; ip++) {
            // initialize mappings and allocation
            for (i=0; i < INS.numservices; i++) {
                INS.mapping[i] = -1;
                INS.allocation[i] = 0.0;
            }

            // call the GREEDY scheduler
            status = GREEDY_scheduler(sorting[is], picking[ip], NULL);

            // If worked and improved, then great
            if (status == RESOURCE_ALLOCATION_SUCCESS) {
                metastatus = RESOURCE_ALLOCATION_SUCCESS;

                // compute minimum yield
                minyield = compute_minimum_yield();

                // optimize and compute maximum yield
                maximize_average_yield();
                aveyield = compute_average_yield();

                if ((maxminyield == -1) || (minyield > maxminyield) || 
                    ((minyield > maxminyield - EPSILON) && 
                    (aveyield > maxaveyield))) {
                    maxminyield = minyield;
                    maxaveyield = aveyield;
                    // allocate memory
                    if (mapping == NULL) {
                        mapping = (int *)calloc(INS.numservices, sizeof(int));
                        allocation = 
                            (float *)calloc(INS.numservices, sizeof(float));
                    }
                    // save the returned answer
                    for (i=0; i < INS.numservices; i++) {
                        mapping[i] = INS.mapping[i];
                        allocation[i] = INS.allocation[i];
                    }
                }
            }
        }
    }

    if (metastatus == RESOURCE_ALLOCATION_SUCCESS) {
        for (i=0; i < INS.numservices; i++) {
            INS.mapping[i] = mapping[i];
            INS.allocation[i] = allocation[i];
        }
    }

    free(mapping);
    free(allocation);

    return metastatus;
}

int METAGREEDYLIGHT_scheduler(char *S, char *P, char *ignore)
{
  char *sorting[] = {"S3","S2","S5","S7","S6","S2","S3","S6","S5",NULL};
  char *picking[] = {"P1","P1","P1","P1","P1","P2","P2","P2","P4",NULL};
  int i, isp;
  int status;
  int metastatus = RESOURCE_ALLOCATION_FAILURE;
  float minyield, maxminyield = -1;
  float aveyield, maxaveyield = -1;
  int *mapping = NULL;
  float *allocation = NULL;

  for (isp=0; sorting[isp]; isp++) {
    // initialize mappings and allocation
    for (i=0; i < INS.numservices; i++) {
      INS.mapping[i] = -1;
      INS.allocation[i] = 0.0;
    }
    // call the GREEDY scheduler
    status = GREEDY_scheduler(sorting[isp], picking[isp], NULL);
    // If worked and improved, then great
    if (status == RESOURCE_ALLOCATION_SUCCESS) {
      metastatus = RESOURCE_ALLOCATION_SUCCESS;

      // compute minimum yield
      minyield = compute_minimum_yield();
      // optimize and compute maximum yield
      maximize_average_yield();
      aveyield = compute_average_yield();

      maxminyield = minyield;
      maxaveyield = aveyield;
      // allocate memory
      if (mapping == NULL) {
        mapping = (int *)calloc(INS.numservices, sizeof(int));
        allocation = (float *)calloc(INS.numservices, sizeof(float));
      }
      // save the returned answer
      for (i=0; i < INS.numservices; i++) {
        mapping[i] = INS.mapping[i];
        allocation[i] = INS.allocation[i];
      }
      break;
    }
  }

  if (metastatus == RESOURCE_ALLOCATION_SUCCESS) {
    for (i=0; i < INS.numservices; i++) {
      INS.mapping[i] = mapping[i];
      INS.allocation[i] = allocation[i];
    }
  }
    free(mapping);
    free(allocation);

  return metastatus;
}


