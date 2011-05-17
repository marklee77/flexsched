#include "flexsched.h"

/* Keeping the services in their original order */
int GREEDY_compare_S1(const void *x, const void *y)
{
    return 0;
}

/* Sorting the services by decreasing of max fluid needs */
int GREEDY_compare_S2(const void *x, const void *y)
{
    return CMP(
        array_max(flex_prob->fluid_needs[*(int *)x], flex_prob->num_fluid),
        array_max(flex_prob->fluid_needs[*(int *)y], flex_prob->num_fluid));
}

/* Sorting the services by decreasing sum of fluid needs */
int GREEDY_compare_S3(const void *x, const void *y)
{
    return CMP(
        array_sum(flex_prob->fluid_needs[*(int *)x], flex_prob->num_fluid),
        array_sum(flex_prob->fluid_needs[*(int *)y], flex_prob->num_fluid));
}

/* Sorting the services by decreasing max of rigid and constrained fluid need */
int GREEDY_compare_S4(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float max_rigid_x = 
        array_max(flex_prob->rigid_needs[ix], flex_prob->num_rigid);
    float max_rigid_y = 
        array_max(flex_prob->rigid_needs[iy],flex_prob->num_rigid);
    float max_constrainedfluid_x = flex_prob->slas[ix] * 
        array_max(flex_prob->fluid_needs[ix], flex_prob->num_fluid);
    float max_constrainedfluid_y = flex_prob->slas[iy] * 
        array_max(flex_prob->fluid_needs[iy], flex_prob->num_fluid);

    return CMP(MAX(max_rigid_x, max_constrainedfluid_x), 
        MAX(max_rigid_y, max_constrainedfluid_y));
}

/* Sorting the services by decreasing sum of rigid and constrained fluid need */
int GREEDY_compare_S5(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float sum_rigid_x = 
        array_sum(flex_prob->rigid_needs[ix], flex_prob->num_rigid);
    float sum_rigid_y = 
        array_sum(flex_prob->rigid_needs[iy], flex_prob->num_rigid);
    float sum_constrainedfluid_x = flex_prob->slas[ix] * 
        array_sum(flex_prob->fluid_needs[ix], flex_prob->num_fluid);
    float sum_constrainedfluid_y = flex_prob->slas[iy] * 
        array_sum(flex_prob->fluid_needs[iy], flex_prob->num_fluid);

    return CMP(sum_rigid_x + sum_constrainedfluid_x, 
        sum_rigid_y + sum_constrainedfluid_y);
}

/* Sorting the services by decreasing max of all needs */
int GREEDY_compare_S6(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float max_rigid_x = 
        array_max(flex_prob->rigid_needs[ix], flex_prob->num_rigid);
    float max_rigid_y = 
        array_max(flex_prob->rigid_needs[iy], flex_prob->num_rigid);
    float max_fluid_x = 
        array_max(flex_prob->fluid_needs[ix], flex_prob->num_fluid);
    float max_fluid_y = 
        array_max(flex_prob->fluid_needs[iy], flex_prob->num_fluid);

    return CMP(MAX(max_rigid_x, max_fluid_x), MAX(max_rigid_y, max_fluid_y));
}

/* Sorting the services by decreasing sum of all needs */
int GREEDY_compare_S7(const void *x, const void *y)
{
    int ix = *((int *)x);
    int iy = *((int *)y);

    float sum_x = array_sum(flex_prob->rigid_needs[ix], flex_prob->num_rigid) + 
        array_sum(flex_prob->fluid_needs[ix],flex_prob->num_fluid);
    float sum_y = array_sum(flex_prob->rigid_needs[iy],flex_prob->num_rigid) +
        array_sum(flex_prob->fluid_needs[iy],flex_prob->num_fluid);

    return CMP(sum_x, sum_y);
}

void GREEDY_sort_services(const char *S, int *sortmap)
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
    for (i = 0; i < flex_prob->num_services; i++) sortmap[i] = i;

    // Sort the jobs
    qsort(sortmap, flex_prob->num_services, sizeof(int), compar);

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
    for (i=0; i < flex_prob->num_servers; i++) {
        if (!service_can_fit_on_server_fast(service, i)) continue;
        load = compute_server_load_in_dimension_fast(i, "fluid", array_argmax(
            flex_prob->fluid_needs[service], flex_prob->num_fluid));
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
    for (i = 0; i < flex_prob->num_servers; i++) {
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
    for (i = 0; i < flex_prob->num_servers; i++) {
        if (!service_can_fit_on_server_fast(service, i)) continue;

        if (array_max(flex_prob->rigid_needs[service], flex_prob->num_rigid) > 
            flex_prob->slas[service] * 
            array_max(flex_prob->fluid_needs[service], flex_prob->num_fluid)) {
            load = compute_server_load_in_dimension_fast(i, "rigid", 
                array_argmax(flex_prob->rigid_needs[service], 
                    flex_prob->num_rigid));
        } else {
            load = compute_server_load_in_dimension_fast(i, "fluidmin", 
                array_argmax(flex_prob->fluid_needs[service],
                    flex_prob->num_fluid));
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

    for (i=0; i < flex_prob->num_servers; i++) {

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

    for (i = 0; i < flex_prob->num_servers; i++) {
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

void GREEDY_compute_mapping(const char *P, int *sortmap, 
        flexsched_solution flex_soln)
{
    int i, j; 
    int service, server;


    for (i = 0; i < flex_prob->num_services; i++) {
        service = sortmap[i];
        server = GREEDY_pick_server(P, service);

        if (-1 == server) return;
        
        flex_soln->mapping[service] = server;

        // only needed if we're using _fast
        add_service_load_to_server(service, server);

    }

    flex_soln->success = 1;
    return;

}

flexsched_solution GREEDY_scheduler(char *S, char *P, char *ignore)
{
    int sortmap[flex_prob->num_services];
    char scheduler_name[20];
    sprintf(scheduler_name, "GREEDY_%s_%s", S, P);
    flexsched_solution flex_soln = new_flexsched_solution(scheduler_name);

    // needed for _fast
    initialize_global_server_loads();

    // Sort the services appropriately
    GREEDY_sort_services(S, sortmap);

    // Map each service to a particular host
    GREEDY_compute_mapping(P, sortmap, flex_soln);

    // greedy has to do this since it doesn't have a concept of allocation
    if (flex_soln->success) {
        maximize_minimum_then_average_yield(flex_soln);
    }

    // cleanup
    free_global_server_loads();

    return flex_soln;
}

flexsched_solution METAGREEDY_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution("METAGREEDY");

    char *sorting[] = {"S1","S2","S3","S4","S5","S6","S7",NULL};
    char *picking[] = {"P1","P2","P3","P4","P5","P6","P7",NULL};
    int i, is, ip;

    flexsched_solution curr_soln = NULL;
    float minyield, maxminyield = -1;
    float avgyield, maxavgyield = -1;

    for (is = 0; sorting[is]; is++) {
        for (ip = 0; picking[ip]; ip++) {

            // call the GREEDY scheduler
            curr_soln = GREEDY_scheduler(sorting[is], picking[ip], NULL);

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
                    for (i=0; i < flex_prob->num_services; i++) {
                        flex_soln->mapping[i] = curr_soln->mapping[i];
                        flex_soln->scaled_yields[i] = 
                            curr_soln->scaled_yields[i];
                    }
                }
            }

            free_flexsched_solution(curr_soln);
        }
    }

    return flex_soln;
}

flexsched_solution METAGREEDYLIGHT_scheduler(char *S, char *P, char *ignore)
{
    flexsched_solution flex_soln = new_flexsched_solution("METAGREEDYLIGHT");

    char *sorting[] = {"S3","S2","S5","S7","S6","S2","S3","S6","S5",NULL};
    char *picking[] = {"P1","P1","P1","P1","P1","P2","P2","P2","P4",NULL};
    int i, isp;

    flexsched_solution curr_soln = NULL;
    float minyield, maxminyield = -1;
    float avgyield, maxavgyield = -1;

    for (isp = 0; sorting[isp] && picking[isp]; isp++) {

        // call the GREEDY scheduler
        curr_soln = GREEDY_scheduler(sorting[isp], picking[isp], NULL);

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
                for (i=0; i < flex_prob->num_services; i++) {
                    flex_soln->mapping[i] = curr_soln->mapping[i];
                    flex_soln->scaled_yields[i] = curr_soln->scaled_yields[i];
                }
            }
        }

        free_flexsched_solution(curr_soln);
    }

    return flex_soln;
}
