#include "flexsched.h"

struct implemented_scheduler_s {
    char *name;
    flexsched_solution_t (*func)(flexsched_problem_t, char *, char **);
    int sanitycheck;
};

struct implemented_scheduler_s implemented_schedulers[] = {
    {"GREEDY",          GREEDY_scheduler,          1},
    {"METAGREEDY",      METAGREEDY_scheduler,      1},
    {"METAGREEDYLIGHT", METAGREEDYLIGHT_scheduler, 1},
    {"MILP",            MILP_scheduler,            1},
    {"LPBOUND",         LPBOUND_scheduler,         0},
    {"RRND",            LPROUNDING_scheduler,      1},
    {"RRNZ",            LPROUNDING_scheduler,      1},
/*
    {"VP",              VP_scheduler,              1},
    {"HVP",             HVP_scheduler,             1},
    {"METAHVP",         HVP_scheduler,             1},
    {"OLDMETAHVP",      METAHVP_scheduler,         1},
*/
    {NULL, NULL, 0}
};

typedef struct call_scheduler_s {
    char *string;
    flexsched_solution_t (*func)(flexsched_problem_t, char *, char **);
    char *name;
    char **options;
    int active;
    int sanitycheck;
} *call_scheduler_t;

/* parse_scheduler_list() 
 *    Create a NULL-terminated list of scheduler names from
 *    the list command-line argument
 */
call_scheduler_t * parse_scheduler_list(char *scheduler_list) {
    const char *ssep = " \t", *osep = "_";
    call_scheduler_t *schedulers = NULL, *sched;
    char *s;
    int i, j;

    schedulers = (call_scheduler_t *)malloc(sizeof(call_scheduler_t));
    *schedulers = NULL;
    i = 0;
    for (s = strtok(strdup(scheduler_list), ssep); s; s = strtok(NULL, ssep)) {
        schedulers[i] = 
            (call_scheduler_t)malloc(sizeof(struct call_scheduler_s));
        schedulers[i]->string = s;
        i++;
        schedulers = (call_scheduler_t *)realloc(schedulers, 
            (i+1) * sizeof(call_scheduler_t));
        schedulers[i] = NULL;
    }
    // kind of stilly, but not all systems offer reentrant strtok
    for (sched = schedulers; *sched; sched++) {
        (*sched)->name = strtok(strdup((*sched)->string), osep);
        (*sched)->func = NULL;
        (*sched)->sanitycheck = 0;
        (*sched)->active = 0;
        for (j = 0; implemented_schedulers[j].name; j++) {
            if (!strcmp((*sched)->name, implemented_schedulers[j].name)) {
                (*sched)->func = implemented_schedulers[j].func;
                (*sched)->sanitycheck = implemented_schedulers[j].sanitycheck;
                (*sched)->active = 1;
                break;
            }
        }
        if (!(*sched)->active) {
            fprintf(stderr, "Warning: Scheduler %s not found!\n", 
                (*sched)->name);
        }
        (*sched)->options = (char **)malloc(sizeof(char *));
        i = 0;
        for ((*sched)->options[0] = strtok(NULL, osep); (*sched)->options[i];
                (*sched)->options[i] = strtok(NULL, osep)) {
            i++;
            (*sched)->options = 
                (char **)realloc((*sched)->options, (i+1) * sizeof(char *));
        }
    }

    return schedulers;
}

void deactivate_unneeded_schedulers(call_scheduler_t *schedulers, char *file)
{
    FILE *f;
    char buffer[1024];
    call_scheduler_t *sched;

    // If not output file, then fine
    if (!(f = fopen(file,"r"))) return; 

    while (fgets(buffer, 1024, f)) {
        for (sched = schedulers; *sched; sched++) {
            if (!strcmp(strtok(buffer, "|"), (*sched)->string)) {
                (*sched)->active = 0;
                fprintf(stderr, "Warning: no need to re-run algorithm '%s'\n", 
                    (*sched)->string);
            }
        }
    }

    fclose(f);

    return;
}

server_t new_server(int num_resources) {
    server_t server;
    int i;

    if (!(server = (server_t)malloc(sizeof(struct server_s)))) {
        fprintf(stderr, "Insufficient memory to allocate server!\n");
        exit(1);
    }

    if (!(server->unit_capacities = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate server unit capacities!\n");
        exit(1);
    }

    if (!(server->total_capacities = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate server total capacities!\n");
        exit(1);
    }

    return server;
}

service_t new_service(int num_resources) {
    service_t service;
    int i;

    if (!(service = (service_t)malloc(sizeof(struct service_s)))) {
        fprintf(stderr, "Insufficient memory to allocate service!\n");
        exit(1);
    }

    if (!(service->unit_rigid_requirements = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate service unit rigid_requirements!\n");
        exit(1);
    }

    if (!(service->unit_fluid_needs = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate service unit fluid_needs!\n");
        exit(1);
    }

    if (!(service->total_rigid_requirements = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate service total rigid_requirements!\n");
        exit(1);
    }

    if (!(service->total_fluid_needs = (float *)calloc(num_resources,
                    sizeof(float)))) {
        fprintf(stderr, "Insufficient memory to allocate service total fluid_needs!\n");
        exit(1);
    }

    return service;
}
    
flexsched_problem_t new_flexsched_problem(FILE *input) {
    flexsched_problem_t flex_prob;
    int i, j;

    // initialize flex_problem
    if (!(flex_prob = 
        (flexsched_problem_t)malloc(sizeof(struct flexsched_problem_s)))) {
        fprintf(stderr, 
            "Insufficient memory to allocate flexsched problem structure!\n");
        exit(1);
    }
    /* Parse instance file */
    fscanf(input,"%d", &flex_prob->num_resources);
    fscanf(input,"%d", &flex_prob->num_servers);
    fscanf(input,"%d", &flex_prob->num_services);

    flex_prob->servers = 
        (server_t *)calloc(flex_prob->num_servers, sizeof(server_t));

    for (i = 0; i < flex_prob->num_servers; i++) {
        flex_prob->servers[i] = new_server(flex_prob->num_resources);
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", &(flex_prob->servers[i]->unit_capacities[j]));
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", &(flex_prob->servers[i]->total_capacities[j]));
    }

    flex_prob->services = 
        (service_t *)calloc(flex_prob->num_services, sizeof(service_t));

    for (i = 0; i < flex_prob->num_services; i++) {
        flex_prob->services[i] = new_service(flex_prob->num_resources);
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", 
                &(flex_prob->services[i]->unit_rigid_requirements[j]));
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", &(flex_prob->services[i]->unit_fluid_needs[j]));
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", 
                &(flex_prob->services[i]->total_rigid_requirements[j]));
        for (j = 0; j < flex_prob->num_resources; j++) 
            fscanf(input,"%f", &(flex_prob->services[i]->total_fluid_needs[j]));
    }

    return flex_prob;
}

void free_server(server_t server) {
    if (!server) return;
    free(server->unit_capacities);
    free(server->total_capacities);
    free(server);
    return;
}

void free_service(service_t service) {
    if (!service) return;
    free(service->unit_rigid_requirements);
    free(service->unit_fluid_needs);
    free(service->total_rigid_requirements);
    free(service->total_fluid_needs);
    free(service);
    return;
}

void free_flexsched_problem(flexsched_problem_t flex_prob) {
    int i;

    if (!flex_prob) return;

    for (i = 0; i < flex_prob->num_servers; i++) {
        free_server(flex_prob->servers[i]);
    }
    free(flex_prob->servers);
    for (i = 0; i < flex_prob->num_services; i++) {
        free_service(flex_prob->services[i]);
    }
    free(flex_prob->services);
    free(flex_prob);
}

int main(int argc, char *argv[])
{
    FILE *input, *output;
    call_scheduler_t *schedulers, *sched;

    int i, j;
    struct timeval time1, time2;
    flexsched_problem_t flex_prob = NULL;
    flexsched_solution_t flex_soln = NULL;
    double elasped_seconds;
    double non_optimized_average_yield;
    double non_optimized_utilization;
    double total;

    if ((argc < 2) || (argc > 4)) {
        fprintf(stderr,
            "Usage: %s <scheduler list> [<instance file> [result file]]\n",
            argv[0]);
        fprintf(stderr,"  Known schedulers:  "); 
        for (i = 0; implemented_schedulers[i].name; i++) {
            fprintf(stderr, "%s ", implemented_schedulers[i].name);
        }
        fprintf(stderr, "\n  See README file for scheduler descriptions\n");
	    exit(1);
    }

    // schedulers
    schedulers = parse_scheduler_list(argv[1]);

    // input
    if (argc < 3) {
	    input = stdin;
    } else if (!(input = fopen(argv[2], "r"))) {
	    fprintf(stderr, "Can't open file '%s' for reading\n", argv[2]);
	    exit(1);
    }

    // output 
    if (argc < 4) {
	    output = stdout;
    } else {
        deactivate_unneeded_schedulers(schedulers, argv[3]);
	    if (!(output = fopen(argv[3],"a"))) {
	        fprintf(stderr, "Can't open file '%s' for writing/appending\n", 
                argv[3]);
	        exit(1);
	    }
    }

    flex_prob = new_flexsched_problem(input);

    fclose(input);

#if 0
    printf("FlexSched Solver\n");
    printf("servers: %d services: %d\n", flex_prob->num_servers,
        flex_prob->num_services);
    printf("total resources:");
    for (j = 0; j < flex_prob->num_rigid; j++) {
        total = 0.0;
        for(i = 0; i < flex_prob->num_servers; i++) {
            total += flex_prob->rigid_capacities[i][j];
        }
        printf(" %.3f", total);
    }
    for (j = 0; j < flex_prob->num_fluid; j++) {
        total = 0.0;
        for(i = 0; i < flex_prob->num_servers; i++) {
            total += flex_prob->fluid_capacities[i][j];
        }
        printf(" %.3f", total);
    }
    printf("\n");
    printf("total requirements:");
    for (j = 0; j < flex_prob->num_rigid; j++) {
        total = 0.0;
        for(i = 0; i < flex_prob->num_services; i++) {
            total += flex_prob->rigid_needs[i][j];
        }
        printf(" %.3f", total);
    }
    for (j = 0; j < flex_prob->num_fluid; j++) {
        total = 0.0;
        for(i = 0; i < flex_prob->num_services; i++) {
            total += flex_prob->fluid_needs[i][j];
        }
        printf(" %.3f", total);
    }
    printf("\n");
#endif

    for (sched = schedulers; *sched; sched++) {

        // if not needed, then skip
        if (!(*sched)->active) continue;

        // if we did find it, then run it
        gettimeofday(&time1, NULL);

        // Call the scheduler
        flex_soln = (*sched)->func(flex_prob, (*sched)->name, (*sched)->options);

        gettimeofday(&time2, NULL);

        if (flex_soln->success) {

            // Sanity check the allocation
            if ((*sched)->sanitycheck && sanity_check(flex_soln)) {
                fprintf(stderr,"Invalid allocation\n");
                exit(1);
            }

            non_optimized_average_yield = compute_average_yield(flex_soln);
            non_optimized_utilization = compute_utilization(flex_soln);
            maximize_minimum_then_average_yield(flex_soln);

            // Re-Sanity check the allocation
            if ((*sched)->sanitycheck && sanity_check(flex_soln)) {
                fprintf(stderr,"Invalid allocation\n");
                exit(1);
            }

        }

        // Compute elapsed time
        elasped_seconds  = (double)(time2.tv_sec - time1.tv_sec);
        elasped_seconds += (double)(time2.tv_usec - time1.tv_usec) / 1000000.0;
  
        // Print output
        fprintf(output, "%s|", (*sched)->string);
        if (!strcmp((*sched)->name,"LPBOUND")) {
            fprintf(output, "%.3f|", compute_minimum_yield(flex_soln));
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
        } else if (flex_soln->success) {
            fprintf(output, "%.3f|", compute_minimum_yield(flex_soln));
            fprintf(output, "%.3f|", non_optimized_average_yield);
            fprintf(output, "%.3f|", compute_average_yield(flex_soln));
            fprintf(output, "%.3f|", non_optimized_utilization);
            fprintf(output, "%.3f|", compute_utilization(flex_soln));
        } else {
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);
            fprintf(output, "%.3f|", -1.0);

        }
        fprintf(output, "%.3f|", elasped_seconds);
        fprintf(output, "%s\n", flex_soln->misc_output);

        // clean up
        free_flexsched_solution(flex_soln);

    }

    fclose(output);

    free_flexsched_problem(flex_prob);

    if (*schedulers) {
        free((*schedulers)->string);
    }

    for (sched = schedulers; *sched; sched++) {
        free((*sched)->name);
        free((*sched)->options);
        free(*sched);
    }
    free(schedulers);

    return 0;
}
