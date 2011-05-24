#include "flexsched.h"

flexsched_problem flex_prob = NULL;

struct implemented_scheduler_t {
    char *name;
    flexsched_solution (*func)(char *, char **);
    int sanitycheck;
};

struct implemented_scheduler_t implemented_schedulers[] = {
    {"GREEDY",          GREEDY_scheduler,          1},
    {"METAGREEDY",      METAGREEDY_scheduler,      1},
    {"METAGREEDYLIGHT", METAGREEDYLIGHT_scheduler, 1},
    {"MILP",            MILP_scheduler,       1},
    {"LPBOUND",         LPBOUND_scheduler,    0},
    {"RRND",            LPROUNDING_scheduler, 1},
    {"RRNZ",            LPROUNDING_scheduler, 1},
    {"VP",              VP_scheduler,         1},
/*
    {"HVP",             HVP_scheduler,        1},
    {"VP_CHEKURI",   VP_scheduler,         "CHEKURI", NULL, NULL},
    {"GA",           GA_scheduler,         "N", "N", "N"},
    {"GAF",          GA_scheduler,         "F", "N", "N"},
    {"GAFF",         GA_scheduler,         "F", "F", "N"},
    {"GAFFF",        GA_scheduler,         "F", "F", "F"},
*/
    {NULL, NULL, 0}
};

struct scheduler_t {
    char *string;
    flexsched_solution (*func)(char *, char **);
    char *name;
    char **options;
    int active;
    int sanitycheck;
};

/* parse_scheduler_list() 
 *    Create a NULL-terminated list of scheduler names from
 *    the list command-line argument
 */
struct scheduler_t **parse_scheduler_list(char *schedulers_string) {
    struct scheduler_t **list = NULL;
    struct scheduler_t **sched;
    char *s;
    char *ssep = " \t";
    char *osep = "_";
    int i, j;

    list = (struct scheduler_t **)malloc(sizeof(struct scheduler_t *));
    *list = NULL;
    i = 0;
    for (s = strtok(strdup(schedulers_string), ssep); s; s = strtok(NULL, ssep))
    {
        list[i] = (struct scheduler_t *)malloc(sizeof(struct scheduler_t));
        list[i]->string = s;
        i++;
        list = (struct scheduler_t**)
            REALLOC(list, (i+1) * sizeof(struct scheduler_t *));
        list[i] = NULL;
    }
    for (sched = list; *sched; sched++) {
        (*sched)->name = strtok(strdup((*sched)->string), osep);
        (*sched)->func = NULL;
        (*sched)->active = 1;
        for (j = 0; implemented_schedulers[j].name; j++) {
            if (!strcmp((*sched)->name, implemented_schedulers[j].name)) {
                (*sched)->func = implemented_schedulers[j].func;
                (*sched)->sanitycheck = implemented_schedulers[j].sanitycheck;
                break;
            }
        }
        (*sched)->options = (char **)malloc(sizeof(char *));
        i = 0;
        for ((*sched)->options[0] = strtok(NULL, osep); (*sched)->options[i];
                (*sched)->options[i] = strtok(NULL, osep)) 
        {
            i++;
            (*sched)->options = 
                (char **)REALLOC((*sched)->options, (i+1) * sizeof(char *));
        }
    }

    return list;
}

void deactivate_unneeded_schedulers(struct scheduler_t **list, char *file)
{
    FILE *f;
    char buffer[1024];
    struct scheduler_t **sched;

    // If not output file, then fine
    if (!(f = fopen(file,"r"))) return; 

    while (fgets(buffer, 1024, f)) {
        for (sched = list; *sched; sched++) {
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

int main(int argc, char *argv[])
{
    FILE *input, *output;
    struct scheduler_t **schedulers, **sched;

    int i, j;
    struct timeval time1, time2;
    flexsched_solution flex_soln = NULL;
    double elasped_seconds;
    double non_optimized_average_yield;
    double non_optimized_utilization;

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

    // initialize problem
    if (!(flex_prob = 
        (flexsched_problem)calloc(1, sizeof(struct flexsched_problem_struct)))) 
    {
        fprintf(stderr, 
            "Insufficient memory to allocate flexsched problem structure!\n");
        exit(1);
    }

    /* Parse instance file */
    fscanf(input,"%d", &flex_prob->num_rigid);
    fscanf(input,"%d", &flex_prob->num_fluid);
    fscanf(input,"%d", &flex_prob->num_servers);
    fscanf(input,"%d", &flex_prob->num_services);

    flex_prob->rigid_capacities = 
        (float **)calloc(flex_prob->num_servers, sizeof(float *));
    flex_prob->fluid_capacities = 
        (float **)calloc(flex_prob->num_servers, sizeof(float *));
    for (i = 0; i < flex_prob->num_servers; i++) {
        flex_prob->rigid_capacities[i] = 
            (float *)calloc(flex_prob->num_rigid, sizeof(float *));
        flex_prob->fluid_capacities[i] = 
            (float *)calloc(flex_prob->num_fluid, sizeof(float *));
    }
    flex_prob->slas = 
        (float *)calloc(flex_prob->num_services, sizeof(float));
    flex_prob->rigid_needs = 
        (float **)calloc(flex_prob->num_services, sizeof(float*));
    flex_prob->fluid_needs = 
        (float **)calloc(flex_prob->num_services, sizeof(float*));
    for (i = 0; i < flex_prob->num_services; i++) {
      flex_prob->rigid_needs[i] = 
          (float *)calloc(flex_prob->num_rigid, sizeof(float));
      flex_prob->fluid_needs[i] = 
          (float *)calloc(flex_prob->num_fluid,sizeof(float));
    }

    for (i = 0; i < flex_prob->num_servers; i++) {
        for (j=0; j < flex_prob->num_rigid; j++)
            fscanf(input,"%f", &(flex_prob->rigid_capacities[i][j]));
        for (j=0; j < flex_prob->num_fluid; j++)
            fscanf(input,"%f", &(flex_prob->fluid_capacities[i][j]));
    }
    for (i = 0; i < flex_prob->num_services; i++) {
      fscanf(input, "%f", &(flex_prob->slas[i]));
      for (j=0; j<flex_prob->num_rigid; j++) 
        fscanf(input, "%f", &(flex_prob->rigid_needs[i][j]));
      for (j=0; j<flex_prob->num_fluid; j++) 
        fscanf(input, "%f", &(flex_prob->fluid_needs[i][j]));
    }

    fclose(input);

    flex_prob->lpbound = compute_LP_bound();

    for (sched = schedulers; *sched; sched++) {

        // if not needed, then skip
        if (!(*sched)->active) continue;

        // if we did find it, then run it
        gettimeofday(&time1, NULL);

        // Call the scheduler
        flex_soln = (*sched)->func((*sched)->name, (*sched)->options);

        gettimeofday(&time2, NULL);

        if (flex_soln->success && (*sched)->sanitycheck) {

            // Sanity check the allocation
            if (sanity_check(flex_soln)) {
                fprintf(stderr,"Invalid allocation\n");
                exit(1);
            }

            non_optimized_average_yield = compute_average_yield(flex_soln);
            non_optimized_utilization = compute_utilization(flex_soln);
            maximize_minimum_then_average_yield(flex_soln);

            // Re-Sanity check the allocation
            if (sanity_check(flex_soln)) {
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
            fprintf(output, "%.3f|", flex_prob->lpbound);
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

    for (i = 0; i < flex_prob->num_servers; i++) {
        free(flex_prob->rigid_capacities[i]);
        free(flex_prob->fluid_capacities[i]);
    }
    free(flex_prob->rigid_capacities);
    free(flex_prob->fluid_capacities);
    for (i = 0; i < flex_prob->num_services; i++) {
        free(flex_prob->rigid_needs[i]);
        free(flex_prob->fluid_needs[i]);
    }
    free(flex_prob->rigid_needs);
    free(flex_prob->fluid_needs);
    free(flex_prob->slas);
    free(flex_prob);

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
