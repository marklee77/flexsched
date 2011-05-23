#include "flexsched.h"

flexsched_problem flex_prob = NULL;

struct scheduler_t implemented_schedulers[] = {
    {"GREEDY_S1_P1", GREEDY_scheduler,     "S1", "P1", NULL},
    {"GREEDY_S1_P2", GREEDY_scheduler,     "S1", "P2", NULL},
    {"GREEDY_S1_P3", GREEDY_scheduler,     "S1", "P3", NULL},
    {"GREEDY_S1_P4", GREEDY_scheduler,     "S1", "P4", NULL},
    {"GREEDY_S1_P5", GREEDY_scheduler,     "S1", "P5", NULL},
    {"GREEDY_S1_P6", GREEDY_scheduler,     "S1", "P6", NULL},
    {"GREEDY_S1_P7", GREEDY_scheduler,     "S1", "P7", NULL},
    {"GREEDY_S2_P1", GREEDY_scheduler,     "S2", "P1", NULL},
    {"GREEDY_S2_P2", GREEDY_scheduler,     "S2", "P2", NULL},
    {"GREEDY_S2_P3", GREEDY_scheduler,     "S2", "P3", NULL},
    {"GREEDY_S2_P4", GREEDY_scheduler,     "S2", "P4", NULL},
    {"GREEDY_S2_P5", GREEDY_scheduler,     "S2", "P5", NULL},
    {"GREEDY_S2_P6", GREEDY_scheduler,     "S2", "P6", NULL},
    {"GREEDY_S2_P7", GREEDY_scheduler,     "S2", "P7", NULL},
    {"GREEDY_S3_P1", GREEDY_scheduler,     "S3", "P1", NULL},
    {"GREEDY_S3_P2", GREEDY_scheduler,     "S3", "P2", NULL},
    {"GREEDY_S3_P3", GREEDY_scheduler,     "S3", "P3", NULL},
    {"GREEDY_S3_P4", GREEDY_scheduler,     "S3", "P4", NULL},
    {"GREEDY_S3_P5", GREEDY_scheduler,     "S3", "P5", NULL},
    {"GREEDY_S3_P6", GREEDY_scheduler,     "S3", "P6", NULL},
    {"GREEDY_S3_P7", GREEDY_scheduler,     "S3", "P7", NULL},
    {"GREEDY_S4_P1", GREEDY_scheduler,     "S4", "P1", NULL},
    {"GREEDY_S4_P2", GREEDY_scheduler,     "S4", "P2", NULL},
    {"GREEDY_S4_P3", GREEDY_scheduler,     "S4", "P3", NULL},
    {"GREEDY_S4_P4", GREEDY_scheduler,     "S4", "P4", NULL},
    {"GREEDY_S4_P5", GREEDY_scheduler,     "S4", "P5", NULL},
    {"GREEDY_S4_P6", GREEDY_scheduler,     "S4", "P6", NULL},
    {"GREEDY_S4_P7", GREEDY_scheduler,     "S4", "P7", NULL},
    {"GREEDY_S5_P1", GREEDY_scheduler,     "S5", "P1", NULL},
    {"GREEDY_S5_P2", GREEDY_scheduler,     "S5", "P2", NULL},
    {"GREEDY_S5_P3", GREEDY_scheduler,     "S5", "P3", NULL},
    {"GREEDY_S5_P4", GREEDY_scheduler,     "S5", "P4", NULL},
    {"GREEDY_S5_P5", GREEDY_scheduler,     "S5", "P5", NULL},
    {"GREEDY_S5_P6", GREEDY_scheduler,     "S5", "P6", NULL},
    {"GREEDY_S5_P7", GREEDY_scheduler,     "S5", "P7", NULL},
    {"GREEDY_S6_P1", GREEDY_scheduler,     "S6", "P1", NULL},
    {"GREEDY_S6_P2", GREEDY_scheduler,     "S6", "P2", NULL},
    {"GREEDY_S6_P3", GREEDY_scheduler,     "S6", "P3", NULL},
    {"GREEDY_S6_P4", GREEDY_scheduler,     "S6", "P4", NULL},
    {"GREEDY_S6_P5", GREEDY_scheduler,     "S6", "P5", NULL},
    {"GREEDY_S6_P6", GREEDY_scheduler,     "S6", "P6", NULL},
    {"GREEDY_S6_P7", GREEDY_scheduler,     "S6", "P7", NULL},
    {"GREEDY_S7_P1", GREEDY_scheduler,     "S7", "P1", NULL},
    {"GREEDY_S7_P2", GREEDY_scheduler,     "S7", "P2", NULL},
    {"GREEDY_S7_P3", GREEDY_scheduler,     "S7", "P3", NULL},
    {"GREEDY_S7_P4", GREEDY_scheduler,     "S7", "P4", NULL},
    {"GREEDY_S7_P5", GREEDY_scheduler,     "S7", "P5", NULL},
    {"GREEDY_S7_P6", GREEDY_scheduler,     "S7", "P6", NULL},
    {"GREEDY_S7_P7", GREEDY_scheduler,     "S7", "P7", NULL},
    {"METAGREEDY",   METAGREEDY_scheduler, NULL, NULL, NULL},
    {"METAGREEDYLIGHT",  METAGREEDYLIGHT_scheduler, NULL, NULL, NULL},
    {"MILP",         MILP_scheduler,       NULL, NULL, NULL},
    {"LPBOUND",      LPBOUND_scheduler,    NULL, NULL, NULL},
    {"RRND",         RRND_scheduler,       NULL, NULL, NULL},
    {"RRNZ",         RRNZ_scheduler,       NULL, NULL, NULL},
/*
    {"SLOWDIVING",   SLOWDIVING_scheduler, NULL, NULL, NULL},
    {"FASTDIVING",   FASTDIVING_scheduler, NULL, NULL, NULL},
*/
    {"VP_PPMAXDIFF", VP_scheduler,         "PPMAXDIFF", NULL, NULL},
    {"VP_PPMAXRATIO",VP_scheduler,         "PPMAXRATIO", NULL, NULL},
    {"VP_PPSUM",     VP_scheduler,         "PPSUM", NULL, NULL},
    {"VP_PPMAX",     VP_scheduler,         "PPMAX", NULL, NULL},
    {"VP_CPMAXDIFF", VP_scheduler,         "CPMAXDIFF", NULL, NULL},
    {"VP_CPMAXRATIO", VP_scheduler,         "CPMAXRATIO", NULL, NULL},
    {"VP_CPSUM",     VP_scheduler,         "CPSUM", NULL, NULL},
    {"VP_CPMAX",     VP_scheduler,         "CPMAX", NULL, NULL},
    {"VP_BFDSUM",    VP_scheduler,         "BFDSUM", NULL, NULL},
    {"VP_BFDMAX",    VP_scheduler,         "BFDMAX", NULL, NULL},
    {"VP_BFDLEX",    VP_scheduler,         "BFDLEX", NULL, NULL},
    {"VP_FFDSUM",    VP_scheduler,         "FFDSUM", NULL, NULL},
    {"VP_FFDMAX",    VP_scheduler,         "FFDMAX", NULL, NULL},
    {"VP_FFDLEX",    VP_scheduler,         "FFDLEX", NULL, NULL},

    {"HVP_PPMAX",     HVP_scheduler,        "PPMAX", NULL, NULL},
    {"HVP_CPMAX",     HVP_scheduler,        "CPMAX", NULL, NULL},
    {"HVP_CPMAXR",    HVP_scheduler,       "CPMAXR", NULL, NULL},
    {"HVP_CPMAXMAX",  HVP_scheduler,     "CPMAXMAX", NULL, NULL},
    {"HVP_CPMAXMAXR", HVP_scheduler,    "CPMAXMAXR", NULL, NULL},
    {"HVP_CPMAXMAXS", HVP_scheduler,    "CPMAXMAXS", NULL, NULL},
    {"HVP_CPMAXMAXRS",HVP_scheduler,   "CPMAXMAXRS", NULL, NULL},
/*
    {"VP_CHEKURI",   VP_scheduler,         "CHEKURI", NULL, NULL},
    {"GA",           GA_scheduler,         "N", "N", "N"},
    {"GAF",          GA_scheduler,         "F", "N", "N"},
    {"GAFF",         GA_scheduler,         "F", "F", "N"},
    {"GAFFF",        GA_scheduler,         "F", "F", "F"},
*/
    {NULL,           NULL,                 NULL,      NULL, NULL}
};

/* parse_scheduler_list() 
 *    Create a NULL-terminated list of scheduler names from
 *    the list command-line argument
 */
char **parse_scheduler_list(char *s) {
    char **list = NULL;
    char *sep=" \t";
    char *scheduler;
    int count;

    count = 0;
    for (scheduler = strtok(s, sep); scheduler; scheduler = strtok(NULL, sep)) {
        count++;
        list = (char **)REALLOC(list, count * sizeof(char*));
        list[count-1] = strdup(scheduler);
    }
    list = (char **)REALLOC(list, (count + 1) * sizeof(char*));
    list[count] = NULL;
    return list;
}

/* remove_existing_scheduler() 
 *    open the output file, if it exists, and remove schedulers for which 
 *    results are already available
 */
void remove_existing_schedulers(char **list, char *file)
{
    FILE *f;
    char buffer[1024];
    char **schedulerptr;

    // If not output file, then fine
    if (!(f = fopen(file,"r"))) return; 

    while (fgets(buffer, 1024, f)) {
        for (schedulerptr = list; *schedulerptr; schedulerptr++) {
            if (strcmp(*schedulerptr, "not_needed")) {
                if (strstr(buffer, *schedulerptr)) {
                    fprintf(stderr, 
                        "Warning: no need to re-run algorithm '%s'\n",
                        *schedulerptr);
                    free(*schedulerptr);
	                *schedulerptr = (char*)strdup("not_needed");
                }
            }
        }
    }

    fclose(f);

    return;
}

int main(int argc, char *argv[])
{
    char *scheduler;
    FILE *input, *output;
    char **schedulers;
    char **schedulerptr;

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
	    remove_existing_schedulers(schedulers, argv[3]);
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

    for (schedulerptr = schedulers; *schedulerptr; schedulerptr++) {

        // if not needed, then skip
        if (!strcmp(*schedulerptr,"not_needed")) continue;

        // search for the scheduler in the list
        for (i = 0; implemented_schedulers[i].name && 
            strcmp(*schedulerptr, implemented_schedulers[i].name); i++);
 
        // if we didn't find it, exit
        if (!implemented_schedulers[i].name) {
            fprintf(stderr, "Error: no scheduler called \"%s\".\n", 
                *schedulerptr);
            continue;
        }

        // if we did find it, then run it
        gettimeofday(&time1, NULL);

        // Call the scheduler
        flex_soln = implemented_schedulers[i].func(
            implemented_schedulers[i].arg1, implemented_schedulers[i].arg2, 
            implemented_schedulers[i].arg3);

        gettimeofday(&time2, NULL);

        if (flex_soln->success && strcmp(*schedulerptr, "LPBOUND")) {

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
        fprintf(output, "%s|", *schedulerptr);
        if (!strcmp(*schedulerptr,"LPBOUND")) {
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

    for (schedulerptr = schedulers; *schedulerptr; schedulerptr++) {
        free(*schedulerptr);
    }
    free(schedulers);

    return 0;
}
