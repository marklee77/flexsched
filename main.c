#include "flexsched.h"

/* Global variable that keeps track of the current problem instance */
struct flexsched_instance INS;

/* This is a global variable that keeps track of server load
*     and that is used to speed up a bunch of the naively implemented
*        algorithms (especially greedy ones) */
float **global_server_fluid_loads;
float **global_server_rigid_loads;
float **global_server_fluidmin_loads;

/* Global for some VP stuff acceleration*/
float **global_vp_bin_loads;

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
    {"SLOWDIVING",   SLOWDIVING_scheduler, NULL, NULL, NULL},
    {"FASTDIVING",   FASTDIVING_scheduler, NULL, NULL, NULL},
    {"VP_PPMAXDIFF", VP_scheduler,         "PPMAXDIFF", NULL, NULL},
    {"VP_PPMAXRATIO",VP_scheduler,         "PPMAXRATIO", NULL, NULL},
    {"VP_PPSUM",     VP_scheduler,         "PPSUM", NULL, NULL},
    {"VP_PPMAX",     VP_scheduler,         "PPMAX", NULL, NULL},
    {"VP_CPMAXDIFF", VP_scheduler,         "CPMAXDIFF", NULL, NULL},
    {"VP_CPMAXRATIO",VP_scheduler,         "CPMAXRATIO", NULL, NULL},
    {"VP_CPSUM",     VP_scheduler,         "CPSUM", NULL, NULL},
    {"VP_CPMAX",     VP_scheduler,         "CPMAX", NULL, NULL},
    {"VP_BFDSUM",    VP_scheduler,         "BFDSUM", NULL, NULL},
    {"VP_BFDMAX",    VP_scheduler,         "BFDMAX", NULL, NULL},
    {"VP_BFDLEX",    VP_scheduler,         "BFDLEX", NULL, NULL},
    {"VP_FFDSUM",    VP_scheduler,         "FFDSUM", NULL, NULL},
    {"VP_FFDMAX",    VP_scheduler,         "FFDMAX", NULL, NULL},
    {"VP_FFDLEX",    VP_scheduler,         "FFDLEX", NULL, NULL},
    {"VP_CHEKURI",   VP_scheduler,         "CHEKURI", NULL, NULL},
    {"GA",           GA_scheduler,         "N", "N", "N"},
    {"GAF",          GA_scheduler,         "F", "N", "N"},
    {"GAFF",         GA_scheduler,         "F", "F", "N"},
    {"GAFFF",        GA_scheduler,         "F", "F", "F"},
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

  for (count=0, scheduler = strtok(s, sep); scheduler; 
                  count++,scheduler = strtok(NULL, sep)) {
    list = (char **)REALLOC(list,(count+1)*sizeof(char*));
    list[count] = strdup(scheduler);
  }
  list = (char **)REALLOC(list,(count+1)*sizeof(char*));
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
  if (!(f = fopen(file,"r"))) {
    return; 
  } 

  while (fgets(buffer, 1024, f)) {
    for (schedulerptr=list; *schedulerptr; schedulerptr++) {
      if (strcmp(*schedulerptr,"not_needed")) {
        if (strstr(buffer, *schedulerptr)) {
          fprintf(stderr,"Warning: no need to re-run algorithm '%s'\n",*schedulerptr);
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

    int i,j;
    int sched;
    struct timeval time1, time2, time3;
    double elasped_seconds;
    double non_optimized_average_yield;
    double non_optimized_utilization;

    if ((argc < 2) || (argc > 4)) {
        fprintf(stderr,"Usage: %s <scheduler list> [<instance file> [result file]]\n",argv[0]);
        fprintf(stderr,"  Known schedulers:  "); 
        sched = 0;
        while (implemented_schedulers[sched].name) {
          fprintf(stderr,"%s ",
              implemented_schedulers[sched].name);
          sched++; 
        }
        fprintf(stderr,
           "\n  See README file for scheduler descriptions\n");
	exit(1);
    } 

    // schedulers
    schedulers = parse_scheduler_list(argv[1]);

    // input
    if (argc < 3) {
	input = stdin;
    } else {
	if (!(input = fopen(argv[2],"r"))) {
	  fprintf(stderr,"Can't open file '%s' for reading\n",argv[2]);
	  exit(1);
   	} 
    }

    // output 
    if (argc < 4) {
	output = stdout;
    } else {
	remove_existing_schedulers(schedulers, argv[3]);
	if (!(output = fopen(argv[3],"a"))) {
	  fprintf(stderr,"Can't open file '%s' for writing/appending\n",argv[3]);
	  exit(1);
	}
    }

    /* Parse instance file */
    fscanf(input,"%d", &INS.numrigid);
    fscanf(input,"%d", &INS.numfluid);
    fscanf(input,"%d", &INS.numservers);
    fscanf(input,"%d", &INS.numservices);

    INS.rigidcapacities = (float **)calloc(INS.numservers, sizeof(float *));
    INS.fluidcapacities = (float **)calloc(INS.numservers, sizeof(float *));
    for (i=0; i < INS.numservers; i++) {
        INS.rigidcapacities[i] = (float *)calloc(INS.numrigid,sizeof(float *));
        INS.fluidcapacities[i] = (float *)calloc(INS.numfluid,sizeof(float *));
    }
    INS.slas = (float *)calloc(INS.numservices,sizeof(float));
    INS.rigidneeds = (float **)calloc(INS.numservices,sizeof(float*));
    INS.fluidneeds = (float **)calloc(INS.numservices,sizeof(float*));
    for (i=0; i<INS.numservices; i++) {
      INS.rigidneeds[i] = (float *)calloc(INS.numrigid,sizeof(float));
      INS.fluidneeds[i] = (float *)calloc(INS.numfluid,sizeof(float));
    }
    INS.mapping = (int *)calloc(INS.numservices,sizeof(int));
    INS.allocation = (float *)calloc(INS.numservices,sizeof(float));

    for (i = 0; i < INS.numservers; i++) {
        for (j=0; j < INS.numrigid; j++)
            fscanf(input,"%f", &(INS.rigidcapacities[i][j]));
        for (j=0; j < INS.numfluid; j++)
            fscanf(input,"%f", &(INS.fluidcapacities[i][j]));
    }
    for (i = 0; i < INS.numservices; i++) {
      fscanf(input,"%f", &(INS.slas[i]));
      for (j=0; j<INS.numrigid; j++) 
        fscanf(input,"%f", &(INS.rigidneeds[i][j]));
      for (j=0; j<INS.numfluid; j++) 
        fscanf(input,"%f", &(INS.fluidneeds[i][j]));
    }

    // Allocated memory for global_server_loads
    global_server_rigid_loads =	 	(float **)calloc(INS.numservers,sizeof(float*));
    global_server_fluid_loads = 	(float **)calloc(INS.numservers,sizeof(float*));
    global_server_fluidmin_loads = 	(float **)calloc(INS.numservers,sizeof(float*));
    for (i=0; i < INS.numservers; i++) {
      global_server_rigid_loads[i] = (float *)calloc(INS.numrigid,sizeof(float));
      global_server_fluid_loads[i] = (float *)calloc(INS.numfluid,sizeof(float));
      global_server_fluidmin_loads[i] = (float *)calloc(INS.numfluid,sizeof(float));
    }

    int status;

    for (schedulerptr = schedulers; *schedulerptr; schedulerptr++) {

      INS.misc_output[0] = '\0';

      // if not needed, then skip
      if (!strcmp(*schedulerptr,"not_needed"))  {
        continue;
      }

      // initialize mapping and allocations
      for (i=0; i < INS.numservices; i++) {
        INS.mapping[i] = -1;
        INS.allocation[i] = 0.0;
      }

      // search for the scheduler in the list
      sched = 0;
      while (implemented_schedulers[sched].name &&
         strcmp(*schedulerptr, implemented_schedulers[sched].name)) {
        sched++; 
      }
 
      // if we didn't find it, exit
      if (!implemented_schedulers[sched].name) {
        fprintf(stderr,
          "Error: no scheduler called \"%s\".\n", *schedulerptr);
	exit(1);
      }

      // Initialize the server loads 
      initialize_global_server_loads();

      // if we did find it, then run it
      gettimeofday(&time1, NULL);

      // Call the scheduler
      status = implemented_schedulers[sched].func(
                 implemented_schedulers[sched].arg1,
                 implemented_schedulers[sched].arg2,
                 implemented_schedulers[sched].arg3);

      gettimeofday(&time2, NULL);

      if ((status == RESOURCE_ALLOCATION_SUCCESS) && 
          (strcmp(*schedulerptr,"LPBOUND"))) {

        // Sanity check the allocation
        if (sanity_check() == RESOURCE_ALLOCATION_FAILURE) {
          fprintf(stderr,"Invalid allocation\n");
          exit(1);
        }

        // maximize_average_yield()
        non_optimized_average_yield = compute_average_yield();
        non_optimized_utilization = compute_utilization();
        maximize_average_yield();

        // Re-Sanity check the allocation
        if (sanity_check() == RESOURCE_ALLOCATION_FAILURE) {
          fprintf(stderr,"Invalid allocation\n");
          exit(1);
        }
      }

      // Compute elapsed time
      elasped_seconds  = (double)(time2.tv_sec - time1.tv_sec);
      elasped_seconds += (double)(time2.tv_usec - time1.tv_usec) / 1000000.0;
  
      // Print output
      fprintf(output,"%s|",*schedulerptr);
      if (status == RESOURCE_ALLOCATION_FAILURE) {
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
      } else if (!strcmp(*schedulerptr,"LPBOUND")) {
        fprintf(output,"%.3f|",INS.lpbound);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
        fprintf(output,"%.3f|",-1.0);
      } else {
        fprintf(output,"%.3f|",compute_minimum_yield());
        fprintf(output,"%.3f|",non_optimized_average_yield);
        fprintf(output,"%.3f|",compute_average_yield());
        fprintf(output,"%.3f|",non_optimized_utilization);
        fprintf(output,"%.3f|",compute_utilization());
      }
      fprintf(output,"%.3f|",elasped_seconds);
      fprintf(output,"%s\n",INS.misc_output);
    }
#if 0
    for (i = 0; i < INS.numservices; i++) {
        printf("%d\t", INS.mapping[i]);
    }
    printf("\n");
    for (i = 0; i < INS.numservices; i++) {
        printf("%f\t", INS.allocation[i]);
    }
    printf("\n");
#endif
//    fprintf(output,"##################\n");

    for (schedulerptr = schedulers; *schedulerptr; schedulerptr++) 
      free(*schedulerptr);
    free(schedulers);
    
    return 0;
}
