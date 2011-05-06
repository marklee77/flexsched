#include <glpk.h>
#include "flexsched.h"


/* Global variable that keeps track of the current problem instance */
struct flexsched_instance INS;

/* This is a global variable that keeps track of server load
   and that is used to speed up a bunch of the naively implemented
   algorithms (especially greedy ones) */
float **global_server_fluid_loads;
float **global_server_rigid_loads;
float **global_server_fluidmin_loads;

/* Global for some VP stuff acceleration*/
float **global_vp_bin_loads;

/* Scheduler function prototypes */
int GREEDY_scheduler(char*,char*,char*);
int METAGREEDY_scheduler(char*,char*,char*);
int METAGREEDYLIGHT_scheduler(char*,char*,char*);
int MILP_scheduler(char*,char*,char*);
int LPBOUND_scheduler(char*,char*,char*);
int RRND_scheduler(char*,char*,char*);
int RRNZ_scheduler(char*,char*,char*);
int SLOWDIVING_scheduler(char*,char*,char*);
int FASTDIVING_scheduler(char*,char*,char*);
int VP_scheduler(char*,char*,char*);
int GA_scheduler(char*,char*,char*);

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




/* Local Prototype (sloppy for now) */
float compute_minimum_yield();
float compute_average_yield();
void maximize_average_yield();

/* Utility function */
void free_2D_int_array(int **a, int n)
{
  int i;
  for (i=0; i<n; i++)
    free(a[i]);
  free(a);
  return;
}

void initialize_global_server_loads()
{
  int i,j;

  for (i=0; i < INS.numservers; i++) {
    for (j=0; j < INS.numrigid; j++)  {
      global_server_rigid_loads[i][j] = 0.0;
    }
    for (j=0; j < INS.numfluid; j++)  {
      global_server_fluid_loads[i][j] = 0.0;
    }
    for (j=0; j < INS.numfluid; j++)  {
      global_server_fluidmin_loads[i][j] = 0.0;
    }
  }
}

/* Wrapperss around GLPK 
 * Returns 0 on success
 * Returns 1 if couldn't solve
 * Returns 2 if time limit was reached
 */
int solve_linear_program(int rational, int verbose, glp_prob *prob)
{
  int solver_status, solution_status;

  if (verbose) { // print out the program
    fprintf(stderr,"Writing LP to file 'linearprogram'..");
    glp_write_lp(prob,NULL,"linearprogram");
    fprintf(stderr,"done\n");
  }

  if (rational) {
    if (verbose) {fprintf(stderr,"Solving a rational LP ... ");}
    solver_status = glp_simplex(prob,NULL);
    solution_status = (glp_get_status(prob) != GLP_OPT);
    if (verbose) {fprintf(stderr,"solver_stat: %d, ",solver_status);}
    if (verbose) {fprintf(stderr,"solution_stat: %d\n",solution_status);}
  } else {
    if (verbose) {fprintf(stderr,"Solving a MILP ... "); }
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    parm.tm_lim = 10*60*1000; // 10 minutes time limit!

    solver_status = glp_intopt(prob,&parm);
    solution_status = (glp_mip_status(prob) != GLP_OPT);
    if (verbose) {fprintf(stderr,"solver_status = %d, ",solver_status);}
    if (verbose) {fprintf(stderr,"solution_status = %d\n",solution_status);}
  }
 
  if (verbose) { // print out the solution
    fprintf(stderr,"Writing LP solution to file 'solution'..");
//    lpx_print_mip(prob,"solution_other");
    glp_write_sol(prob,"solution");
    fprintf(stderr,"done\n");
  }

  // Reached the time limit
  if (solver_status == GLP_ETMLIM) {
    return 2;
  }

  // Failed for some other reason
  if ((solver_status) || (solution_status)) {
    return 1;
  }

  return 0;
}


int glp_stderr_out(void *info, const char *s) {
//    fprintf(stderr, s);
    return 1;
}

glp_prob *create_placement_lp(int rational)
{
    glp_prob *prob;
    int N = INS.numservices;
    int H = INS.numservers;
    int R = INS.numrigid;
    int F = INS.numfluid;
    int row, i, j, h, k;

    // Create GLPK problem
    glp_term_hook(&glp_stderr_out, NULL);
    prob = glp_create_prob();
    glp_set_prob_name(prob, "task placement");

    // Number of variables (called "columns"): 
    //  e_ih: 1 .. NH
    //  y_ih: NH+1 .. 2NH
    //  Y:    2NH+1
//    fprintf(stderr,"calling glp_add_cols with %d\n",2*N*H+1);
    glp_add_cols(prob, 2*N*H+1);

    // Defining convenient variables to reference 
    // columns without too many index calculations
    int **eih;
    eih = (int **)calloc(N+1,sizeof(int*));
    for (i=1; i<=N; i++)
      eih[i] = (int *)calloc(H+1,sizeof(int));
    k = 1;
    for (i=1; i<=N; i++) {
      for (h=1; h<=H; h++) {
        eih[i][h] = k++;
      }
    }
    int **yih;
    yih = (int **)calloc(N+1,sizeof(int*));
    for (i=1; i<=N; i++)
      yih[i] = (int *)calloc(H+1,sizeof(int));
    k = N*H+1;
    for (i=1; i<=N; i++) {
      for (h=1; h<=H; h++) {
        yih[i][h] = k++;
      }
    }
    int Y = 2*N*H+1;
    glp_set_col_name(prob, Y, "Y");

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_name(prob, Y, "Y");
    glp_set_col_bnds(prob, Y, GLP_DB, 0.0, 1.0);


    // Set bounds for the e_ij: binary if !rational
    k = 1;
    for (i=1; i<= N; i++) {
      for (j=1; j<= H; j++) {
         char buffer[10];
         sprintf(buffer,"e_{%d,%d}",i,j);
         if (!rational) {
           glp_set_col_kind(prob, k, GLP_BV);
         } else {
           glp_set_col_bnds(prob, k, GLP_DB, 0.0, 1.0);
         }
         k++;
      }
    }

    // Set bounds for the y_ij: positive
    k = N*H+1;
    for (i=1; i<= N; i++) {
      for (j=1; j<= H; j++) {
         char buffer[10];
         sprintf(buffer,"y_{%d,%d}",i,j);
         glp_set_col_name(prob, k, buffer);
         glp_set_col_bnds(prob, k, GLP_LO,0.0,-666);
         k++;
      }
    }

    
    // objective is minyield
    glp_set_obj_name(prob, "minimum scaled yield");
    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, Y, 1.0);

    // Number of rows  
    //  (A) for all i 	   sum_h e_ih = 1:       		  N rows
    //  (B) for all i,h	   y_ih <= e_ih				  NH rows
    //  (C) for all i      sum_h y_ih >= sla_i			  N rows
    //  (D) for all h, jr  sum_i r_ij*e_ih <= 1			  HR rows
    //  (E) for all h, jf  sum_i r_ij*y_ih <= 1			  HF rows
    //  (F) for all i	   (sum_h y_ih - sla_i)/(1-sla_i) >= Y	  N rows
    glp_add_rows(prob, N+N*H+N+H*(F+R)+N);
//    fprintf(stderr,"Adding %d rows to the problem\n",N+N*H+N+H*(F+R)+N);

    // compute the number of non-zero matrix elements
    int ne = N*H + N*H*2 + N*H + H*R*N + H*F*N + N*(H+1);
    int *ia = (int *)calloc(ne+1,sizeof(int));
    int *ja = (int *)calloc(ne+1,sizeof(int));
    double *ra = (double *)calloc(ne+1,sizeof(double));
  
    k = 1;
    // Constraint (A)
    //  (A) for all i 	   sum_h e_ih = 1:       	1 .. N
    row = 1;
//    fprintf(stderr,"Constraints (A): k = %d\n",k);
    for (i=1; i <= N; i++) {
      glp_set_row_bnds(prob,row,GLP_FX,1.0,1.0);
      for (h=1; h <= H; h++) {
        ia[k] = row; ja[k] = eih[i][h]; ra[k] = 1.0; k++;
      }
      row++;
    }
   
    // Constraint (A)
    // Constraint (B)
    //  (B) for all i,h	   y_ih <= e_ih			N+1 .. N(H+1)+1
    row = N+1;
//    fprintf(stderr,"Constraints (B): k = %d\n",k);
    for (i=1; i <= N; i++) {
      for (h=1; h <= H; h++) {
        glp_set_row_bnds(prob,row,GLP_UP,-666,0.0);
        ia[k] = row; ja[k] = yih[i][h]; ra[k] = 1.0;  k++;
        ia[k] = row; ja[k] = eih[i][h]; ra[k] = -1.0; k++;
        row++;
      }
    }


    // Constraint (C)
    //  (C) for all i      sum_h y_ih >= sla_i		N(H+1)+2 .. N(H+2)+2
    row = N + N*H + 1;
//    fprintf(stderr,"Constraints (C): k = %d\n",k);
    for (i=1; i<=N; i++) {
      glp_set_row_bnds(prob,row,GLP_LO,INS.slas[i-1],666);
      for (h=1; h <= H; h++) {
        ia[k] = row; ja[k] = yih[i][h]; ra[k] = 1.0; k++;
      }
      row++;
    }  

    // Constraint (D)
    //  (D) for all h, jr  sum_i r_ij*e_ih <= 1		N(H+2)+3 .. N(H+2)+3+HR
    row = N + N*H + N + 1;
//    fprintf(stderr,"Constraints (D): k = %d\n",k);
    for (h=1; h<=H; h++) {
      for (j=1; j <= R; j++) {
        glp_set_row_bnds(prob,row,GLP_UP,666,1.0);
        for (i=1; i <= N; i++) {
          ia[k] = row; ja[k] = eih[i][h]; ra[k] = INS.rigidneeds[i-1][j-1]; k++;
        } 
        row++;
      }
    }

    // Constraint (E)
    //  (E) for all h, jf  sum_i r_ij*y_ih <= 1		N(H+2)+4+HR .. N(H+2)+4+HR+HF
    row = N + N*H + N + H*R + 1;
//    fprintf(stderr,"Constraints (E): k = %d\n",k);
    for (h=1; h <= H; h++) {
      for (j=1; j <= F; j++) {
        glp_set_row_bnds(prob,row,GLP_UP,666,1.0);
        for (i=1; i <= N; i++) {
          ia[k] = row; ja[k] = yih[i][h]; ra[k] = INS.fluidneeds[i-1][j-1]; k++;
        } 
        row++;
      }
    }


    // Constraint (F)
    //  (F) for all i	   (sum_h y_ih - sla_i)/(1-sla_i) >= Y	
    //                   sum_h y_ih  + Y (sla_i - 1) >= sla_i
    //						N(H+2)+4+HR+HF+1 .. N(H+3)+4+HR+HF+1
    row = N + N*H + N + H*R + H*F + 1;
//    fprintf(stderr,"Constraints (F): k = %d\n",k);
    for (i=1; i <= N; i++) {
      glp_set_row_bnds(prob,row,GLP_LO,INS.slas[i-1],666);
      for (h = 1; h <= H; h++) {
        ia[k] = row; ja[k] = yih[i][h]; ra[k] = 1.0; k++;
      } 
      ia[k] = row; ja[k] = Y; ra[k] = INS.slas[i-1] - 1.0; k++;
      row++;
    }

    //fprintf(stderr,"k = %d\n",k-1);
    //fprintf(stderr,"ne = %d\n",ne);
    //fprintf(stderr,"row = %d\n",row-1);


    glp_load_matrix(prob, ne, ia, ja, ra);

    free_2D_int_array(eih, N+1);
    free_2D_int_array(yih, N+1);

    free(ia);
    free(ja);
    free(ra);
   
    return prob;
}

#if 1
double compute_LP_bound()
{
  int i,j; 
  float minratio, ratio, sum1, sum2;

  // Compute the bound directly
  minratio = -1.0;
  for (j=0; j<INS.numfluid; j++) {
    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<INS.numservices; i++) {
      sum1 += INS.slas[i] * INS.fluidneeds[i][j];
      sum2 += (1-INS.slas[i]) * INS.fluidneeds[i][j];
    }
    ratio = (INS.numservers - sum1) / sum2;
    if ((minratio == -1.0) || (ratio < minratio)) {
      minratio = ratio;
    }
  }
  if (minratio < 0.0)
    return -1.0;
  else 
    return MIN(1.0,minratio);
}
#endif


int LPBOUND_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  int i,j; 
  float minratio, ratio, sum1, sum2;

  // Compute the bound directly
  INS.lpbound = compute_LP_bound();

  // Validate with the solution of the Rational LP
  if (1) {
    int status;
    glp_prob *prob = create_placement_lp(RATIONAL);

    if (solve_linear_program(RATIONAL, NOT_VERBOSE, prob)) {
      return RESOURCE_ALLOCATION_FAILURE;
    } 

    if (fabs(INS.lpbound - glp_get_col_prim(prob, 2*INS.numservices*INS.numservers+1)) > EPSILON) {
      fprintf(stderr,"Error: Immediate bound = %.3f, LP bound = %.3f\n",
                     INS.lpbound,
                     glp_get_col_prim(prob, 2*INS.numservices*INS.numservers+1));
      exit(1);
    } else {
      fprintf(stderr,"AGREEMENT\n");
    }
  }

  return RESOURCE_ALLOCATION_SUCCESS;
}

int MILP_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  int i, h, k;
  int status;
  glp_prob *prob;

  // Create the problem (in MILP form)
  prob = create_placement_lp(INTEGER);

  // Solve it
  if (status = solve_linear_program(INTEGER, NOT_VERBOSE, prob)) {
    if (status == 2) {
      strcat(INS.misc_output,"T");
    } 
    return RESOURCE_ALLOCATION_FAILURE;
  } 

  // Retrieve the mapping
  k = 1;
  for (i = 0; i < INS.numservices; i++) {
    for (h = 0; h < INS.numservers; h++) {
      if (glp_mip_col_val(prob,k++) >= 1.0 - EPSILON) 
        INS.mapping[i] = h;
    }
  }

  // Retrieve the allocation
  for (i = 0; i < INS.numservices; i++) {
    for (h = 0; h < INS.numservers; h++) {
      if (INS.mapping[i] != h) {
        k++;
      } else {
        // y_ih in the linear program is the unscaled yield!
        INS.allocation[i] = (glp_mip_col_val(prob,k++) - INS.slas[i]) /
                        (1 - INS.slas[i]);
//        fprintf(stderr,"Mapping[%d] = %d\n", i, mapping[i]);
//        fprintf(stderr,"Allocation[%d] = %f\n", i, allocation[i]);
      }
    }
  }

  // For debugging, print the minimum yield achieved
//  fprintf(stderr,"MILP solution: minimum yield = %.3f\n",
//        glp_mip_col_val(prob,2*numservices*numservers+1));
  return RESOURCE_ALLOCATION_SUCCESS;
}

/* array_sum(): Utility function */
float array_sum(float *array, int size)
{
  int i;
  float sum = 0.0;
  for (i=0; i<size; i++)
    sum += array[i];
  return sum;
}

/* array_max(): Utility function */
float array_max(float *array, int size)
{
  int i;
  float max = array[0];
  for (i=1; i<size; i++)
    if (max < array[i])
      max = array[i];
  return max;
}

/* array_min(): Utility function */
float array_min(float *array, int size)
{
  int i;
  float min = array[0];
  for (i=1; i<size; i++)
    if (min > array[i])
      min = array[i];
  return min;
}

/* array_argmax(): Utility function */
int array_argmax(float *array, int size)
{
  int i;
  int argmax = 0;
  float max = array[0];
  for (i=1; i<size; i++)
    if (max < array[i]) {
      argmax = i;
      max = array[i];
    }
  return argmax;
}

/* Keeping the services in their original order */
int GREEDY_sort_S1(const void *x, const void *y)
{
  return 0;
}

/* Sorting the services by decreasing of max fluid needs */
int GREEDY_sort_S2(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float delta = array_max(INS.fluidneeds[ix],INS.numfluid) - 
                array_max(INS.fluidneeds[iy],INS.numfluid);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}

/* Sorting the services by decreasing sum of fluid needs */
int GREEDY_sort_S3(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float delta = array_sum(INS.fluidneeds[ix],INS.numfluid) - 
                array_sum(INS.fluidneeds[iy],INS.numfluid);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}

/* Sorting the services by decreasing max of rigid and constrained fluid needs */
int GREEDY_sort_S4(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float max_rigid_x;
  float max_rigid_y;
  float max_constrainedfluid_x;
  float max_constrainedfluid_y;

  max_rigid_x = array_max(INS.rigidneeds[ix],INS.numrigid);
  max_rigid_y = array_max(INS.rigidneeds[iy],INS.numrigid);

  max_constrainedfluid_x = INS.slas[ix] * array_max(INS.fluidneeds[ix],INS.numfluid);
  max_constrainedfluid_y = INS.slas[iy] * array_max(INS.fluidneeds[iy],INS.numfluid);

  float delta = MAX(max_rigid_x, max_constrainedfluid_x) -
                MAX(max_rigid_y, max_constrainedfluid_y);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}

/* Sorting the services by decreasing sum of rigid and constrained fluid needs */
int GREEDY_sort_S5(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float sum_rigid_x;
  float sum_rigid_y;
  float sum_constrainedfluid_x;
  float sum_constrainedfluid_y;

  sum_rigid_x = array_sum(INS.rigidneeds[ix],INS.numrigid);
  sum_rigid_y = array_sum(INS.rigidneeds[iy],INS.numrigid);

  sum_constrainedfluid_x = INS.slas[ix] * array_sum(INS.fluidneeds[ix],INS.numfluid);
  sum_constrainedfluid_y = INS.slas[iy] * array_sum(INS.fluidneeds[iy],INS.numfluid);

  float delta = (sum_rigid_x +  sum_constrainedfluid_x) -
                (sum_rigid_y +  sum_constrainedfluid_y);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}

/* Sorting the services by decreasing max of all needs */
int GREEDY_sort_S6(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float max_rigid_x;
  float max_rigid_y;
  float max_fluid_x;
  float max_fluid_y;

  max_rigid_x = array_max(INS.rigidneeds[ix],INS.numrigid);
  max_fluid_x = array_max(INS.fluidneeds[ix],INS.numfluid);
  max_rigid_y = array_max(INS.rigidneeds[iy],INS.numrigid);
  max_fluid_y = array_max(INS.fluidneeds[iy],INS.numfluid);

  float delta = MAX(max_rigid_x, max_fluid_x) -
                MAX(max_rigid_y, max_fluid_y);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}

/* Sorting the services by decreasing sum of all needs */
int GREEDY_sort_S7(const void *x, const void *y)
{
  int ix = *((int *)x); 
  int iy = *((int *)y); 

  float sum_rigid_x;
  float sum_rigid_y;
  float sum_fluid_x;
  float sum_fluid_y;

  sum_rigid_x = array_sum(INS.rigidneeds[ix],INS.numrigid);
  sum_fluid_x = array_sum(INS.fluidneeds[ix],INS.numfluid);
  sum_rigid_y = array_sum(INS.rigidneeds[iy],INS.numrigid);
  sum_fluid_y = array_sum(INS.fluidneeds[iy],INS.numfluid);

  float delta = (sum_rigid_x +  sum_fluid_x) -
                (sum_rigid_y +  sum_fluid_y);

  if (delta < 0)
    return 1;
  else if (delta > 0)
    return -1;
  else 
    return 0;
}


void GREEDY_sort_services(const char *S, int *sorted)
{
  int i;
  int (*compar)(const void *, const void *);

  // Select the appropriate sorting routine
  if (!strcmp(S,"S1")) {
    compar = GREEDY_sort_S1;
  } else if (!strcmp(S,"S2")) {
    compar = GREEDY_sort_S2;
  } else if (!strcmp(S,"S3")) {
    compar = GREEDY_sort_S3;
  } else if (!strcmp(S,"S4")) {
    compar = GREEDY_sort_S4;
  } else if (!strcmp(S,"S5")) {
    compar = GREEDY_sort_S5;
  } else if (!strcmp(S,"S6")) {
    compar = GREEDY_sort_S6;
  } else if (!strcmp(S,"S7")) {
    compar = GREEDY_sort_S7;
  } else {
    fprintf(stderr,"Greedy algorithm: unknown sorting procedure '%s'\n",S);
    exit(1);
  }

  // Initialize index array with the original order
  for (i=0; i<INS.numservices; i++)
    sorted[i]=i;

  // Sort the jobs
  qsort(sorted, INS.numservices, sizeof(int), compar);

  return;

}

/*
 * compute_sum_server_load:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_sum_server_load(
           int server, const char *type)
{
  int i; 

  float sumload = 0.0;

  for (i=0; i<INS.numservices; i++) {
    if (INS.mapping[i] != server)
      continue;
    if (!strcmp(type,"rigid")) {
      sumload += array_sum(INS.rigidneeds[i],INS.numrigid);
    } else if (!strcmp(type,"fluid")) {
      sumload += array_sum(INS.fluidneeds[i],INS.numfluid);
    } else if (!strcmp(type,"fluidmin")) {
      sumload += INS.slas[i] * array_sum(INS.fluidneeds[i],INS.numfluid);
    } else {
      fprintf(stderr,"compute_sum_server_load(): unknown type '%s'\n",type);
      exit(1);
    }
  }
  return sumload;
}

/*
 * compute_sum_server_load_fast:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_sum_server_load_fast(
           int server, const char *type)
{
  float sumload = 0.0;

  if (!strcmp(type,"rigid")) {
    sumload = array_sum(global_server_rigid_loads[server], INS.numrigid);
  } else if (!strcmp(type,"fluid")) {
    sumload = array_sum(global_server_fluid_loads[server], INS.numfluid);
  } else if (!strcmp(type,"fluidmin")) {
    sumload = array_sum(global_server_fluidmin_loads[server], INS.numfluid);
  } else {
    fprintf(stderr,"compute_sum_server_load(): unknown type '%s'\n",type);
    exit(1);
  }
  return sumload;
}

/*
 * compute_server_load_in_dimension:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_server_load_in_dimension(
           int server, const char *type, int dim)
{
  int i; 

  float load = 0.0;

  for (i=0; i<INS.numservices; i++) {
    if (INS.mapping[i] != server)
      continue;
    if (!strcmp(type,"rigid")) {
      load += INS.rigidneeds[i][dim];           
    } else if (!strcmp(type,"fluid")) {
      load += INS.fluidneeds[i][dim];           
    } else if (!strcmp(type,"fluidmin")) {
      load += INS.slas[i] * INS.fluidneeds[i][dim];           
    } else {
      fprintf(stderr,"compute_server_load_in_dimension: unknown type '%s'\n",type);
      exit(1);
    }
  }
  return load;
}

/*
 * compute_server_load_in_dimension_fast:  (used for GREEDY)
 *  type: "fluid", "rigid", "fluidmin"
 */
float compute_server_load_in_dimension_fast(
           int server, const char *type, int dim)
{
  int i; 

  float load;

  if (!strcmp(type,"rigid")) {
    load = global_server_rigid_loads[server][dim];
  } else if (!strcmp(type,"fluid")) {
    load = global_server_fluid_loads[server][dim];
  } else if (!strcmp(type,"fluidmin")) {
    load = global_server_fluidmin_loads[server][dim];
  } else {
    fprintf(stderr,"compute_server_load_in_dimension: unknown type '%s'\n",type);
    exit(1);
  }
  return load;
}

int service_can_fit_on_server(int service, int server)
{
  int i,j;
  

  // Rigid needs
  for (j=0; j < INS.numrigid; j++) {
    float load;
    load = 0.0;
    for (i=0; i<INS.numservices; i++) {
      if (i == service)
        continue;
      if (INS.mapping[i] != server)
        continue;
      load += INS.rigidneeds[i][j];
    } 
    // BEWARE OF THE EPSILON
    if (1.0 - load  + EPSILON < INS.rigidneeds[service][j]) {
      return 0;
    }
  }
  // Fluid needs
  for (j=0; j < INS.numfluid; j++) {
    float load;
    load = 0.0;
    for (i=0; i<INS.numservices; i++) {
      if (i == service)
        continue;
      if (INS.mapping[i] != server)
        continue;
      load += INS.slas[i] * INS.fluidneeds[i][j];
    } 
    if (1.0 - load + EPSILON < INS.slas[service] * INS.fluidneeds[service][j])
      return 0;
  }

  return 1;
}

int service_can_fit_on_server_fast(int service, int server)
{
  int i,j;
  
  // Rigid needs
  for (j=0; j < INS.numrigid; j++) {
    // BEWARE OF THE EPSILON
    if (1.0 - global_server_rigid_loads[server][j] + EPSILON < INS.rigidneeds[service][j]) {
      return 0;
    }
  }
  // Rigid needs
  for (j=0; j < INS.numfluid; j++) {
    // BEWARE OF THE EPSILON
    if (1.0 - global_server_fluidmin_loads[server][j] + EPSILON < INS.slas[service] * INS.fluidneeds[service][j]) {
      return 0;
    }
  }

  return 1;
}

/* TO UNDO THE FAST THING: just remove all the _fast below */

int GREEDY_pick_server_P1(int service)
{
  int i;
  float load;
  float minload;
  int picked;
  
  minload = -1.0;
  picked = -1;
  for (i=0; i<INS.numservers; i++) {
    if (!service_can_fit_on_server_fast(service,i)) {
      continue;
    }
    load = compute_server_load_in_dimension_fast(
               i, "fluid",array_argmax(INS.fluidneeds[service],INS.numfluid)); 
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
  for (i=0; i<INS.numservers; i++) {
    if (!service_can_fit_on_server_fast(service,i)) {
      continue;
    }
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
  int i,picked;
  float objload;

  picked = -1;
  //fprintf(stderr,"PLACING Service %d\n",service);
  for (i=0; i<INS.numservers; i++) {
    float load;
    //fprintf(stderr,"Server %d: ",i);
    if (!service_can_fit_on_server_fast(service,i)) {
        //fprintf(stderr,"Can't fit\n");
      continue;
   }
    if (array_max(INS.rigidneeds[service],INS.numrigid) > 
           INS.slas[service]*array_max(INS.fluidneeds[service],INS.numfluid)) {
      load = compute_server_load_in_dimension_fast(i, "rigid",
                           array_argmax(INS.rigidneeds[service],INS.numrigid));
    } else {
      load = compute_server_load_in_dimension_fast(i, "fluidmin",
                           array_argmax(INS.fluidneeds[service],INS.numfluid));
    }
    //fprintf(stderr,"load(%d) = %.2f\n",i,load);
    if (!strcmp(mode,"firstfit")) {
       picked = i;
       //fprintf(stderr,"  firstfit picked %d\n",picked);
       break;
    } else if (!strcmp(mode,"bestfit")) {
      if ((picked == -1) || (load > objload)) {
        objload = load;
        picked = i;
        //fprintf(stderr,"  bestfit picked %d\n",picked);
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

int GREEDY_pick_server_P4_P6(int service, const char *mode)
{
  int i,picked;
  float objload = 0.0;

  picked = -1;
  //fprintf(stderr,"PLACING Service %d\n",service);
  for (i=0; i<INS.numservers; i++) {
    float load;
    //fprintf(stderr,"Server %d: ",i);
    if (!service_can_fit_on_server_fast(service,i)) {
      //fprintf(stderr,"Can't fit\n");
      continue;
    }
    load = compute_sum_server_load_fast(i, "rigid") +
           compute_sum_server_load_fast(i, "fluidmin");

    //fprintf(stderr,"load(%d) = %.2f\n",i,load);

    if (!strcmp(mode,"firstfit")) {
       picked = i;
       //fprintf(stderr,"  firstfit picked %d\n",picked);
       break;
    } else if (!strcmp(mode,"bestfit")) {
      if ((picked == -1) || (load > objload)) {
        objload = load;
        picked = i;
        //fprintf(stderr,"  bestfit picked %d\n",picked);
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

int GREEDY_pick_server_P3(int service) 
{
  return GREEDY_pick_server_P3_P5(service, "bestfit");
}

int GREEDY_pick_server_P4(int service) 
{
  return GREEDY_pick_server_P4_P6(service, "bestfit");
}

int GREEDY_pick_server_P5(int service) 
{
  return GREEDY_pick_server_P3_P5(service, "worstfit");
}

int GREEDY_pick_server_P6(int service) 
{
  return GREEDY_pick_server_P4_P6(service, "worstfit");
}

int GREEDY_pick_server_P7(int service) 
{
  int i,picked;
  float objload = 0.0;

  picked = -1;
  for (i=0; i<INS.numservers; i++) {
    float load;
    if (!service_can_fit_on_server_fast(service,i)) {
      continue;
    }
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
  } else {
    fprintf(stderr,"Greedy algorithm: unknown picking procedure '%s'\n",P);
    exit(1);
  }
}

int GREEDY_compute_mapping(const char *P, int *sorted)
{
  int i,j;

  /* Initialize the server load (for speed of computation) */
  initialize_global_server_loads();

  /* Deal with the services */
  for (i=0; i<INS.numservices; i++) {
    INS.mapping[sorted[i]] = GREEDY_pick_server(P, sorted[i]);
    if (INS.mapping[sorted[i]] == -1) {
      return RESOURCE_ALLOCATION_FAILURE;
    }
    /* Update the server loads */
    for (j=0; j < INS.numrigid; j++) {
      global_server_rigid_loads[INS.mapping[sorted[i]]][j] += INS.rigidneeds[sorted[i]][j];
    }
    for (j=0; j < INS.numfluid; j++) {
      global_server_fluid_loads[INS.mapping[sorted[i]]][j] += INS.fluidneeds[sorted[i]][j];
    }
    for (j=0; j < INS.numfluid; j++) {
      global_server_fluidmin_loads[INS.mapping[sorted[i]]][j] += INS.slas[sorted[i]]  * INS.fluidneeds[sorted[i]][j];
    }
  }
  return RESOURCE_ALLOCATION_SUCCESS;
}

/*
 * Used for GREEDY algorithms
 * Right now, gives every services on the same server
 * the same yield, which is the highest feasible such yield.
 */
void compute_allocations_given_mapping(int server)
{
  int i,j;
  float yield, yield_min;

  yield_min = -1.0;
  for (j=0; j<INS.numfluid; j++) {
    float sum1=0.0, sum2=0.0;
    for (i=0; i<INS.numservices; i++) {
      if (INS.mapping[i] != server)
        continue;
      sum1 += INS.slas[i] * INS.fluidneeds[i][j];
      sum2 += (1 - INS.slas[i]) * INS.fluidneeds[i][j];
    }
    yield = MIN(1,(1 - sum1) / sum2);
    if ((yield_min == -1.0) || (yield < yield_min))
      yield_min = yield;
  }

  for (i=0; i<INS.numservices; i++) {
    if (INS.mapping[i] != server)
      continue;
    INS.allocation[i] = yield_min;
  }
  return;
}

int GREEDY_scheduler(char *S, char *P, char *ignore)
{
  int i;
  int sorted[INS.numservices];

  // Sort the services appropriately
  GREEDY_sort_services(S,sorted);
#if 0
  fprintf(stderr,"SORTED SERVICES:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"%d, ",sorted[i]);
  }
  fprintf(stderr,"\n");
#endif


  // Map each service to a particular host
  if (GREEDY_compute_mapping(P,sorted) == RESOURCE_ALLOCATION_FAILURE) {
    return RESOURCE_ALLOCATION_FAILURE; 
  }
#if 0
  fprintf(stderr,"MAPPING:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d on server %d\n",i,INS.mapping[i]);
  }
#endif

  // For each host, compute its services' allocations
  for (i=0; i<INS.numservers; i++) {
    compute_allocations_given_mapping(i);
  } 
#if 0
  fprintf(stderr,"ALLOCATIONS:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d has an allocation (yield) of %f\n",i,INS.allocation[i]);
  }

#endif
  
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

        if ( (maxminyield == -1) || 
             (minyield > maxminyield) ||
             ((minyield > maxminyield - EPSILON) && (aveyield > maxaveyield)) ) {
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
  if (mapping)
    free(mapping);
  if (allocation)
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
  if (mapping)
    free(mapping);
  if (allocation)
    free(allocation);

  return metastatus;
}

/*
 *  Wrapper for RRND and RRNZ
 *  	- RRND: LPROUNDING(0.0)
 *  	- RRNZ: LPROUNDING(xxx)
 */
int LPROUNDING_scheduler(float epsilon)
{
  glp_prob *prob;
  int status;
  int i,h,k;
  float rational_mapping[INS.numservices][INS.numservers];

  srand(RANDOM_SEED);

  // Create the placement problem in rational mode
  // and solve it to compute rational mappings,
  // adding an epsilon to zero values
  prob = create_placement_lp(RATIONAL);
  if (solve_linear_program(RATIONAL, NOT_VERBOSE, prob)) {
    return RESOURCE_ALLOCATION_FAILURE;
  } 

  k = 1;
  for (i = 0; i < INS.numservices; i++) {
    for (h = 0; h < INS.numservers; h++) {
      rational_mapping[i][h] =  glp_get_col_prim(prob,k++);
      if (rational_mapping[i][h] <= 0.0)
        rational_mapping[i][h] = epsilon;
    }
  }

  // For each service, pick on which server it lands
  for (i=0; i < INS.numservices; i++) {

    int feasible[INS.numservers];
    int num_feasible_servers = 0;

    for (h=0; h < INS.numservers; h++)
      feasible[h] = 0;

    // Mark hosts that can't work and count how
    // Many are possible
    for (h=0; h < INS.numservers; h++) {
      if (service_can_fit_on_server(i,h) && 
             (rational_mapping[i][h] > 0.0)) {
        feasible[h] = 1;
        num_feasible_servers++;
      }
    }
 
    // If nobody works, forget it
    if (num_feasible_servers == 0)
      return RESOURCE_ALLOCATION_FAILURE;

    // Compute the sum of the probabilities of the feasible servers
    float total_alloc = 0.0;
    for (h=0; h<INS.numservers; h++) {
      if (feasible[h]) {
        total_alloc += rational_mapping[i][h];
      }
    }

    // Pick a probability
    float rnd = total_alloc*(rand() / (RAND_MAX + 1.0));
    int picked_server = -1;
    float total = 0.0;
    for (h=0; h < INS.numservers; h++) {
      // skip unfeasible servers
      if (!feasible[h])
        continue;
      total += rational_mapping[i][h];
      if (total >= rnd) {
        picked_server = h;
        break;
      }
    } 
    
    // set the mappings appropriately
     INS.mapping[i] = picked_server;
  }


  // Compute the allocations
  for (h=0; h<INS.numservers; h++) {
    compute_allocations_given_mapping(h);
  }

  return RESOURCE_ALLOCATION_SUCCESS;
}

int RRND_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  return LPROUNDING_scheduler(0.0);
}

int RRNZ_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  return LPROUNDING_scheduler(0.01);
}


int DIVING_scheduler(int mode)
{
  int i,h;
  glp_prob *prob;
  int status;

  // Create the placement problem in rational mode
  prob = create_placement_lp(RATIONAL);

  // Create the convenient index arrays
  int **eih = (int **)calloc(INS.numservices+1,sizeof(int*));
  for (i=1; i<=INS.numservices; i++)
    eih[i] = (int *)calloc(INS.numservers+1,sizeof(int));
  int k = 1;
  for (i=1; i<=INS.numservices; i++) {
    for (h=1; h<=INS.numservers; h++) {
      eih[i][h] = k++;
    }
  }

  int iter;
  int num_iterations;
  
  if (mode == SLOW_DIVE)
    num_iterations = INS.numservers*INS.numservices;
  else
    num_iterations = INS.numservices;

  for (iter = 0; iter < num_iterations; iter++) {

    // fprintf(stderr,"diving iteration %d\n",iter);
    // Solve the problem
    if (solve_linear_program(RATIONAL, NOT_VERBOSE, prob)) {
      free_2D_int_array(eih,INS.numservices+1);
      return RESOURCE_ALLOCATION_FAILURE;
    } 

    // Find the closest-to-integer
    float mindelta = -1.0;
    int i_selected, h_selected;
    float new_value;
    for (i=1; i <= INS.numservices; i++) {
      for (h=1; h <= INS.numservers; h++) {
        // skip over things that were already fixed
        if (glp_get_col_type(prob,eih[i][h]) == GLP_FX)
          continue;
        // compute the delta and update mindelta
        float value = glp_get_col_prim(prob,eih[i][h]);
        float delta;
        if (mode == SLOW_DIVE)       
          delta = MIN(1.0 - value, value);
        else
          delta = 1.0 - value;
//        fprintf(stderr,"e_{%d,%d} = %.3f\n",i,h,value);
        if ((mindelta == -1.0) || (delta < mindelta)) {
          mindelta = delta;
          i_selected = i;
          h_selected = h; 
          if (mode == SLOW_DIVE)
            new_value = (1.0 - value > value ? 0.0 : 1.0);
          else
            new_value = 1.0;
        }
      }
    }

    // Modify the linear program accordingly
//    fprintf(stderr,"Fixing e_{%d,%d} to %.2f\n",
//                i_selected, h_selected, new_value);
    glp_set_col_bnds(prob, eih[i_selected][h_selected], 
                     GLP_FX, new_value, new_value);

    
    // In fast dive, we always just set something to 1.0
    // In slow dive we could also do some optimization, but
    // it's a bit of a pain due to bookkeeping (remembering
    // what has been fixed in the past, etc.)
    if (mode == FAST_DIVE) {
      for (h = 1; h <= INS.numservers; h++) {
        if (h != h_selected) {
          glp_set_col_bnds(prob, eih[i_selected][h], GLP_FX, 0.0, 0.0);
        }
      }
    }
  }

// Solve the program one last time, for the final allocation
  if (solve_linear_program(RATIONAL, NOT_VERBOSE, prob)) {
    free_2D_int_array(eih,INS.numservices+1);
    return RESOURCE_ALLOCATION_FAILURE;
  } 

//  for (i=1; i <= numservices; i++) {
//    for (h=1; h <= numservers; h++) {
//      fprintf(stderr,"e_{%d,%d} = %.3f\n",
//              i,h,glp_get_col_prim(prob, eih[i][h]));
//    }
//  }

  // Create the mapping
//  fprintf(stderr,"Creating the mapping\n");
  for (i=1; i <= INS.numservices; i++) {
    for (h=1; h <= INS.numservers; h++) {
      if (glp_get_col_prim(prob,eih[i][h]) >= 1.0-EPSILON) {
        INS.mapping[i-1] = h-1;
      }
    }
  }

  // Compute the allocations
//  fprintf(stderr,"Computing the allocations\n");
  for (h=0; h<INS.numservers; h++) {
    compute_allocations_given_mapping(h);
  }

  free_2D_int_array(eih,INS.numservices+1);
  return RESOURCE_ALLOCATION_SUCCESS;
}

int SLOWDIVING_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  return DIVING_scheduler(SLOW_DIVE);
}

int FASTDIVING_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
  return DIVING_scheduler(FAST_DIVE);
}

/* generate_vp_instance() */
struct vp_instance *generate_vp_instance(float yield) 
{
  int i,j,dim_counter;
  struct vp_instance *instance;

  // INPUT
  instance = (struct vp_instance *)calloc(1, sizeof(struct vp_instance));
  instance->num_vectors = INS.numservices;
  instance->num_dims = INS.numrigid + INS.numfluid;
  instance->vectors = (struct vp_vector *)calloc(
                   instance->num_vectors,sizeof(struct vp_vector));
  instance->mapping = (int *)calloc(instance->num_vectors,sizeof(int));
  for (i=0; i<instance->num_vectors; i++) {
    instance->vectors[i].service = i;
    instance->vectors[i].num_dims = INS.numrigid + INS.numfluid;
    instance->vectors[i].x = (float *)calloc(
             instance->vectors[i].num_dims, sizeof(float));
    dim_counter = 0;
    for (j=0; j < INS.numrigid; dim_counter++, j++) {
      instance->vectors[i].x[dim_counter] = INS.rigidneeds[i][j]; 
    }
    for (j=0; j < INS.numfluid; dim_counter++, j++) {
      instance->vectors[i].x[dim_counter] = 
              ((1.0 - INS.slas[i])*yield + INS.slas[i]) * INS.fluidneeds[i][j]; 
    }
  }

  // OUTPUT
  instance->num_bins = 0;
  for (i=0; i < instance->num_vectors; i++)
    instance->mapping[i] = -1; // initialized to "not mapped"

#if 0
  fprintf(stderr,"VP instance:\n");
  for (i=0; i<instance->num_vectors; i++) {
    int j;
    for(j=0; j<instance->num_dims; j++) 
      fprintf(stderr,"%.2f ",instance->vectors[i].x[j]);
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"----------------\n");
#endif

  return instance;
}

/* free_vp_instance */
void free_vp_instance(struct vp_instance *instance)
{
  int i;
  for (i=0; i<instance->num_vectors; i++) {
    free(instance->vectors[i].x);
  }
  free(instance->vectors);
  free(instance->mapping);
  free(instance);
  return;
}

/* Lexicographical comparison of vectors */
int VP_sort_LEX(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
  
  for (i=0; i < vx->num_dims; i++) {
    if (vx->x[i] < vy->x[i])
      return 1;
    if (vx->x[i] > vy->x[i])
      return -1;
  }
  return 0;
}

/* Max comparison of vectors */
int VP_sort_MAX(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
  float max_x, max_y;
  
  max_x = vx->x[0];
  max_y = vy->x[0];

  for (i=0; i < vx->num_dims; i++) {
    if (vx->x[i] > max_x)
      max_x = vx->x[i];
    if (vy->x[i] > max_y)
      max_y = vy->x[i];
  }
  if (max_x > max_y)
    return -1;
  if (max_x < max_y)
    return 1;
  return 0;
}

/* Sum comparison of vectors */
int VP_sort_SUM(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
  float sum_x, sum_y;
  
  sum_x = 0.0;
  sum_y = 0.0;

  for (i=0; i < vx->num_dims; i++) {
    sum_x += vx->x[i]; 
    sum_y += vy->x[i]; 
  }
 
  if (sum_x > sum_y)
    return -1;
  if (sum_x < sum_y)
    return 1;
  return 0;
}

/* Misc comparison of vectors  */
/* By decreasing order of MISC */
int VP_sort_MISC(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
 
  if (vx->misc > vy->misc)
    return -1;
  if (vx->misc < vy->misc)
    return 1;
  return 0;
}


int vp_vector_can_fit_in_bin(struct vp_instance *vp, int v, int b) 
{
  int i,j;
  float load[vp->num_dims]; 

  for (j=0; j < vp->num_dims; j++) 
    load[j] = 0.0;

  // Compute the load
  for (i=0; i < vp->num_vectors; i++) {
    if (vp->mapping[i] != b)
      continue;
    for (j=0; j < vp->num_dims; j++) {
      load[j] += vp->vectors[i].x[j];
    }
  }

  // Check that the load can accomodate the new vector
  for (j=0; j < vp->num_dims; j++) {
    if (load[j] + vp->vectors[v].x[j] > 1.0)
      return 0;
  }
  return 1;
}

// load computed as the sum of all dims of all objects
int vp_vector_can_fit_in_bin_fast(struct vp_instance *vp, int v, int b) 
{
  int i,j;
  float load[vp->num_dims]; 

  for (j=0; j < vp->num_dims; j++)  {
    if (global_vp_bin_loads[b][j] + vp->vectors[v].x[j] > 1.0) {
      return 0;
    }
  }
  return 1;
}

// load computed as the sum of all dims of all objects
float vp_compute_bin_load(struct vp_instance *vp, int b)
{
  float load = 0.0;
  int i,j;

  for (i=0; i < vp->num_vectors; i++) {
    if (vp->mapping[i] != b)
      continue;
    for (j=0; j < vp->num_dims; j++)
      load += vp->vectors[i].x[j];
  }
  return load;
}

// load computed as the sum of all dims of all objects
float vp_compute_bin_load_fast(struct vp_instance *vp, int b)
{
  float load = 0.0;
  int j;

  for (j=0; j < vp->num_dims; j++) {
    load += global_vp_bin_loads[b][j];
  }
  return load;
}

/* A helper function to compute the thingies from
 * the Maruyama article, and puts them into the misc 
 * field of the vp->vectors
 */
void vp_compute_degrees_of_dominance(struct vp_instance *vp)
{
  int k,j,i, status;
  float *sums, *degrees;

  degrees = (float *)calloc(vp->num_dims, sizeof(float));
  sums = (float *)calloc(vp->num_bins, sizeof(float));

  // For each dimension compute its degree of dominance
  for (j=0; j < vp->num_dims; j++) {
    // sets the sums to zero
    for (k=0; k < vp->num_bins; k++)
      sums[k] = 0.0;
    // compute the sums
    for (i=0; i < vp->num_vectors; i++) {
      if (vp->mapping[i] < 0) {
        free(degrees); 
        free(sums);
        fprintf(stderr,"Warning: can't compute degrees of dominance due to unmapped vector!\n");
        return;
      }
      sums[vp->mapping[i]] += vp->vectors[i].x[j];
    }
    // compute the ceil
    for (k=0; k < vp->num_bins; k++)
      sums[k] = ceilf(sums[k] - 1.0);
    // compute the degree
    degrees[j] = 0.0;
    for (k=0; k < vp->num_bins; k++)
      degrees[j] += sums[k];
    degrees[j] /= vp->num_bins;
  }

  for (i=0; i < vp->num_vectors; i++) {
    vp->vectors[i].misc = 0.0;   
    for (j=0; j < vp->num_dims; j++)
      vp->vectors[i].misc += degrees[j] * vp->vectors[i].x[j];
  }

  free(sums);
  free(degrees);
  return;
}

/* Implementation of The standard "Fit" algorithms for
 * solving binpacking:
 *    "FIRST"   or "BEST":  fitting policy
 *    "LEX", "MAX", or "SUM": sorting policy
 *    "MARUYAMA" sorting: the Maruyama optimization
 *
 * Returns the number of used bins
 */
//#define FIT_DEBUG 0
int solve_vp_instance_FITD(const char *fit_type, 
                           const char *sort_type,
                           struct vp_instance *vp)
{
  int i,j,k;
  int (*compar)(const void *, const void *);
  int max_bin;

  // Sort the vectors in the instance according to the sort type
  if (!strcmp(sort_type,"LEX")) {
    compar = VP_sort_LEX;
  } else if (!strcmp(sort_type,"MAX")) {
    compar = VP_sort_MAX;
  } else if (!strcmp(sort_type,"SUM")) {
    compar = VP_sort_SUM;
  } else if (!strcmp(sort_type,"MARUYAMA")) { 
    compar = VP_sort_MISC;
  } else {
    fprintf(stderr,"Invalid VP sort type '%s'\n",sort_type);
    exit(1);
  }
  qsort(vp->vectors, vp->num_vectors, 
        sizeof(struct vp_vector), compar);

#if 0
  fprintf(stderr,"SORTED VP instance:\n");
  for (i=0; i<vp->num_vectors; i++) {
    int j;
    for(j=0; j<vp->num_dims; j++) 
      fprintf(stderr,"%.2f ",vp->vectors[i].x[j]);
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"----------------\n");
#endif

  // Allocate and initialize vp_bin_loads 
  global_vp_bin_loads = (float **)calloc(vp->num_vectors,sizeof(float*));
  for (k=0; k < vp->num_vectors; k++) {
    global_vp_bin_loads[k] = (float *)calloc(vp->num_dims, sizeof(float));
  }

  max_bin = 0;

  // Place vectors into bins
  for (i=0; i < vp->num_vectors; i++) {
    if (!strcmp(fit_type, "FIRST")) { // First Fit
      for (j=0; j<vp->num_vectors; j++) { // At worst, one bin per job
        if (vp_vector_can_fit_in_bin_fast(vp, i, j)) {
          vp->mapping[i] = j;        
          for (k=0; k < vp->num_dims; k++) {
            global_vp_bin_loads[j][k] += vp->vectors[i].x[k];
          }
	  if (j > max_bin)
            max_bin = j;
          break;
        } 
      }
    } else if (!strcmp(fit_type, "BEST")) { // Best Fit
      float load, max_load = -1.0;
#ifdef FIT_DEBUG
      fprintf(stderr,"Trying to fit vector ");
      for (j=0; j<vp->num_dims; j++)
        fprintf(stderr,"%.2f ",vp->vectors[i].x[j]);
      fprintf(stderr,"(svc=%d) into a bin\n",i);
#endif
      for (j=0; j<vp->num_vectors; j++) { // At worst one bin per vector
        if (!vp_vector_can_fit_in_bin_fast(vp, i, j)) {
#ifdef FIT_DEBUG
          fprintf(stderr,"  can't fit in bin %d\n",j);
#endif
          continue;
        }
        load = vp_compute_bin_load_fast(vp, j);
#ifdef FIT_DEBUG
        fprintf(stderr,"  bin %d has load %.3f\n",j,load);
#endif
        if ((max_load == -1) || (load > max_load)) {
          max_load = load;
          vp->mapping[i] = j;
          for (k=0; k < vp->num_dims; k++) {
            global_vp_bin_loads[j][k] += vp->vectors[i].x[k];
          }
	  if (j > max_bin)
            max_bin = j;
        }
      } 
    } else {
      fprintf(stderr,"Invalid VP fit type '%s'\n",fit_type);
      exit(1);
    }
#ifdef FIT_DEBUG
    fprintf(stderr,"vector %d (%.3f %.3f) was put in bin %d\n",
           i, 
           vp->vectors[i].x[0],
           vp->vectors[i].x[1],
           vp->mapping[i]);
#endif
  }
  vp->num_bins = (max_bin+1);

  for (k=0; k < vp->num_vectors; k++) {
    free(global_vp_bin_loads[k]);
  }
  free(global_vp_bin_loads);

  // Apply the Marumaya optimization in case it helps, if
  // we didn't just do it
  if (vp->num_bins > INS.numservers) {
    if (strcmp(sort_type,"MARUYAMA")) {
      int status;
      int old_num_bins = vp->num_bins;
      vp_compute_degrees_of_dominance(vp);
      status = solve_vp_instance_FITD(fit_type,"MARUYAMA",vp);
      if (vp->num_bins <= INS.numservers) {
        strcat(INS.misc_output,"M");
        // fprintf(stderr,"MARUYAMA HELPED! (%d -> %d)\n", old_num_bins, vp->num_bins);
      }
      return status;
    }
  }

  return VP_SOLVED;
}

void find_top_two_dims(float *array, int n,
                       int *d1, int *d2) 
{
  int i;
  float max2;
  
  *d1 = array_argmax(array,n);
  max2 = -1.0;
  for (i=0; i < n; i++) {
    if (i == *d1)
      continue; 
    if ((max2 == -1.0) || (array[i] > max2)) {
      max2 = array[i];
      *d2 = i;
    }
  }
 
  return;
}

/* MCB_compar_MAX */
int MCB_compar_MAX(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);

  if (max1 > max2) {
    return -1;
  } else if (max2 > max1) {
    return 1;
  } else {
    return 0;
  }
}

/* MCB_compar_SUM */
int MCB_compar_SUM(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float sum1, sum2;

  sum1 = array_sum(v1->x, v1->num_dims);
  sum2 = array_sum(v2->x, v2->num_dims);
  
  if (sum1 > sum2)
    return -1;
  else if (sum2 > sum1)
    return 1;
  else
    return 0;
}

/* MCB_compar_MAXRATIO */
int MCB_compar_MAXRATIO(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;
  float min1, min2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);
  min1 = array_min(v1->x, v1->num_dims);
  min2 = array_min(v2->x, v2->num_dims);
  
  if (max1/min1 > max2/min2)
    return -1;
  else if (max2/min2 > max1/min1)
    return 1;
  else
    return 0;
}

/* MCB_compar_MAXDIFF */
int MCB_compar_MAXDIFF(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;
  float min1, min2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);
  min1 = array_min(v1->x, v1->num_dims);
  min2 = array_min(v2->x, v2->num_dims);
  
  if (max1-min1 > max2-min2)
    return -1;
  else if (max2-min2 > max1-min1)
    return 1;
  else
    return 0;
}

int int_float_compar(const void *x, const void *y)
{
  struct int_float ifx = *(struct int_float *)x;
  struct int_float ify = *(struct int_float *)y;

  if (ifx.f > ify.f)
    return 1;
  else if (ifx.f < ify.f)
    return -1;
  else
    return 0;
}

/* helper function 
 *   Takes a bin_load and returns an array of dimension indices
 *   sorted from the lowest loaded dim to the highest loaded dim 
 */
void sort_dims_by_increasing_load(int num_dims, int *sorted_dims, float *bin_load)
{
  struct int_float *sorted;
  int i;

  sorted = (struct int_float*)calloc(num_dims,sizeof(struct int_float));
  for (i=0; i<num_dims; i++) {
    sorted[i].i = i;
    sorted[i].f = bin_load[i];
  }
  
  qsort(sorted,num_dims,sizeof(struct int_float),int_float_compar);
  for (i=0; i<num_dims; i++)
    sorted_dims[i] = sorted[i].i;
  free(sorted);
}

void add_to_list_of_lists(int first_dim, int second_dim, 
                         int *num_lists, int **first_dims, int **second_dims)
{
  int i;

  // If already in the lists, foret it
  for (i=0; i < *num_lists; i++) {
    if ((first_dim == (*first_dims)[i]) && 
        (second_dim == (*second_dims)[i])) {
        return;
    }
  }

  // Add to the list 
  (*num_lists)++;
  *first_dims = (int *)REALLOC(*first_dims, (*num_lists) * sizeof(int));
  *second_dims = (int *)REALLOC(*second_dims, (*num_lists) * sizeof(int));
  (*first_dims)[(*num_lists)-1] = first_dim;
  (*second_dims)[(*num_lists)-1] = second_dim;

  return;
}

/* MCB_find_and_extract_feasible_vector_in_list() 
 *  Given a vector list, and a bin load, find the first vector
 *  that fits, and remove it from the list (return it)
 *
 *  The extraction just erases the end of the list (with shifting)
 *  just to be a bit more efficient
 */
struct vp_vector *MCB_find_and_extract_feasible_vector_in_list(
                  int num_dims, float *bin_load, 
                  struct vp_vector **list, int *list_size)
{
  struct vp_vector *vector;
  int i,j;

  for (i=0; i < *list_size; i++) {
    // fprintf(stderr,"  Looking at: ");
    // for (j=0; j < num_dims; j++) 
    //   fprintf(stderr,"%.2f ",list[i]->x[j]);
    // fprintf(stderr,"\n");
    
    for (j=0; j < num_dims; j++) {
      if (1.0 - bin_load[j] < list[i]->x[j])
        break;
    }
    if (j >= num_dims) 
      break;
  }
  
  if (i >= *list_size) 
    return NULL;

  // extract it
  vector = list[i];
  for (j=i; j < *list_size-1; j++)
    list[j] = list[j+1];
  (*list_size)--;
    
  return vector; 
}

//#define MCB_DEBUG 0


/* MCB_pick_vector(): The core of the BCB implementation 
 *   right now implements CP only
 */
struct vp_vector *MCB_pick_vector(int isPP,
                                  int num_dims, float *bin_load,
                                  struct vp_vector ****lists, 
                                  int **list_sizes)
{
  int i,j;
  int *sorted_dims;
  int bottom_d, top_d;

  //fprintf(stderr,"IN MCB_PICK_VECTOR\n");

  // sort the bin dims by increasing load
  sorted_dims = (int *)calloc(num_dims,sizeof(int));
  sort_dims_by_increasing_load(num_dims, sorted_dims, bin_load);

  // Go through the grinds...
  //   **.....  (best case)
  //   *.*....  (d1 ok, and other)
  //   .**....  (d2 ok, and other)
  //   *..*...  etc..
  //   .*.*...
  //   *...*..
  //   .*..*..

  // Construct an array of list indices to look at, in the
  // right order 
  int num_lists = 0;
  int *first_dims=NULL, *second_dims=NULL;

  for (bottom_d = 1; bottom_d < num_dims; bottom_d++) {
    int first_dim, second_dim;
    for (top_d = 0; top_d < bottom_d; top_d++) {
      // Top dim and d
      first_dim  = sorted_dims[bottom_d];
      second_dim = sorted_dims[top_d];
      if (isPP) {
        if (first_dim < second_dim) {
         add_to_list_of_lists(first_dim, second_dim, 
                         &num_lists, &first_dims, &second_dims);
         add_to_list_of_lists(second_dim, first_dim, 
                         &num_lists, &first_dims, &second_dims);
        } else {
         add_to_list_of_lists(second_dim, first_dim, 
                         &num_lists, &first_dims, &second_dims);
         add_to_list_of_lists(first_dim, second_dim, 
                         &num_lists, &first_dims, &second_dims);
        }
      } else {
        first_dim  = MAX(sorted_dims[bottom_d],sorted_dims[top_d]);
        second_dim = MIN(sorted_dims[bottom_d],sorted_dims[top_d]);
        add_to_list_of_lists(first_dim, second_dim, 
                         &num_lists, &first_dims, &second_dims);
      }
    }
  }

  free(sorted_dims);

#ifdef MCB_DEBUG
  fprintf(stderr,"Lists order: ");
  for (i=0; i< num_lists; i++)
    fprintf(stderr,"[%d][%d],",first_dims[i],second_dims[i]);
  fprintf(stderr,"\n");
#endif

  // Go through the lists in order and look for a feasible vector to pick
  struct vp_vector *vector = NULL;
  for (i=0; i < num_lists; i++) {
    struct vp_vector **list = lists[first_dims[i]][second_dims[i]];
    int *list_size = &(list_sizes[first_dims[i]][second_dims[i]]);
#ifdef MCB_DEBUG
    fprintf(stderr,"  Looking for a feasible vector in LIST[%d][%d] (size=%d):",
                    first_dims[i],second_dims[i],*list_size);
    for (j=0; j < *list_size; j++)
      fprintf(stderr,"%d ",list[j]->service);
    fprintf(stderr,"\n");
#endif
    vector = MCB_find_and_extract_feasible_vector_in_list(
                  num_dims, bin_load, list, list_size); 
    if (vector) {// We found a vector in the list
//      fprintf(stderr,"    Picked vector corr to svc. %d in LIST[%d][%d]\n",
//                vector->service, first_dims[i],second_dims[i]);
      // fprintf(stderr,"  Found vector: ");
      // int j;
      // for (j=0; j < vector->num_dims; j++)
      //   fprintf(stderr,"%.2f ",vector->x[j]);
      // fprintf(stderr,"\n");
      break;
    }
  }

  free(first_dims); free(second_dims);
 
  // If we found a vector, then great
  if (vector) 
    return vector;
  else
    return NULL;
}

/*
 * The MCB implementation
 *   "PP" (permutation pack) or "CP" (choose pack)
 *   sorting type: "MAX", "SUM", "MAXDIFF", "MAXRATIO"
 * 
 *  Everything implemented for w=2
 */
int solve_vp_instance_MCB(char *pack, char *sort_type, 
                          struct vp_instance *vp)
{
  int num_lists,i,j,d1,d2;
  int (*compar)(const void *, const void *);
  char isPP,isCP;

  if (!strcmp(pack,"PP")) {
    isPP = 1; isCP = 0;
  } else if (!strcmp(pack,"CP")) {
    isPP = 0; isCP = 1;
  } else {
    fprintf(stderr,"MCB: unknown pack type '%s'\n",pack);
    exit(1);
  }

#ifdef MCB_DEBUG
  fprintf(stderr,"#############################################\n");
  fprintf(stderr,"########## SOLVING AN INSTANCE ##############\n");
  fprintf(stderr,"#############################################\n");
#endif

  // Initialize lists
  struct vp_vector ****lists;
  int **list_sizes;

  lists = (struct vp_vector ****)calloc(vp->num_dims,
             sizeof(struct vp_vector ***));
  list_sizes = (int **)calloc(vp->num_dims, 
             sizeof(int *));
  for (d1=0; d1 < vp->num_dims; d1++) {
    lists[d1] = (struct vp_vector ***)calloc(vp->num_dims,
             sizeof(struct vp_vector **));
    list_sizes[d1] = (int *)calloc(vp->num_dims, 
             sizeof(int));
  }

  // For CP we need only to initialize the triangular
  // matrix but what the heck.. this doesn't hurt
  for (d1=0; d1 < vp->num_dims; d1++) {
    for (d2=0; d2 < vp->num_dims; d2++) {
      lists[d1][d2] = NULL;
      list_sizes[d1][d2] = 0;
    }
  }

  // Put vectors in lists
  for (i=0; i < vp->num_vectors; i++) {
    int d1, d2;

    // Find the indices of the two largest coordinates of the vector
    // with d1 being the index of the dimension with the 
    // largest value and d2 the index of the dimension with
    // second largest value
    find_top_two_dims(vp->vectors[i].x, 
                      vp->vectors[i].num_dims,&d1,&d2);

    // If CP, then sort the dims indices, since we only
    // need a triangular matrix
    if (isCP) {
      if (d1 < d2) {
        int tmp = d1;
        d1 = d2;
        d2 = tmp;
      }
    }

    // Put the vector in the list
#ifdef MCB_DEBUG
    fprintf(stderr,"Putting vector ");
    for (j=0; j < vp->num_dims; j++)
      fprintf(stderr,"%.2f ",vp->vectors[i].x[j]);
    fprintf(stderr," (svc %d) in LIST[%d][%d]\n",
                  vp->vectors[i].service,
                  d1,d2);
#endif
    list_sizes[d1][d2]++;
    lists[d1][d2] = (struct vp_vector **)REALLOC(lists[d1][d2],
                         list_sizes[d1][d2] * sizeof (struct vp_vector *));     
    lists[d1][d2][list_sizes[d1][d2] - 1] = &(vp->vectors[i]);
  }

  // Sort the lists according to the criteria

  // Select the sorting function
  if (!strcmp(sort_type,"MAX")) {
    compar = MCB_compar_MAX;
  } else if (!strcmp(sort_type,"SUM")) {
    compar = MCB_compar_SUM;
  } else if (!strcmp(sort_type,"MAXDIFF")) {
    compar = MCB_compar_MAXDIFF;
  } else if (!strcmp(sort_type,"MAXRATIO")) {
    compar = MCB_compar_MAXRATIO;
  } else {
    fprintf(stderr,"Invalid MCB compare type '%s'\n",sort_type);
    exit(1);
  }

  // Sort all the lists
  for (d1=0; d1 < vp->num_dims; d1++) {
    int d2_loop_bound;
    if (isCP)
      d2_loop_bound = d1;  // only triangular
    else
      d2_loop_bound = vp->num_dims; // full matrix
    for (d2=0; d2 < d2_loop_bound; d2++) {
      qsort(lists[d1][d2],list_sizes[d1][d2], sizeof(struct vp_vector *), compar);
#ifdef MCB_DEBUG
      fprintf(stderr,"  Sorted LIST[%d][%d] (size=%d):",
                      d1,d2,list_sizes[d1][d2]);
      for (j=0; j < list_sizes[d1][d2]; j++)
        fprintf(stderr,"%d ",lists[d1][d2][j]->service);
      fprintf(stderr,"\n");
#endif
    }
  }

  // Go through bins
  int bin=0;
  int num_vectors_picked=0;

  // This is conservative: at most one bin per service
  // Could be replaced with a dynamic array, but I opted
  // for the simple, if a bit wasteful, solution
  float load[vp->num_vectors][vp->num_dims];

  // Set the loads to zero
  for (i=0; i < vp->num_vectors; i++) 
    for (j=0; j < vp->num_dims; j++)
      load[i][j] = 0.0;
  
  while (1) {
    struct vp_vector *vector;
    // Pick the vector to put in the bin
    // This call removes the vector from lists
#ifdef MCB_DEBUG
    fprintf(stderr,"** BIN %d load= ",bin);
    for (j=0; j < vp->num_dims; j++)
      fprintf(stderr,"%.2f ",load[bin][j]);
    fprintf(stderr,"\n");
#endif
 
    vector = MCB_pick_vector(isPP, vp->num_dims, load[bin],  lists, list_sizes);
    if (!vector) { // Couldn't find a vector, go to next bin
#ifdef MCB_DEBUG
      fprintf(stderr,"Can't find a vector, adding a new bin\n");
#endif
      bin++;
      continue;
    } 
#ifdef MCB_DEBUG
    fprintf(stderr,"Picked vector: ");
    for (j=0; j < vp->num_dims; j++)
      fprintf(stderr,"%.2f ",vector->x[j]);
    fprintf(stderr," (svc %d)\n",vector->service);
#endif
    // Update the mapping
    vp->mapping[vector->service] = bin;
    // Update the load of the bin
    for (j=0; j < vp->num_dims; j++)
      load[bin][j] += vector->x[j];
    // Are we done?
    num_vectors_picked++;
    if (num_vectors_picked == vp->num_vectors)
      break; 
  }

  // Free some memory
  for (d1=0; d1 < vp->num_dims; d1++) {
    free(list_sizes[d1]);
    for (d2=0; d2 < vp->num_dims; d2++)  {
      if (lists[d1][d2]) 
        free(lists[d1][d2]);
    }
    free(lists[d1]);
  }
  free(lists);
  free(list_sizes);

  // Return the number of bins used
  vp->num_bins = bin+1;
  return VP_SOLVED;
}

//#define CHEKURI_DEBUG 1


glp_prob *CHEKURI_create_lp(struct vp_instance *vp, int m)
{
  glp_prob *prob;
  int i,j,k;
  int n = vp->num_vectors;
  int d = vp->num_dims;
  int counter;
  int obj;
  int row;

  // We use the indices from the Cheruki paper:
  //   - i: vectors
  //   - j: bins
  //   - k: resource dimensions

  // Create GLPK problem
  glp_term_hook(&glp_stderr_out, NULL);
  prob = glp_create_prob();
  glp_set_prob_name(prob, "CHEKURI task placement");

  // Number of variables (called "columns"): 
  //  x_ij: 1 .. n*m
  //  plus a bogus objective
//  fprintf(stderr,"Adding %d columns\n",n*m+1);
  glp_add_cols(prob, n*m+1);

  // Defining convenient variables to reference 
  // columns without too many index calculations
  int **xij;
  xij = (int **)calloc(n+1,sizeof(int*));
  for (i=1; i<=n; i++)
    xij[i] = (int *)calloc(m+1,sizeof(int));
  counter = 1;
  for (i=1; i<=n; i++) {
    for (j=1; j<=m; j++) {
      xij[i][j] = counter++;
    }
  }

  // Set bounds for the e_ij: binary if !rational
  counter = 1;
  for (i=1; i<= n; i++) {
    for (j=1; j<= m; j++) {
      char buffer[10];
      sprintf(buffer,"x_{%d,%d}",i,j);
      glp_set_col_name(prob, counter, buffer);
      glp_set_col_bnds(prob, counter, GLP_DB, 0.0, 1.0);
      counter++;
    }
  }

  //fprintf(stderr,"obj = %d\n",counter);
  obj = counter;
  glp_set_col_name(prob, obj, "bogus objective");
  glp_set_col_bnds(prob, obj, GLP_DB, 0.0, 1.0);

  // objective is minyield
  glp_set_obj_name(prob, "bogus objective (we only need feasible)");
  glp_set_obj_dir(prob, GLP_MAX);
  glp_set_obj_coef(prob, obj, 1.0);

  // Number of rows  
  //  (A) for all i 	   sum_j e_ij = 1:       		  n rows
  //  (B) for all k      sum_i p_ik x_ij <= 1		          d*m rows
  //fprintf(stderr,"adding %d rows\n",n+n*m);
  glp_add_rows(prob, n + n*m);

  // compute the number of non-zero matrix elements
  int ne = n*m + m*d*n;
  //fprintf(stderr,"non-empty: %d\n",ne);
  int *ia = (int *)calloc(ne+1,sizeof(int));
  int *ja = (int *)calloc(ne+1,sizeof(int));
  double *ra = (double *)calloc(ne+1,sizeof(double));
  
  counter = 1;
  // Constraint (A)
  row = 1;
  for (i=1; i <= n; i++) {
    glp_set_row_bnds(prob,row,GLP_FX,1.0,1.0);
    for (j=1; j <= m; j++) {
      ia[counter] = row; 
      ja[counter] = xij[i][j];
      ra[counter] = 1.0; 
      //fprintf(stderr,"%d: %d %d %f\n",counter,ia[counter],ja[counter],ra[counter]);
      counter++;
    }
    row++;
  }
   
  // Constraint (B)
  row = n+1;
  for (k=0; k < d; k++) {
    for (j=1; j <= m; j++) {
      glp_set_row_bnds(prob,row,GLP_UP,-666,1.0);
      for (i=1; i <= n; i++) {
        ia[counter] = row; 
        ja[counter] = xij[i][j]; 
        ra[counter] = vp->vectors[i-1].x[k];  
      //fprintf(stderr,"%d: %d %d %f\n",counter,ia[counter],ja[counter],ra[counter]);
        counter++;
      }
      row++;
    }
  }

  glp_load_matrix(prob, ne, ia, ja, ra);

  free_2D_int_array(xij, n+1);
  free(ia);
  free(ja);
  free(ra);

  return prob;
}

/* CHEKURI_process_integer_assignments() */
int CHEKURI_process_integer_assignments(
         struct vp_instance *vp, glp_prob *prob)
{
  int i,j,k;
  int counter;
  int num_ints;

  // Defining convenient variables to reference 
  // columns without too many index calculations
  int **xij;
  xij = (int **)calloc(vp->num_vectors+1,sizeof(int*));
  for (i=1; i <= vp->num_vectors; i++)
    xij[i] = (int *)calloc(vp->num_bins+1,sizeof(int));
  counter = 1;
  for (i=1; i <= vp->num_vectors; i++) {
    for (j=1; j <= vp->num_bins; j++) {
      xij[i][j] = counter++;
    }
  }

  // Go through the vectors
  num_ints = 0;
  for (i=0; i < vp->num_vectors; i++) {
    for (j=0; j < vp->num_bins; j++) {
      float value = glp_get_col_prim(prob, xij[i+1][j+1]);
//      fprintf(stderr,"[%d,%d]=%.3f\n",i+1,j+1,glp_get_col_prim(prob, xij[i+1][j+1]));
      if (value > 1 - EPSILON) {
        vp->mapping[i] = j;
        num_ints++;
      } else if (value < EPSILON) {
        num_ints++;
      }
    }
  }

  free_2D_int_array(xij,vp->num_vectors + 1);

#if 0
  fprintf(stderr,"CHEKURI: num integers = %d\n",num_ints);
  fprintf(stderr,"CHEKURI: num non-integers = %d\n",vp->num_vectors * vp->num_bins - num_ints);
  fprintf(stderr,"CHEKURI: Max num of non-integers = %d\n",vp->num_bins * vp->num_dims);
#endif
  // Sanity check 
  if (vp->num_vectors * vp->num_bins - num_ints > vp->num_dims * vp->num_bins) {
//    fprintf(stderr,"CHEKURI: There are too many non-integer variables values!\n");
//    return 1;
  }

  return 0;
}

/* CHEKURI_get_unmapped_vectors() */
int CHEKURI_get_unmapped_vectors(struct vp_instance *vp, 
                                 int **remaining_unmapped_vectors) 
{
  int i;
  int count=0;
  int *list = NULL;

  for (i=0; i < vp->num_vectors; i++) {
    if (vp->mapping[i] != -1)
      continue;
    count++;
    list = (int *)REALLOC(list, count * sizeof(int));
    list[count-1] = i;
  }
  *remaining_unmapped_vectors = list;
  return count;
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

/* CHEKURI_process_non_integer_assignments() */
int CHEKURI_process_non_integer_assignments(struct vp_instance *vp)
{
  int *remaining_unmapped_vectors;
  int num_remaining_unmapped_vectors;
  int k;
  int last_max_subset_size = -1;


  while (num_remaining_unmapped_vectors = 
             CHEKURI_get_unmapped_vectors(vp, &remaining_unmapped_vectors)) {
    int subset_size;
    int *subset;

#ifdef CHEKURI_DEBUG
     fprintf(stderr,"** remaining unmapped vectors: ");
     for (k = 0; k < num_remaining_unmapped_vectors; k++)
       fprintf(stderr,"%d ",remaining_unmapped_vectors[k]);
     fprintf(stderr,"\n");
#endif

    // find the maximum subset
    subset_size = find_maximum_subset(vp,remaining_unmapped_vectors,
                        num_remaining_unmapped_vectors, &subset, 
                        last_max_subset_size);
    last_max_subset_size = subset_size;
    free(remaining_unmapped_vectors);
#ifdef CHEKURI_DEBUG
    fprintf(stderr,"Found a maximum subset of size %d: ",subset_size);
    for (k=0; k < subset_size; k++)
      fprintf(stderr,"%d ",subset[k]);
    fprintf(stderr,"\n");
#endif
  

    if (subset_size == 0) {
      fprintf(stderr,"CANNOT FIND A SUBSET!!\n");
      return 1;
    }

    // update the mappings
    //fprintf(stderr,"Updating the mappings\n");
    for (k=0; k < subset_size; k++) {
      vp->mapping[subset[k]] = vp->num_bins;
      //fprintf(stderr,"  setting vp->mapping[%d] to %d\n", subset[k], vp->num_bins);
    }
    (vp->num_bins)++;
    free(subset);
  }
  
  return 0;
}

//#define CHEKURI_DEBUG 1

/* Cheruki's guaranteed algorithm */
int solve_vp_instance_CHEKURI(struct vp_instance *vp)
{
  glp_prob *prob = NULL;
  glp_prob *last_solved_prob = NULL;
  int m, m_lo, m_hi;
  int last_successful_num_bins;

  m_hi = vp->num_vectors;
  m_lo = 1;
  m = -1;

  #ifdef CHEKURI_DEBUG
  fprintf(stderr,"Starting the binary search\n");
  #endif
  while(m_hi != m_lo) {
   
    if ((prob) && (prob != last_solved_prob))
      glp_delete_prob(prob);

    if (m == ceilf(((float)m_hi + (float)m_lo)/2.0)) {
      m = m_lo;
      m_hi = m_lo;
    } else {
      m = ceilf(((float)m_hi + (float)m_lo)/2.0);
    }
    #ifdef CHEKURI_DEBUG
    fprintf(stderr,"Trying with %d bins (lo=%d  hi=%d)\n",m,m_lo,m_hi);
    #endif
    
    // Create the LP problem
    prob = CHEKURI_create_lp(vp,m);

    //fprintf(stderr,"Writing CHEKURI LP to file 'linearprogram'..");
    glp_write_lp(prob,NULL,"linearprogram");
    //fprintf(stderr,"done\n");

    // Solve it
    if (solve_linear_program(RATIONAL, NOT_VERBOSE, prob)) {
      m_lo = m;
    } else {
      m_hi = m; 
      last_solved_prob = prob;
      last_successful_num_bins = m;
    }
    
  }

#ifdef CHEKURI_DEBUG
  fprintf(stderr,"Solvable with %d bins minimum\n",
                 last_successful_num_bins);
#endif

  // Compute the "output" of the vp instance
  vp->num_bins = last_successful_num_bins;

  // integer assignments
  if (CHEKURI_process_integer_assignments(vp,last_solved_prob)) {
    fprintf(stderr,"CHEKURI: cannot process integer assignments\n");
    return VP_NOT_SOLVED;
  }

  // Done with the linear program
  glp_delete_prob(last_solved_prob);
  

  // non-integer assignments
  if (CHEKURI_process_non_integer_assignments(vp)) {
    fprintf(stderr,"CHEKURI: cannot process integer assignments\n");
    return VP_NOT_SOLVED;
  }

  return VP_SOLVED;
}

/* solve_vp_instance() 
 *  Returns the number of bins used
 */
int solve_vp_instance(struct vp_instance *vp, 
                      char *vp_algorithm)
{
  int status;

  if (!strcmp(vp_algorithm,"FFDLEX")) {
    status = solve_vp_instance_FITD("FIRST","LEX",vp);
  } else if (!strcmp(vp_algorithm,"FFDMAX")) {
    status = solve_vp_instance_FITD("FIRST","MAX",vp);
  } else if (!strcmp(vp_algorithm,"FFDSUM")) {
    status = solve_vp_instance_FITD("FIRST","SUM",vp);
  } else if (!strcmp(vp_algorithm,"BFDLEX")) {
    status = solve_vp_instance_FITD("BEST","LEX",vp);
  } else if (!strcmp(vp_algorithm,"BFDMAX")) {
    status = solve_vp_instance_FITD("BEST","MAX",vp);
  } else if (!strcmp(vp_algorithm,"BFDSUM")) {
    status = solve_vp_instance_FITD("BEST","SUM",vp);
  } else if (!strcmp(vp_algorithm,"CPMAX")) {
    status = solve_vp_instance_MCB("CP","MAX",vp);
  } else if (!strcmp(vp_algorithm,"CPSUM")) {
    status = solve_vp_instance_MCB("CP","SUM",vp);
  } else if (!strcmp(vp_algorithm,"CPMAXRATIO")) {
    status = solve_vp_instance_MCB("CP","MAXRATIO",vp);
  } else if (!strcmp(vp_algorithm,"CPMAXDIFF")) {
    status = solve_vp_instance_MCB("CP","MAXDIFF",vp);
  } else if (!strcmp(vp_algorithm,"PPMAX")) {
    status = solve_vp_instance_MCB("PP","MAX",vp);
  } else if (!strcmp(vp_algorithm,"PPSUM")) {
    status = solve_vp_instance_MCB("PP","SUM",vp);
  } else if (!strcmp(vp_algorithm,"PPMAXRATIO")) {
    status = solve_vp_instance_MCB("PP","MAXRATIO",vp);
  } else if (!strcmp(vp_algorithm,"PPMAXDIFF")) {
    status = solve_vp_instance_MCB("PP","MAXDIFF",vp);
  } else if (!strcmp(vp_algorithm,"CHEKURI")) {
    status = solve_vp_instance_CHEKURI(vp);
  } else {
    fprintf(stderr,"Unknown vp_algorithm '%s'\n",vp_algorithm);
    exit(1);
  }
  return status;
}

/* 
 * VP_scheduler() 
 */
int VP_scheduler(char *vp_algorithm, char *ignore2, char *ignore3)
{
  double yield, yieldlb, yieldub, best_yield;
  struct vp_instance *vp = NULL;
  int i, status;
  int num_bins;

  yieldlb = 0.0;
  yieldub = compute_LP_bound();

  best_yield = -1.0;
  yield = 0.0;


  while((yield < 1.0 - EPSILON) && (yieldub - yieldlb > 0.001)) {
    yield = (yieldub + yieldlb) / 2.0;
#if 0
    fprintf(stderr,"Trying yield=%.6f (up=%.4f down=%.4f)\n",yield,yieldub,yieldlb);
#endif

    // Generate the VP instance
    if (vp)
      free_vp_instance(vp);
    vp = generate_vp_instance(yield);

    // Solve the VP instance
    if ((status = solve_vp_instance(vp, vp_algorithm)) != VP_SOLVED) {
      fprintf(stderr,"FAILED\n");
      return RESOURCE_ALLOCATION_FAILURE;
    }
     
    if (vp->num_bins <= INS.numservers) {
      yieldlb = yield; 
      best_yield = yield;
#if 0
      fprintf(stderr,"Instance was solved with %d bins\n",vp->num_bins);
#endif
      // Save the computed mapping
      for (i=0; i < INS.numservices; i++) {
        INS.mapping[vp->vectors[i].service] = vp->mapping[i];
        INS.allocation[vp->vectors[i].service] = best_yield;
#if 0
        fprintf(stderr,"SERVICE %d : on server %d with yield %.2f\n",
                vp->vectors[i].service, 
                mapping[vp->vectors[i].service], 
                allocation[vp->vectors[i].service]);
#endif
      }
    } else {
      yieldub = yield;
#if 0
      fprintf(stderr,"Instance was not solved (needed %d bins, and we have %d servers)\n",
              vp->num_bins,INS.numservers);
#endif
    }
  }
  if (vp) 
    free_vp_instance(vp); // free the last instance,vp->num_bins

  // We never solved the VP instance
  if (best_yield <= 0)
    return RESOURCE_ALLOCATION_FAILURE;

  // We set the yield and mapping of all services

  return RESOURCE_ALLOCATION_SUCCESS;
}

int GA_scheduler(char *arg1, char *arg2, char *arg3)
{
  int i;
  int initial = !strcmp(arg1,"F");
  int mutated = !strcmp(arg2,"F");
  int crossed = !strcmp(arg3,"F");

  // Map each service to a particular host
  if (GA_compute_mapping(initial,mutated,crossed) == RESOURCE_ALLOCATION_FAILURE) {
    return RESOURCE_ALLOCATION_FAILURE; 
  }
#if 0
  fprintf(stderr,"GA MAPPING:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d on server %d\n",i,INS.mapping[i]);
  }
#endif

  // For each host, compute its services' allocations
  for (i=0; i<INS.numservers; i++) {
    compute_allocations_given_mapping(i);
  } 
#if 0
  fprintf(stderr,"GA ALLOCATIONS:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d has an allocation (yield) of %f\n",i,INS.allocation[i]);
  }

#endif
  
  return RESOURCE_ALLOCATION_SUCCESS;
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
  int i,j;
  
  for (j=0; j<INS.numfluid; j++) {
    float free_resource = 1.0;
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
      if (load > 1.0 + EPSILON) {
        fprintf(stderr,"Error: Rigid Capacity %d of server %d exceeded (%f)\n", j, server, load);
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
        load += INS.fluidneeds[service][j] * 
           (INS.allocation[service] + INS.slas[service]*(1 - INS.allocation[service]));
//        fprintf(stderr,"load due to service %d: %f * (%f + %f * (1 - %f)) = %f\n",service,fluidneeds[service][j], allocation[service], slas[service],allocation[service], fluidneeds[service][j] * (allocation[service] + slas[service]*(1 - allocation[service])));
      }
      if (load > 1.0 + EPSILON) {
        fprintf(stderr,"Error: Fluid Capacity %d of server %d exceeded (%f)\n",
                   j, server, load);
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

    INS.slas = (float *)calloc(INS.numservices,sizeof(float));
    INS.rigidneeds = (float **)calloc(INS.numservices,sizeof(float*));
    for (i=0; i<INS.numservices; i++)
      INS.rigidneeds[i] = (float *)calloc(INS.numrigid,sizeof(float));
    INS.fluidneeds = (float **)calloc(INS.numservices,sizeof(float*));
    for (i=0; i<INS.numservices; i++)
      INS.fluidneeds[i] = (float *)calloc(INS.numfluid,sizeof(float));
    INS.mapping = (int *)calloc(INS.numservices,sizeof(int));
    INS.allocation = (float *)calloc(INS.numservices,sizeof(float));

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

      // if we did find it, then run it
      gettimeofday(&time1, NULL);

      // Initialize the server loads 
      initialize_global_server_loads();

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
//    fprintf(output,"##################\n");

    for (schedulerptr = schedulers; *schedulerptr; schedulerptr++) 
      free(*schedulerptr);
    free(schedulers);
    
    return 0;
}
