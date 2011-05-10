#include "flexsched.h"

int glp_stderr_null_out(void *info, const char *s) 
{
    return 1;
}

/* GLPK  solver
 * Returns 0 on success
 * Returns 1 if couldn't solve
 * Returns 2 if time limit was reached
 */
int solve_linear_program(glp_prob *prob, int rational)
{
    int solver_status, solution_status;

#ifdef DEBUG
    fprintf(stderr, "Writing LP to file 'linearprogram'...");
    glp_write_lp(prob, NULL, "linearprogram");
    fprintf(stderr, "done\n");
#else
    // kill lp output
    glp_term_hook(&glp_stderr_null_out, NULL);
#endif

    if (rational) {
#ifdef DEBUG
        fprintf(stderr, "Solving a rational LP...\n");}
#endif
        solver_status = glp_simplex(prob, NULL);
        solution_status = (glp_get_status(prob) != GLP_OPT);
    } else {
#ifdef DEBUG
        fprintf(stderr, "Solving a MILP...\n");
#endif
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        parm.tm_lim = GLPK_TIME_LIMIT;
        solver_status = glp_intopt(prob, &parm);
        solution_status = (glp_mip_status(prob) != GLP_OPT);
    }
#ifdef DEBUG
    fprintf(stderr, "solver_stat: %d, ", solver_status);
    fprintf(stderr, "solution_stat: %d\n", solution_status);
    fprintf(stderr, "Writing LP solution to file 'solution'...");
    glp_write_sol(prob, "solution");
    fprintf(stderr, "done\n");
#endif

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

// Number of variables (called "columns"): 
//  e_ij = H*(i-1)+j: 1 .. NH
//  y_ih = N*H+(i-1)+j: NH+1 .. 2NH
//  Y: 2NH+1
#define LP_E_IJ_COL (H*i+j)
#define LP_Y_IJ_COL (N*H+H*i+j)
#define LP_OBJ_COL  (2*N*H+1)
#define LP_NUM_COLS (2*N*H+1)

// Number of rows  
//  (A) for all i      sum_h e_ih = 1:                N rows
//  (B) for all i,h    y_ih <= e_ih                   NH rows
//  (C) for all i      sum_h y_ih >= sla_i            N rows
//  (D) for all h, jr  sum_i r_ij*e_ih <= 1           HR rows
//  (E) for all h, jf  sum_i r_ij*y_ih <= 1           HF rows
//  (F) for all i      (sum_h y_ih - sla_i)/(1-sla_i) >= Y    N rows
#define LP_NUM_ROWS (N+N*H+N+H*R+H*F+N)

// non-zero matrix elements
// Constraint A: N rows of H elements -- N*H
// Constraint B: N*H rows of 2 elements -- N*H*2
// Constraint C: N rows of H elements -- N*H
// Constraint D: H*R rows of N elements -- H*R*N
// Constraint E: H*F rows of N elements -- H*F*N
// Constraint F: N rows of H+1 elements -- N*(H+1);
#define LP_NUM_ELTS (N*H + N*H*2 + N*H + H*R*N + H*F*N + N*(H+1))

glp_prob *create_placement_lp(int rational) 
{
    glp_prob *prob;
    int N = INS.numservices;
    int H = INS.numservers;
    int R = INS.numrigid;
    int F = INS.numfluid;
    int i, j, k, r, d;
    int ia[LP_NUM_ELTS+1];
    int ja[LP_NUM_ELTS+1];
    double ra[LP_NUM_ELTS+1];


    // Create GLPK problem
    prob = glp_create_prob();

    glp_add_cols(prob, LP_NUM_COLS);
    glp_add_rows(prob, LP_NUM_ROWS);

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_bnds(prob, LP_OBJ_COL, GLP_DB, 0.0, 1.0);
    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, LP_OBJ_COL, 1.0);

#ifdef DEBUG
    glp_set_prob_name(prob, "task placement problem");
    glp_set_col_name(prob, LP_OBJ_COL, "Y");
    glp_set_obj_name(prob, "minimum scaled yield");
#endif

    // Set bounds for the e_ij: binary if !rational
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= H; j++) {
#ifdef DEBUG
            char buffer[10];
            sprintf(buffer, "e_{%d,%d}", i, j);
            glp_set_col_name(prob, LP_E_IJ_COL, buffer);
            sprintf(buffer,"y_{%d,%d}", i, j);
            glp_set_col_name(prob, LP_Y_IJ_COL, buffer);
#endif
            if (rational) {
                glp_set_col_kind(prob, LP_E_IJ_COL, GLP_CV);
                glp_set_col_bnds(prob, LP_E_IJ_COL, GLP_DB, 0.0, 1.0);
            } else {
                glp_set_col_kind(prob, LP_E_IJ_COL, GLP_BV);
            }
            glp_set_col_bnds(prob, LP_Y_IJ_COL, GLP_DB, 0.0, 1.0);
        }
    }

    k = 1; // element counter
    r = 1; // row counter

    //  (A) for all i: sum_j e_ij = 1
    for (i = 1; i <= N; i++) {
        glp_set_row_bnds(prob, r, GLP_FX, 1.0, 1.0);
        for (j = 1; j <= H; j++) {
            ia[k] = r; ja[k] = LP_E_IJ_COL; ra[k] = 1.0; k++;
        }
        r++;
    }

    //  (B) for all i, j:  y_ih <= e_ij <=> y_ij - e_ij <= 0.0
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= H; j++) {
            glp_set_row_bnds(prob, r, GLP_UP, -666, 0.0);
            ia[k] = r; ja[k] = LP_Y_IJ_COL; ra[k] = 1.0;  k++;
            ia[k] = r; ja[k] = LP_E_IJ_COL; ra[k] = -1.0; k++;
            r++;
        }
    }
   
    //  (C) for all i: sum_h y_ih >= sla_i
    for (i = 1; i <= N; i++) {
        glp_set_row_bnds(prob, r, GLP_LO, INS.slas[i-1], 666);
        for (j = 1; j <= H; j++) {
            ia[k] = r; ja[k] = LP_Y_IJ_COL; ra[k] = 1.0; k++;
        }
        r++;
    }

    //  (D) for all h, r:  sum_i r_ij*e_ih <= 1
    for (j = 1; j <= H; j++) {
        for (d = 1; d <= R; d++) {
            glp_set_row_bnds(prob, r, GLP_UP, -666, 
                INS.rigidcapacities[j-1][d-1]);
            for (i = 1; i <= N; i++) {
                ia[k] = r; ja[k] = LP_E_IJ_COL; 
                ra[k] = INS.rigidneeds[i-1][d-1]; k++;
            }
            r++;
        }
    }

    //  (E) for all h, f:  sum_i r_ij*y_ih <= 1
    for (j = 1; j <= H; j++) {
        for (d = 1; d <= F; d++) {
            glp_set_row_bnds(prob, r, GLP_UP, 666, 
                INS.fluidcapacities[j-1][d-1]);
            for (i = 1; i <= N; i++) {
                ia[k] = r; ja[k] = LP_Y_IJ_COL; 
                ra[k] = INS.fluidneeds[i-1][d-1]; k++;
            }
            r++;
        }
    }

    //  (F) for all i: (sum_h y_ih - sla_i)/(1-sla_i) >= Y <=>
    //      sum_h y_ih + Y (sla_i - 1) >= sla_i 
    for (i = 1; i <= N; i++) {
        glp_set_row_bnds(prob, r, GLP_LO, INS.slas[i-1], 666);
        for (j = 1; j <= H; j++) {
            ia[k] = r; ja[k] = LP_Y_IJ_COL; ra[k] = 1.0; k++;
        }
        ia[k] = r; ja[k] = LP_OBJ_COL; ra[k] = INS.slas[i-1] - 1.0; k++;
        r++;
    }

    glp_load_matrix(prob, LP_NUM_ELTS, ia, ja, ra);

    return prob;

}

int LPBOUND_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
    // Compute the bound directly
    INS.lpbound = compute_LP_bound();

    glp_prob *prob = create_placement_lp(RATIONAL);

    if (solve_linear_program(prob, RATIONAL)) {
      return RESOURCE_ALLOCATION_FAILURE;
    }

    if (fabs(INS.lpbound - glp_get_col_prim(prob, LP_OBJ_COL)) > EPSILON) {
        fprintf(stderr,"Error: Immediate bound = %.3f, LP bound = %.3f\n", 
            INS.lpbound, glp_get_col_prim(prob, LP_OBJ_COL));
        exit(1);
#ifdef DEBUG
    } else {
        fprintf(stderr,"AGREEMENT\n");
#endif
    }

  return RESOURCE_ALLOCATION_SUCCESS;
}

int MILP_scheduler(char *ignore1, char *ignore2, char *ignore3)
{
    int i, j;
    int status;
    glp_prob *prob;

    // Create the problem (in MILP form)
    prob = create_placement_lp(INTEGER);

    // Solve it
    if (status = solve_linear_program(INTEGER, NOT_VERBOSE, prob)) {
        if (status == 2) {
            strcat(INS.misc_output, "T");
         }
        return RESOURCE_ALLOCATION_FAILURE;
    }

    // Retrieve the mapping
    for (i = 0; i < INS.numservices; i++) {
        for (j = 0; j < INS.numservers; j++) {
            if (glp_mip_col_val(prob, LP_E_IJ_COL) >= 1.0 - EPSILON) {
                INS.mapping[i] = j;
                INS.allocation[i] = (glp_mip_col_val(prob, LP_Y_IJ_COL) - 
                    INS.slas[i]) / (1 - INS.slas[i]);
                break;
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr,"MILP solution: minimum yield = %.3f\n", 
        glp_mip_col_val(prob, LP_OBJ_COL));
#endif

    return RESOURCE_ALLOCATION_SUCCESS;
}

/*
 *  Wrapper for RRND and RRNZ
 *      - RRND: LPROUNDING(0.0)
 *      - RRNZ: LPROUNDING(xxx)
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
