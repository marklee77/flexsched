#include "flexsched.h"

/***********************************************************/
/* Schedulers based on the linear program described in the */
/* JPDC paper. Chekuri uses a different linear program and */
/* is contained in its own file                            */
/***********************************************************/

#define LP_IGNORE 666

// Number of variables (called "columns"): 
//  e_ij = H*(i-1)+j: 1 .. NH
//  y_ih = N*H+(i-1)+j: NH+1 .. 2NH
//  Y: 2NH+1
#define LP_NUM_COLS (2*(flex_prob->num_services)*(flex_prob->num_servers)+1)

// Number of rows  
//  (A) for all i      sum_h e_ih = 1:                N rows
//  (B) for all i,h    y_ih <= e_ih                   NH rows
//  (C) for all i      sum_h y_ih >= sla_i            N rows
//  (D) for all h, jr  sum_i r_ij*e_ih <= 1           HR rows
//  (E) for all h, jf  sum_i r_ij*y_ih <= 1           HF rows
//  (F) for all i      (sum_h y_ih - sla_i)/(1-sla_i) >= Y    N rows
#define LP_NUM_ROWS (((flex_prob->num_services)+(flex_prob->num_rigid)+\
    (flex_prob->num_fluid))*(flex_prob->num_servers)+\
    3*(flex_prob->num_services))

// non-zero matrix elements
// Constraint A: N rows of H elements -- N*H
// Constraint B: N*H rows of 2 elements -- N*H*2
// Constraint C: N rows of H elements -- N*H
// Constraint D: H*R rows of N elements -- H*R*N
// Constraint E: H*F rows of N elements -- H*F*N
// Constraint F: N rows of H+1 elements -- N*(H+1);
#define LP_NUM_ELTS ((flex_prob->num_services)*(flex_prob->num_servers)*\
    (5+(flex_prob->num_rigid)+(flex_prob->num_fluid))+(flex_prob->num_services))

#ifdef CPLEX
// CPLEX does zero indexed rows/cols
#define LP_E_IJ_COL ((flex_prob->num_servers)*i+j)
#define LP_Y_IJ_COL ((flex_prob->num_servers)*((flex_prob->num_services)+i)+j)
#define LP_OBJ_COL  (2*(flex_prob->num_services)*(flex_prob->num_servers))

#include <ilcplex/cplex.h>

int create_placement_lp(int rational, CPXENVptr *retenv, CPXLPptr *retlp) 
{
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int status = 0;
    int i, j, k;
    double obj, lb, ub;
    int row, elt;
    int ia[LP_NUM_ELTS+1];    // array of matrix element row indicies
    int ja[LP_NUM_ELTS+1];    // array of matrix element column indicies
    double ra[LP_NUM_ELTS+1]; // array of matrix element coefficients
    char type;
    char buffer[10];

    *retenv = NULL;
    *retlp = NULL; 

    // create CPLEX problem
    env = CPXopenCPLEX (&status);

    if ( env == NULL ) {
        char  errmsg[1024];
        fprintf (stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        return 1;
    }
  
    /* Turn on output to the screen */
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
        fprintf (stderr,
                 "Failure to turn on screen indicator, error %d.\n", status);
        return 1;
    }
  
    /* Turn on data checking */
    status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
    if ( status ) {
        fprintf (stderr, 
            "Failure to turn on data checking, error %d.\n", status);
        return 1; 
    }
      
    /* Create the problem. */
#ifdef DEBUG
#define LP_PROB_NAME "task placement problem"
#else
#define LP_PROB_NAME NULL
#endif
    lp = CPXcreateprob (env, &status, LP_PROB_NAME);


    if ( lp == NULL ) {
        fprintf (stderr, "Failed to create LP.\n");
        return 1;
    }  

#ifdef DEBUG
#define LP_COL_NAME &buffer
#else
#define LP_COL_NAME NULL
#endif

    // Set bounds for the e_ij: binary if !rational
    obj = 0.0;
    lb = 0.0;
    ub = 1.0;
    if (rational) {
        type = CPX_CONTINUOUS;
    } else {
        type = CPX_BINARY;
    }
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
#ifdef DEBUG
            sprintf(buffer, "e_{%d,%d}", i, j);
#endif
            CPXnewcols(env, lp, 1, &obj, &lb, &ub, &type, LP_COL_NAME);
        }
    }

    if (!rational) {
        type = CPX_CONTINUOUS;
    }
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
#ifdef DEBUG
            sprintf(buffer, "e_{%d,%d}", i, j);
#endif
            CPXnewcols(env, lp, 1, &obj, &lb, &ub, &type, LP_COL_NAME);
        }
    }
#ifdef DEBUG
            sprintf(buffer, "Y");
#endif

    // objective value
    obj = 1.0;
    CPXnewcols(env, lp, 1, &obj, &lb, &ub, &type, LP_COL_NAME);

    elt = 0; // element counter
    row = 0; // row counter

    //  (A) for all i: sum_j e_ij = 1
    type = 'E'
    lb = 1.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        CPXnewrows(env, lp, 1, &lb, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    //  (B) for all i, j:  y_ih <= e_ij <=> y_ij - e_ij <= 0.0
    type = 'L';
    ub = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            CPXnewrows(env, lp, 1, &ub, &type, NULL, NULL);
            glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 0.0);
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0;  elt++;
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    type = 'G';
    //  (C) for all i: sum_h y_ih >= sla_i
    for (i = 0; i < flex_prob->num_services; i++) {
        lb = flex_prob->slas[i];
        CPXnewrows(env, lp, 1, &lb, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    type = 'L';
    //  (D) for all h, r:  sum_i r_ij*e_ih <= 1
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_rigid; k++) {
            ub = flex_prob->rigid_capacities[j][k];
            CPXnewrows(env, lp, 1, &ub, &type, NULL, NULL);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = LP_E_IJ_COL; 
                ra[elt] = flex_prob->rigid_needs[i][k]; elt++;
            }
            row++;
        }
    }

    //  (E) for all h, f:  sum_i r_ij*y_ih <= 1
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            ub = flex_prob->fluid_capacities[j][k];
            CPXnewrows(env, lp, 1, &ub, &type, NULL, NULL);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = LP_Y_IJ_COL; 
                ra[elt] = flex_prob->fluid_needs[i][k]; elt++;
            }
            row++;
        }
    }

    //  (F) for all i: (sum_h y_ih - sla_i)/(1-sla_i) >= Y <=>
    //      sum_h y_ih + Y (sla_i - 1) >= sla_i 
    for (i = 0; i < flex_prob->num_services; i++) {
        lb = flex_prob->slas[i];
        CPXnewrows(env, lp, 1, &lb, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0; elt++;
        }
        ia[elt] = row; ja[elt] = LP_OBJ_COL; 
        ra[elt] = flex_prob->slas[i] - 1.0;
        elt++; row++;
    }

    status = CPXchgcoeflist (env, lp, LP_NUM_ELTS, ia, ja, ra);
    if (status) return 1;

    CPXchgobjsen (env, lp, CPX_MAX);

    *retenv = env;
    *retlp = lp;

    return 0;

}

int solve_linear_program(CPXENVptr env, CPXLPptr lp, int rational)
{
    int status;

    if (rational) {
        status = CPXlpopt(env, lp);
    } else {
        status = CPXmipopt(env, lp);
    }

  return status;
}

#else
#define LP_E_IJ_COL ((flex_prob->num_servers)*i+j+1)
#define LP_Y_IJ_COL ((flex_prob->num_servers)*((flex_prob->num_services)+i)+j+1)
#define LP_OBJ_COL  (2*(flex_prob->num_services)*(flex_prob->num_servers)+1)

#include <glpk.h>
glp_prob *create_placement_lp(int rational) 
{
    glp_prob *prob;
    int i, j, k;
    int row, elt;
    int ia[LP_NUM_ELTS+1];    // array of matrix element row indicies
    int ja[LP_NUM_ELTS+1];    // array of matrix element column indicies
    double ra[LP_NUM_ELTS+1]; // array of matrix element coefficients

    // Create GLPK problem
    prob = glp_create_prob();

    glp_add_cols(prob, LP_NUM_COLS);
    glp_add_rows(prob, LP_NUM_ROWS);

#ifdef DEBUG
    glp_set_prob_name(prob, "task placement problem");
    glp_set_col_name(prob, LP_OBJ_COL, "Y");
    glp_set_obj_name(prob, "minimum scaled yield");
#endif

    // Set bounds for the e_ij: binary if !rational
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
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

    elt = 1; // element counter
    row = 1; // row counter

    //  (A) for all i: sum_j e_ij = 1
    for (i = 0; i < flex_prob->num_services; i++) {
        glp_set_row_bnds(prob, row, GLP_FX, 1.0, 1.0);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    //  (B) for all i, j:  y_ih <= e_ij <=> y_ij - e_ij <= 0.0
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 0.0);
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0;  elt++;
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    //  (C) for all i: sum_h y_ih >= sla_i
    for (i = 0; i < flex_prob->num_services; i++) {
        glp_set_row_bnds(prob, row, GLP_LO, flex_prob->slas[i], LP_IGNORE);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    //  (D) for all h, r:  sum_i r_ij*e_ih <= 1
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_rigid; k++) {
            glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 
                flex_prob->rigid_capacities[j][k]);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = LP_E_IJ_COL; 
                ra[elt] = flex_prob->rigid_needs[i][k]; elt++;
            }
            row++;
        }
    }

    //  (E) for all h, f:  sum_i r_ij*y_ih <= 1
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 
                flex_prob->fluid_capacities[j][k]);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = LP_Y_IJ_COL; 
                ra[elt] = flex_prob->fluid_needs[i][k]; elt++;
            }
            row++;
        }
    }

    //  (F) for all i: (sum_h y_ih - sla_i)/(1-sla_i) >= Y <=>
    //      sum_h y_ih + Y (sla_i - 1) >= sla_i 
    for (i = 0; i < flex_prob->num_services; i++) {
        glp_set_row_bnds(prob, row, GLP_LO, flex_prob->slas[i], LP_IGNORE);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0; elt++;
        }
        ia[elt] = row; ja[elt] = LP_OBJ_COL; 
        ra[elt] = flex_prob->slas[i] - 1.0;
        elt++; row++;
    }

    glp_load_matrix(prob, LP_NUM_ELTS, ia, ja, ra);

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_bnds(prob, LP_OBJ_COL, GLP_DB, 0.0, 1.0);
    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, LP_OBJ_COL, 1.0);

    return prob;

}

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
        fprintf(stderr, "Solving a rational LP...\n");
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
#endif

flexsched_solution LPBOUND_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution("LPBOUND");
    double objval;
#ifdef CPLEX
    CPXENVptr env;
    CPXLPptr lp;
    create_placement_lp(RATIONAL, &env, lp);
#else
    glp_prob *prob = create_placement_lp(RATIONAL);
#endif

#ifdef CPLEX
    if (solve_linear_program(env, lp, RATIONAL)) {
#else
    if (solve_linear_program(prob, RATIONAL)) {
#endif
        return flex_soln;
    }

#ifdef CPLEX
    CPXgetobjval(env, lp, &objval);
#else
    objval = glp_get_col_prim(prob, LP_OBJ_COL);
#endif

    if (fabs(flex_prob->lpbound - objval) > EPSILON)
    {
        fprintf(stderr,"Error: Immediate bound = %.3f, LP bound = %.3f\n", 
            flex_prob->lpbound, objval);
        exit(1);
#ifdef DEBUG
    } else {
        fprintf(stderr,"AGREEMENT\n");
#endif
    }

    flex_soln->success = 1;
    return flex_soln;
}

flexsched_solution MILP_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    int i, j;
    int status;
    double val;
    flexsched_solution flex_soln = new_flexsched_solution("MILP");

#ifdef CPLEX

    CPXENVptr env;
    CPXLPptr lp;
    create_placement_lp(INTEGER, &env, lp);

    if (status = solve_linear_program(env, lp, INTEGER)) {
        return flex_soln;
    }

#else
    // Create the problem (in MILP form)
    glp_prob *prob = create_placement_lp(INTEGER);

    // Solve it
    if (status = solve_linear_program(prob, INTEGER)) {
        if (status == 2) {
            strcat(flex_soln->misc_output, "T");
        }
        return flex_soln;
    }

#endif

    // Retrieve the mapping
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
#ifdef CPLEX
            CPXgetx(env, lp, &val, LP_E_IJ_COL, LP_E_IJ_COL);
#else
            val = glp_mip_col_val(prob, LP_E_IJ_COL);
#endif
            if (val >= 1.0 - EPSILON) {
                flex_soln->mapping[i] = j;
#ifdef CPLEX
            CPXgetx(env, lp &val, LP_Y_IJ_COL, LP_Y_IJ_COL);
#else
            val = (glp_mip_col_val(prob, LP_Y_IJ_COL) - flex_prob->slas[i]) / 
                (1 - flex_prob->slas[i]);
#endif
                flex_soln->scaled_yields[i] = val;
                break;
            }
        }
    }

    flex_soln->success = 1;
    maximize_minimum_then_average_yield(flex_soln);
    return flex_soln;
}

/*
 *      - RRND: LPROUNDING(0.0)
 *      - RRNZ: LPROUNDING(xxx)
 */
void LPROUNDING_scheduler(
    flexsched_solution flex_soln, float min_weight)
{
    glp_prob *prob;
    int i, j;
    float weights[flex_prob->num_servers];
    float x, total_weight, select_over_weight;
    double val;

    // Create the placement problem in rational mode
    // and solve it to compute rational mappings,
    // adding an epsilon to zero values
#ifdef CPLEX

    CPXENVptr env;
    CPXLPptr lp;
    create_placement_lp(RATIONAL, &env, lp);

    if (solve_linear_program(env, lp, INTEGER)) return;

#else

    prob = create_placement_lp(RATIONAL);
    if (solve_linear_program(prob, RATIONAL)) return;

#endif

    srand(RANDOM_SEED);
    initialize_global_server_loads();

    // For each service, pick on which server it lands
    for (i = 0; i < flex_prob->num_services; i++) {
        total_weight = 0.0;
        for (j = 0; j < flex_prob->num_servers; j++) {
#ifdef CPLEX
            CPXgetx(env, lp, &val, LP_E_IJ_COL, LP_E_IJ_COL);
#else
            val = glp_mip_col_val(prob, LP_E_IJ_COL);
#endif
            x = MAX(val, min_weight);
            if (service_can_fit_on_server_fast(i, j) && x > 0.0) {
                weights[j] = x;
                total_weight += x;
            }
        }

        // If nobody works, forget it
        if (total_weight <= 0.0) return;

        // Pick a probability
        select_over_weight = total_weight * (rand() / (RAND_MAX + 1.0));
        x = 0.0;
        for (j = 0; j < flex_prob->num_servers; j++) {
            x += weights[j];
            // need to break to skip unfeasible servers...
            if (x >= select_over_weight) break;
        }

        if (j >= flex_prob->num_servers) return;

        // set the mappings appropriately
        flex_soln->mapping[i] = j;
        add_service_load_to_server(i, j);


    }

    free_global_server_loads();
    flex_soln->success = 1;
    maximize_minimum_then_average_yield(flex_soln);
    return;
}

flexsched_solution RRND_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution("RRND");
    LPROUNDING_scheduler(flex_soln, 0.0);
    return flex_soln;
}

flexsched_solution RRNZ_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution("RRNZ");
    LPROUNDING_scheduler(flex_soln, 0.01);
    return flex_soln;
}

#ifdef CPLEX
// FIXME: obviously not yet implemented...
void maximize_average_yield_given_minimum(
    flexsched_solution flex_soln, float minyield)
{
    int i;

    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->scaled_yields[i] = minyield;
    }

    return;
}
#else
glp_prob *create_allocation_lp(flexsched_solution flex_soln, float minyield) 
{
    glp_prob *prob;
    int i, j, k;
    int ia[(flex_prob->num_fluid+1)*flex_prob->num_services+2];
    int ja[(flex_prob->num_fluid+1)*flex_prob->num_services+2];
    double ra[(flex_prob->num_fluid+1)*flex_prob->num_services+2];
    double rowmin;

    // Create GLPK problem
    prob = glp_create_prob();

    glp_add_cols(prob, flex_prob->num_services + 1);
    glp_add_rows(prob, flex_prob->num_servers * flex_prob->num_fluid + 1);

#ifdef DEBUG
    glp_set_prob_name(prob, "task allocation problem");
    glp_set_obj_name(prob, "average yield");
#endif

    // set fluid capacity constraints...
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            glp_set_row_bnds(prob, flex_prob->num_fluid*j+k+1, 
                GLP_UP, LP_IGNORE, flex_prob->fluid_capacities[j][k] - EPSILON);
        }
    }

    // fill in variables
    for (i = 0; i < flex_prob->num_services; i++) {
#ifdef DEBUG
        char buffer[10];
        sprintf(buffer, "y_{%d}", i);
        glp_set_col_name(prob, i+1, buffer);
#endif
        glp_set_col_kind(prob, i+1, GLP_CV);
        glp_set_col_bnds(prob, i+1, GLP_DB, 
            flex_prob->slas[i] + minyield * (1.0 - flex_prob->slas[i]), 1.0);
        j = flex_soln->mapping[i];
        for (k = 0; k < flex_prob->num_fluid; k++) {
            ia[flex_prob->num_fluid*i+k+1] = flex_prob->num_fluid*j+k+1;
            ja[flex_prob->num_fluid*i+k+1] = i+1;
            ra[flex_prob->num_fluid*i+k+1] = flex_prob->fluid_needs[i][k];
        }
    }

    rowmin = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        ia[flex_prob->num_services*flex_prob->num_fluid+i+1] = 
            flex_prob->num_servers*flex_prob->num_fluid+1;
        ja[flex_prob->num_services*flex_prob->num_fluid+i+1] = i+1; 
        ra[flex_prob->num_services*flex_prob->num_fluid+i+1] = 
            1.0 / (1.0 - flex_prob->slas[i]);
        rowmin += flex_prob->slas[i] / (1.0 - flex_prob->slas[i]);
    }

    glp_set_row_bnds(prob, flex_prob->num_servers * flex_prob->num_fluid + 1,
        GLP_LO, rowmin, LP_IGNORE);

    ia[(flex_prob->num_fluid+1)*flex_prob->num_services+1] =
        flex_prob->num_servers*flex_prob->num_fluid+1;
    ja[(flex_prob->num_fluid+1)*flex_prob->num_services+1] = 
        flex_prob->num_services+1;
    ra[(flex_prob->num_fluid+1)*flex_prob->num_services+1] = 
        -1.0*flex_prob->num_services;

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_name(prob, flex_prob->num_services+1, "Y");
    glp_set_col_bnds(prob, flex_prob->num_services+1, GLP_DB, minyield, 1.0);
    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, flex_prob->num_services+1, 1.0);

    glp_load_matrix(prob, (flex_prob->num_fluid+1)*flex_prob->num_services+1, 
        ia, ja, ra);

    return prob;

}

void maximize_average_yield_given_minimum(
    flexsched_solution flex_soln, float minyield)
{
    int i;

    glp_prob *prob = create_allocation_lp(flex_soln, minyield);

    if (solve_linear_program(prob, RATIONAL)) return;

    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->scaled_yields[i] = glp_get_col_prim(prob, i+1);
    }

    return;
}
#endif
