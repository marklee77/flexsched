#include "flexsched.h"

/***********************************************************/
/* Schedulers based on the linear program described in the */
/* JPDC paper. Chekuri uses a different linear program and */
/* is contained in its own file                            */
/***********************************************************/

#define RATIONAL 1
#define INTEGER 0

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

// CPLEX does zero indexed rows/cols, unlike GLPK
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
    double obj, bound, range;
    int row, elt;
    int ia[LP_NUM_ELTS];    // array of matrix element row indicies
    int ja[LP_NUM_ELTS];    // array of matrix element column indicies
    double ra[LP_NUM_ELTS]; // array of matrix element coefficients
    char type;

    *retenv = NULL;
    *retlp = NULL; 

    // create CPLEX problem
    env = CPXopenCPLEX (&status);
    if (NULL == env) {
        fprintf (stderr, "Could not open CPLEX environment.\n");
        return 1;
    }
  
#ifdef DEBUG
    // turn on output
    status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
        fprintf (stderr,
                 "Failure to turn on screen indicator, error %d.\n", status);
        return 1;
    }
    // Turn on data checking
    status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
    if ( status ) {
        fprintf (stderr, 
            "Failure to turn on data checking, error %d.\n", status);
        return 1; 
    }
#endif

    // create the problem
    lp = CPXcreateprob(env, &status, "task placement problem");

    if ( lp == NULL ) {
        fprintf (stderr, "Failed to create LP.\n");
        return 1;
    }  

    // default objective coefficient and upper/lower bounds for most variables
    obj = 0.0;
    bound = 0.0;
    range = 1.0;

    // Add e_ij columns
    if (rational) {
        type = CPX_CONTINUOUS;
    } else {
        type = CPX_BINARY;
    }
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);
        }
    }

    // Add y_ij columns
    type = CPX_CONTINUOUS;
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);
        }
    }

    // objective column
    obj = 1.0;
    CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);

    elt = 0; // element counter
    row = 0; // row counter

    //  (A) for all i: sum_j e_ij = 1
    type = 'E';
    bound = 1.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    //  (B) for all i, j:  y_ih <= e_ij <=> y_ij - e_ij <= 0.0
    type = 'L';
    bound = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0;  elt++;
            ia[elt] = row; ja[elt] = LP_E_IJ_COL; ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    type = 'G';
    //  (C) for all i: sum_h y_ih >= sla_i
    for (i = 0; i < flex_prob->num_services; i++) {
        bound = flex_prob->slas[i];
        CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = LP_Y_IJ_COL; ra[elt] = 1.0; elt++;
        }
        row++;
    }

    type = 'L';
    //  (D) for all h, r:  sum_i r_ij*e_ih <= capacity
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_rigid; k++) {
            bound = flex_prob->rigid_capacities[j][k];
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
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
            bound = flex_prob->fluid_capacities[j][k];
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = LP_Y_IJ_COL; 
                ra[elt] = flex_prob->fluid_needs[i][k]; elt++;
            }
            row++;
        }
    }

    type = 'G';
    //  (F) for all i: (sum_h y_ih - sla_i)/(1-sla_i) >= Y <=>
    //      sum_h y_ih + Y (sla_i - 1) >= sla_i 
    for (i = 0; i < flex_prob->num_services; i++) {
        bound = flex_prob->slas[i];
        CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
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

#define GLPK_TIME_LIMIT 10*60*1000

// GLPK does not use zero indexing, which seems like a waste.
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

flexsched_solution LPBOUND_scheduler(char *name, char **options)
{
    flexsched_solution flex_soln = new_flexsched_solution();
    double objval;

// create and solve LP, get objective value
#ifdef CPLEX
    CPXENVptr env;
    CPXLPptr lp;
    create_placement_lp(RATIONAL, &env, &lp);
    if (solve_linear_program(env, lp, RATIONAL)) return flex_soln;
    CPXgetobjval(env, lp, &objval);
#else
    glp_prob *prob = create_placement_lp(RATIONAL);
    if (solve_linear_program(prob, RATIONAL)) return flex_soln;
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

flexsched_solution MILP_scheduler(char *name, char **options)
{
    int i, j;
    int status;
    double val;
    flexsched_solution flex_soln = new_flexsched_solution();

#ifdef CPLEX
    CPXENVptr env;
    CPXLPptr lp;

    create_placement_lp(INTEGER, &env, &lp);

    if (status = solve_linear_program(env, lp, INTEGER)) {
        return flex_soln;
    }
#else
    glp_prob *prob = create_placement_lp(INTEGER);

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
            CPXgetx(env, lp, &val, LP_Y_IJ_COL, LP_Y_IJ_COL);
#else
            val = glp_mip_col_val(prob, LP_Y_IJ_COL);
#endif
                flex_soln->scaled_yields[i] = 
                    (val - flex_prob->slas[i]) / (1 - flex_prob->slas[i]);
                break;
            }
        }
    }

    flex_soln->success = 1;
    return flex_soln;
}

/*
 *      - RRND: LPROUNDING(0.0)
 *      - RRNZ: LPROUNDING(xxx)
 */
// FIXME: why do we need to call the optimizer or get fluid resource overload?
void LPROUNDING_solver(flexsched_solution flex_soln, float min_weight)
{
    int i, j;
    float weights[flex_prob->num_servers];
    float x, total_weight, select_over_weight;
    double val;
    int feasible_servers[flex_prob->num_servers];
    float yields[flex_prob->num_servers];
    int num_feasible_servers;

    // Create the placement problem in rational mode
    // and solve it to compute rational mappings,
    // adding an epsilon to zero values
#ifdef CPLEX
    CPXENVptr env;
    CPXLPptr lp;
    create_placement_lp(RATIONAL, &env, &lp);
    if (solve_linear_program(env, lp, INTEGER)) return;
#else
    glp_prob *prob = create_placement_lp(RATIONAL);
    if (solve_linear_program(prob, RATIONAL)) return;
#endif

    srand(RANDOM_SEED);
    initialize_global_server_loads();

    // For each service, pick on which server it lands
    for (i = 0; i < flex_prob->num_services; i++) {
        total_weight = 0.0;
        num_feasible_servers = 0;
        for (j = 0; j < flex_prob->num_servers; j++) {
#ifdef CPLEX
            CPXgetx(env, lp, &val, LP_E_IJ_COL, LP_E_IJ_COL);
#else
            val = glp_get_col_prim(prob, LP_E_IJ_COL);
#endif
            x = MAX(val, min_weight);
            if (service_can_fit_on_server_fast(i, j) && x > 0.0) {
                feasible_servers[num_feasible_servers] = j;
                weights[num_feasible_servers] = x;
#ifdef CPLEX
                CPXgetx(env, lp, &val, LP_Y_IJ_COL, LP_Y_IJ_COL);
#else
                val = glp_get_col_prim(prob, LP_Y_IJ_COL);
#endif
                yields[num_feasible_servers] = val;
                total_weight += x;
                num_feasible_servers++;
            }
        }

        // If nobody works, forget it
        if (!num_feasible_servers) return;

        // Pick a probability
        select_over_weight = total_weight * (rand() / (RAND_MAX + 1.0));
        x = 0.0;
        for (j = 0; j < num_feasible_servers; j++) {
            x += weights[j];
            // need to break to skip unfeasible servers...
            if (x >= select_over_weight) break;
        }

        if (j >= num_feasible_servers) return;

        // set the mappings appropriately
        flex_soln->mapping[i] = feasible_servers[j];
        flex_soln->scaled_yields[i] = yields[j];
        add_service_load_to_server(i, feasible_servers[j]);

    }

    free_global_server_loads();
    flex_soln->success = 1;
    // otherwise rounding gives bad values...
    // FIXME: seems like it *should* work without this...
    maximize_minimum_then_average_yield(flex_soln);
    return;
}

flexsched_solution LPROUNDING_scheduler(char *name, char **options) {
    flexsched_solution flex_soln = new_flexsched_solution();

    if (!strcmp(name, "RRND")) {
        LPROUNDING_solver(flex_soln, 0.0);
    } else if (!strcmp(name, "RRNZ")) {
        LPROUNDING_solver(flex_soln, 0.01);
    } else {
        fprintf(stderr, "LPROUNDING scheduler doesn't recognize %s!\n", name);
    }

    return flex_soln;
}

#define ALLOC_LP_NUM_COLS (flex_prob->num_services+1)
#define ALLOC_LP_NUM_ROWS (flex_prob->num_servers*flex_prob->num_fluid+1)
#define ALLOC_LP_NUM_ELTS ((flex_prob->num_fluid+1)*flex_prob->num_services+1)

#ifdef CPLEX
// CPLEX uses 0 index cols/rows
#define ALLOC_LP_Y_I_COL i
#define ALLOC_LP_OBJ_COL (flex_prob->num_services)
#define ALLOC_LP_ROW_JK (flex_prob->num_fluid*j+k)
#define ALLOC_LP_OBJ_ROW (flex_prob->num_servers*flex_prob->num_fluid)

int create_allocation_lp(flexsched_solution flex_soln, float minyield, 
    CPXENVptr *retenv, CPXLPptr *retlp) 
{
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int status = 0;
    int i, j, k;
    double obj, bound, range;
    char type = CPX_CONTINUOUS;
    int elt;
    int ia[ALLOC_LP_NUM_ELTS+1];
    int ja[ALLOC_LP_NUM_ELTS+1];
    double ra[ALLOC_LP_NUM_ELTS+1];
    double rowmin;

    *retenv = NULL;
    *retlp = NULL; 

    // create CPLEX problem
    env = CPXopenCPLEX (&status);
    if (NULL == env) {
        fprintf (stderr, "Could not open CPLEX environment.\n");
        return 1;
    }
  
#ifdef DEBUG
    // turn on output
    status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
        fprintf (stderr,
                 "Failure to turn on screen indicator, error %d.\n", status);
        return 1;
    }
    // Turn on data checking
    status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
    if ( status ) {
        fprintf (stderr, 
            "Failure to turn on data checking, error %d.\n", status);
        return 1; 
    }
#endif

    // create the problem
    lp = CPXcreateprob(env, &status, "task allocation problem");

    if ( lp == NULL ) {
        fprintf (stderr, "Failed to create LP.\n");
        return 1;
    }  

    // set fluid capacity constraints...
    type = 'L';
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            // not totally happy about the EPSILON, but it seems necessary to
            // keep everything within bounds
            bound = flex_prob->fluid_capacities[j][k] - EPSILON;
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
        }
    }

    elt = 1;

    for (i = 0; i < flex_prob->num_services; i++) {
    }

    // fill in variables
    // for all j,k sum_i (for i mapped to j) need_ik * y_i <= capacity_jk
    // add y_i columns
    obj = 0.0;
    range = 1.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        bound = flex_prob->slas[i] + minyield * (1.0 - flex_prob->slas[i]);
        CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);
        j = flex_soln->mapping[i];
        for (k = 0; k < flex_prob->num_fluid; k++) {
            ia[elt] = ALLOC_LP_ROW_JK; ja[elt] = ALLOC_LP_Y_I_COL; 
            ra[elt] = flex_prob->fluid_needs[i][k]; elt++;
        }
    }

    // add objective column
    obj = 1.0;
    bound = minyield;
    CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);

    rowmin = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        ia[elt] = ALLOC_LP_OBJ_ROW; ja[elt] = ALLOC_LP_Y_I_COL; 
        ra[elt] = 1.0 / (1.0 - flex_prob->slas[i]); elt++;
        rowmin += flex_prob->slas[i] / (1.0 - flex_prob->slas[i]);
    }

    // objective row
    type = 'G';
    CPXnewrows(env, lp, 1, &rowmin, &type, NULL, NULL);

    ia[elt] = ALLOC_LP_OBJ_ROW; ja[elt] = ALLOC_LP_OBJ_COL;
    ra[elt] = -1.0*flex_prob->num_services; elt++;

    status = CPXchgcoeflist (env, lp, ALLOC_LP_NUM_ELTS, ia, ja, ra);
    if (status) return 1;

    CPXchgobjsen (env, lp, CPX_MAX);

    *retenv = env;
    *retlp = lp;

    return 0;
}

#else

#define ALLOC_LP_Y_I_COL (i+1)
#define ALLOC_LP_OBJ_COL (flex_prob->num_services+1)
#define ALLOC_LP_ROW_JK (flex_prob->num_fluid*j+k+1)
#define ALLOC_LP_OBJ_ROW (flex_prob->num_servers*flex_prob->num_fluid+1)

glp_prob *create_allocation_lp(flexsched_solution flex_soln, float minyield) 
{
    glp_prob *prob;
    int i, j, k;
    int elt;
    int ia[ALLOC_LP_NUM_ELTS+1];
    int ja[ALLOC_LP_NUM_ELTS+1];
    double ra[ALLOC_LP_NUM_ELTS+1];
    double rowmin;

    // Create GLPK problem
    prob = glp_create_prob();

    glp_add_cols(prob, ALLOC_LP_NUM_COLS);
    glp_add_rows(prob, ALLOC_LP_NUM_ROWS);

#ifdef DEBUG
    glp_set_prob_name(prob, "task allocation problem");
    glp_set_obj_name(prob, "average yield");
#endif

    // set fluid capacity constraints...
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            // not totally happy about the EPSILON, but it seems necessary to
            // keep everything within bounds
            glp_set_row_bnds(prob, ALLOC_LP_ROW_JK, 
                GLP_UP, LP_IGNORE, flex_prob->fluid_capacities[j][k] - EPSILON);
        }
    }

    elt = 1;

    // fill in variables
    // for all j,k sum_i (for i mapped to j) need_ik * y_i <= capacity_jk
    for (i = 0; i < flex_prob->num_services; i++) {
#ifdef DEBUG
        char buffer[10];
        sprintf(buffer, "y_{%d}", i);
        glp_set_col_name(prob, ALLOC_LP_Y_I_COL, buffer);
#endif
        glp_set_col_kind(prob, ALLOC_LP_Y_I_COL, GLP_CV);
        glp_set_col_bnds(prob, ALLOC_LP_Y_I_COL, GLP_DB, 
            flex_prob->slas[i] + minyield * (1.0 - flex_prob->slas[i]), 1.0);
        j = flex_soln->mapping[i];
        for (k = 0; k < flex_prob->num_fluid; k++) {
            ia[elt] = ALLOC_LP_ROW_JK; ja[elt] = ALLOC_LP_Y_I_COL; 
            ra[elt] = flex_prob->fluid_needs[i][k]; elt++;
        }
    }

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_name(prob, ALLOC_LP_OBJ_COL, "Y");
    glp_set_col_bnds(prob, ALLOC_LP_OBJ_COL, GLP_DB, minyield, 1.0);

    rowmin = 0.0;
    for (i = 0; i < flex_prob->num_services; i++) {
        ia[elt] = ALLOC_LP_OBJ_ROW; ja[elt] = ALLOC_LP_Y_I_COL; 
        ra[elt] = 1.0 / (1.0 - flex_prob->slas[i]); elt++;
        rowmin += flex_prob->slas[i] / (1.0 - flex_prob->slas[i]);
    }

    glp_set_row_bnds(prob, ALLOC_LP_OBJ_ROW, GLP_LO, rowmin, LP_IGNORE);

    ia[elt] = ALLOC_LP_OBJ_ROW; ja[elt] = ALLOC_LP_OBJ_COL;
    ra[elt] = -1.0*flex_prob->num_services; elt++;

    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, ALLOC_LP_OBJ_COL, 1.0);

    glp_load_matrix(prob, ALLOC_LP_NUM_ELTS, ia, ja, ra);

    return prob;

}
#endif

void maximize_average_yield_given_minimum(
    flexsched_solution flex_soln, float minyield)
{
    int i;
    double val;

#ifdef CPLEX
    CPXENVptr env;
    CPXLPptr lp;
    create_allocation_lp(flex_soln, minyield, &env, &lp);
    if (solve_linear_program(env, lp, RATIONAL)) return;
#else
    glp_prob *prob = create_allocation_lp(flex_soln, minyield);
    if (solve_linear_program(prob, RATIONAL)) return;
#endif

    for (i = 0; i < flex_prob->num_services; i++) {
#ifdef CPLEX
        CPXgetx(env, lp, &val, ALLOC_LP_Y_I_COL, ALLOC_LP_Y_I_COL);
#else
        val = glp_get_col_prim(prob, ALLOC_LP_Y_I_COL);
#endif
        flex_soln->scaled_yields[i] = val;
    }

#ifdef CPLEX
    // FIXME: how to delete cplex lp?
#else
    glp_erase_prob(prob);
    glp_delete_prob(prob);
    glp_free_env();
#endif

    return;
}
