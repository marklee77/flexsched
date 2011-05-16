#include <glpk.h>
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
#define LP_E_IJ_COL ((flex_prob->num_servers)*i+j+1)
#define LP_Y_IJ_COL ((flex_prob->num_servers)*((flex_prob->num_services)+i)+j+1)
#define LP_OBJ_COL  (2*(flex_prob->num_services)*(flex_prob->num_servers)+1)
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

flexsched_solution LPBOUND_scheduler(
    char *ignore1, char *ignore2, char *ignore3)
{
    flexsched_solution flex_soln = new_flexsched_solution("LPBOUND");

    glp_prob *prob = create_placement_lp(RATIONAL);

    if (solve_linear_program(prob, RATIONAL)) {
        return flex_soln;
    }

    if (fabs(flex_prob->lpbound - glp_get_col_prim(prob, LP_OBJ_COL)) > EPSILON)
    {
        fprintf(stderr,"Error: Immediate bound = %.3f, LP bound = %.3f\n", 
            flex_prob->lpbound, glp_get_col_prim(prob, LP_OBJ_COL));
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
    glp_prob *prob;
    flexsched_solution flex_soln = new_flexsched_solution("MILP");

    // Create the problem (in MILP form)
    prob = create_placement_lp(INTEGER);

    // Solve it
    if (status = solve_linear_program(prob, INTEGER)) {
        if (status == 2) {
            strcat(flex_soln->misc_output, "T");
        }
        return flex_soln;
    }

    // Retrieve the mapping
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            if (glp_mip_col_val(prob, LP_E_IJ_COL) >= 1.0 - EPSILON) {
                flex_soln->mapping[i] = j;
                flex_soln->scaled_yields[i] = 
                    (glp_mip_col_val(prob, LP_Y_IJ_COL) - flex_prob->slas[i]) / 
                    (1 - flex_prob->slas[i]);
                break;
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr,"MILP solution: minimum yield = %.3f\n", 
        glp_mip_col_val(prob, LP_OBJ_COL));
#endif

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

    // Create the placement problem in rational mode
    // and solve it to compute rational mappings,
    // adding an epsilon to zero values
    prob = create_placement_lp(RATIONAL);
    if (solve_linear_program(prob, RATIONAL)) return;

    srand(RANDOM_SEED);
    initialize_global_server_loads();

    // For each service, pick on which server it lands
    for (i = 0; i < flex_prob->num_services; i++) {
        total_weight = 0.0;
        for (j = 0; j < flex_prob->num_servers; j++) {
            x = MAX(glp_get_col_prim(prob, LP_E_IJ_COL), min_weight);
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

# if 0
glp_prob *create_allocation_lp(flexsched_solution flex_soln, float minyield) 
{
    glp_prob *prob;
    int i, j, k;
    int row, elt;
    int ia[LP_NUM_ELTS+1];    // array of matrix element row indicies
    int ja[LP_NUM_ELTS+1];    // array of matrix element column indicies
    double ra[LP_NUM_ELTS+1]; // array of matrix element coefficients

    // Create GLPK problem
    prob = glp_create_prob();

    glp_add_cols(prob, flex_soln->num_services + 1);
    glp_add_rows(prob, flex_soln->num_servers * flex_soln->num_fluid + 1);

#ifdef DEBUG
    glp_set_prob_name(prob, "task allocation problem");
    glp_set_obj_name(prob, "average yield");
#endif

    // set fluid capacity constraints...
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            glp_set_row_bnds(prob, flex_prob->num_fluid*j+k+1, 
                GLP_UP, LP_IGNORE, flex_prob->fluid_capacities[j][k]);
        }
    }

    // fill in variables
    for (i = 0; i < flex_prob->num_services; i++) {
#ifdef DEBUG
        char buffer[10];
        sprintf(buffer, "y_{%d}", i);
        glp_set_col_name(prob, i, buffer);
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

    glp_set_row_bnds(prob, flex_prob->num_services * flex_prob->num_fluid + 1, 
        GLP_LO, flex_prob->num_services * minyield, LP_IGNORE);

    for (i = 0; i < flex_prob->num_services; i++) {
        ia[flex_prob->num_services * flex_prob->num_fluid + 1] = 
            flex_prob->num_services * flex_prob->num_fluid + 1;
        ja[flex_prob->num_services * flex_prob->num_fluid + 1] = i+1; 
        ra[flex_prob->num_services * flex_prob->num_fluid + 1] = 1.0;
        ia[elt] = row; ja[elt] = LP_OBJ_COL; 
        ra[elt] = flex_prob->slas[i] - 1.0;
    }

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_bnds(prob, flex_prob->num_services + 1, GLP_DB, minyield, 1.0);
    glp_set_obj_dir(prob, GLP_MAX);
    glp_set_obj_coef(prob, flex_prob->num_services + 1, 1.0);

    return prob;

}
#endif

