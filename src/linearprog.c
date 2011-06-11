#include "flexsched.h"
#include "linearprog.h"

// Number of variables (called "columns"): 
//  e_ij = H*(i-1)+j: 1 .. NH
//  y_ih = N*H+(i-1)+j: NH+1 .. 2NH
//  Y: 2NH+1
#define PLACEMENT_LP_NUM_COLS \
    (2*(flex_prob->num_services)*(flex_prob->num_servers)+1)

#define PLACEMENT_LP_COL_E_IJ \
    ((flex_prob->num_servers)*i+j)
#define PLACEMENT_LP_COL_Y_IJ \
    ((flex_prob->num_servers)*((flex_prob->num_services)+i)+j)
#define PLACEMENT_LP_COL_OBJ \
    (2*(flex_prob->num_services)*(flex_prob->num_servers))

// Number of rows  
// (A) for all J: sum_N e_{j,n} = 1                
// (B) for all J: sum_N y_{j,n} >= Y
// (C) for all J,N: y_{j,n} <= e_{j,n}
// (D) for all J,N,R: e_{j,n}*a_{j,r}+y_{j,n}*b_{j,r} <= u_{n,r}
// (E) for all N,R: sum_J e_{j,n}*c_{j,r} + y_{j,n}*b_{j,r} <= t_{n,r}
#define PLACEMENT_LP_NUM_ROWS (2*flex_prob->num_services+\
    flex_prob->num_services*flex_prob->num_servers+\
    flex_prob->num_services*flex_prob->num_servers*flex_prob->num_resources+\
    flex_prob->num_servers*flex_prob->num_resources)

// non-zero matrix elements
// (A) J rows with N elements (J*N)
// (B) J rows with N elements (J*N)
// (C) J*N rows               (J*N)
// (D) J*N*R rows             (J*N*R)
// (E) N*R rows of J elements (J*N*R)
#define PLACEMENT_LP_NUM_ELTS (flex_prob->num_services*flex_prob->num_servers*(\
    3+2*flex_prob->num_resources))

linear_program_t create_placement_lp(
    flexsched_problem_t flex_prob, int rational)
{
    int ia[PLACEMENT_LP_NUM_ELTS];
    int ja[PLACEMENT_LP_NUM_ELTS];
    double ra[PLACEMENT_LP_NUM_ELTS];

    linear_program_t lp = new_linear_program(PLACEMENT_LP_NUM_COLS,
        PLACEMENT_LP_NUM_ROWS);

    int i, j, k;
    int row, elt;

    // Add set column types and bounds for all e_ij and y_ij
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            set_col_params(lp, PLACEMENT_LP_COL_E_IJ, rational, 0.0, 1.0, 0.0);
            set_col_params(lp, PLACEMENT_LP_COL_Y_IJ, RATIONAL, 0.0, 1.0, 0.0);
        }
    }

    // objective column
    set_col_params(lp, PLACEMENT_LP_COL_OBJ, RATIONAL, 0.0, 1.0, 1.0);

    // initialize matrix row and element
    elt = 0;
    row = 0;

    // (A) for all J: sum_N e_{j,n} = 1                
    for (i = 0; i < flex_prob->num_services; i++) {
        set_row_params(lp, row, 1.0, 1.0);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = 1.0; elt++;
        }
        row++;
    }

    // (B) for all J: sum_N y_{j,n} >= Y <==> sum_N y_{j,n} - Y >= 0
    for (i = 0; i < flex_prob->num_services; i++) {
        set_row_params(lp, row, 0.0, LP_IGNORE);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = 1.0; elt++;
        }
        ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_OBJ;
        ra[elt] = -1.0;
        elt++; row++;
    }

    // (C) for all J,N: y_{j,n} <= e_{j,n} <==> e_{j,n} - y_{j,n} >= 0
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            set_row_params(lp, row, 0.0, LP_IGNORE);
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = 1.0; elt++;
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    // (D) for all J,N,R: e_{j,n}*a_{j,r}+y_{j,n}*b_{j,r} <= u_{n,r}
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            for (k = 0; k < flex_prob->num_resources; k++) {
                set_row_params(lp, row, LP_IGNORE, 
                    flex_prob->servers[j]->unit_capacities[k]);
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
                ra[elt] = flex_prob->services[i]->unit_rigid_requirements[k]; 
                elt++;
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
                ra[elt] = flex_prob->services[i]->unit_fluid_needs[k]; elt++;
                row++;
            }
        }
    }

    // (E) for all N,R: sum_J e_{j,n}*c_{j,r} + y_{j,n}*b_{j,r} <= t_{n,r}
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_resources; k++) {
            set_row_params(lp, row, LP_IGNORE, 
                flex_prob->servers[j]->total_capacities[k]);
            for (i = 0; i < flex_prob->num_services; i++) {
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
                ra[elt] = flex_prob->services[i]->total_rigid_requirements[k]; 
                elt++;
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
                ra[elt] = flex_prob->services[i]->total_fluid_needs[k]; elt++;
            }
            row++;
        }
    }

    // now load matrix...
    load_matrix(lp, PLACEMENT_LP_NUM_ELTS, ia, ja, ra);

    return lp;
}

flexsched_solution_t LPBOUND_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    linear_program_t lp = create_placement_lp(flex_prob, RATIONAL);
    double objval;
    int i;

    if (solve_linear_program(lp, RATIONAL)) return flex_soln;
    flex_soln->success = 1;
    objval = get_obj_val(lp);
    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->yields[i] = objval;
    }
    return flex_soln;
}

flexsched_solution_t MILP_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    linear_program_t lp = create_placement_lp(flex_prob, INTEGER);
    int i, j;
    int status;
    double val;

    status = solve_linear_program(lp, INTEGER);

    if (status == 2) {
        strcat(flex_soln->misc_output, "T");
    }

    if (status) return flex_soln;

    flex_soln->success = 1;

    // Retrieve the mapping
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            if(get_mip_col_val(lp, PLACEMENT_LP_COL_E_IJ)) {
                flex_soln->mapping[i] = j;
                flex_soln->yields[i] = get_col_val(lp, PLACEMENT_LP_COL_Y_IJ);
                break;
            }
        }
    }

    return flex_soln;
}

/*
 *      - RRND: LPROUNDING(0.0)
 *      - RRNZ: LPROUNDING(xxx)
 */
flexsched_solution_t LPROUNDING_solver(
    flexsched_problem_t flex_prob, float min_weight)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    int i, j;
    float weights[flex_prob->num_servers];
    float x, total_weight, select_over_weight;
    int feasible_servers[flex_prob->num_servers];
    float yields[flex_prob->num_servers];
    int num_feasible_servers;

    // Create the placement problem in rational mode
    // and solve it to compute rational mappings,
    // adding an epsilon to zero values
    linear_program_t lp = create_placement_lp(flex_prob, RATIONAL);

    if (solve_linear_program(lp, RATIONAL)) return;

    srand(RANDOM_SEED);
    initialize_global_resource_availabilities_and_loads(flex_prob);

    // For each service, pick on which server it lands
    for (i = 0; i < flex_prob->num_services; i++) {
        total_weight = 0.0;
        num_feasible_servers = 0;
        for (j = 0; j < flex_prob->num_servers; j++) {
            x = MAX(get_col_val(lp, PLACEMENT_LP_COL_E_IJ), min_weight);
            if (service_can_fit_on_server_fast(flex_soln, i, j) && x > 0.0) {
                feasible_servers[num_feasible_servers] = j;
                weights[num_feasible_servers] = x;
                yields[num_feasible_servers] = 
                    get_col_val(lp, PLACEMENT_LP_COL_Y_IJ);
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
        put_service_on_server_fast(flex_soln, i, feasible_servers[j]);
        flex_soln->yields[i] = yields[j];
    }

    free_global_resource_availabilities_and_loads(flex_prob);
    flex_soln->success = 1;
    return;
}

flexsched_solution_t LPROUNDING_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options) 
{
    flexsched_solution_t flex_soln = NULL;

    if (!strcmp(name, "RRND")) {
        flex_soln = LPROUNDING_solver(flex_prob, 0.0);
    } else if (!strcmp(name, "RRNZ")) {
        flex_soln = LPROUNDING_solver(flex_prob, 0.01);
    } else {
        fprintf(stderr, "LPROUNDING scheduler doesn't recognize %s!\n", name);
        flex_soln = new_flexsched_solution(flex_prob);
    }

    return flex_soln;
}

#if 0

#define ALLOC_LP_NUM_COLS (flex_prob->num_services+1)
#define ALLOC_LP_Y_I_COL i
#define ALLOC_LP_OBJ_COL (flex_prob->num_services)

#define ALLOC_LP_NUM_ROWS (flex_prob->num_servers*flex_prob->num_fluid+1)
#define ALLOC_LP_NUM_ELTS ((flex_prob->num_fluid+1)*flex_prob->num_services+1)

#define ALLOC_LP_ROW_JK (flex_prob->num_fluid*j+k)
#define ALLOC_LP_OBJ_ROW (flex_prob->num_servers*flex_prob->num_fluid)

linear_program_t create_allocation_lp(
    flexsched_solution flex_soln, float minyield)
{
    int ia[ALLOC_LP_NUM_ELTS];
    int ja[ALLOC_LP_NUM_ELTS];
    double ra[ALLOC_LP_NUM_ELTS];

    int i, j, k;
    int row, elt;

    linear_program_t lp = new_linear_program(ALLOC_LP_NUM_COLS, 
        ALLOC_LP_NUM_ROWS);

    // set fluid capacity constraints...
    type = 'L';
    for (j = 0; j < flex_prob->num_servers; j++) {
        for (k = 0; k < flex_prob->num_fluid; k++) {
            set_row_params(lp, row, LP_IGNORE, flex_soln->prob)
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
    flexsched_solution_t flex_soln, float minyield)
{
#if 0
    int i;
    double val;
    linear_program_t lp = create_allocation_lp(flex_soln, minyield);

    for (i = 0; i < flex_prob->num_services; i++) {
        flex_soln->yields[i] = get_col_val(lp, ALLOC_LP_COL_Y_I);
    }
#endif
    return;
}
