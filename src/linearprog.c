#include "flexsched.h"

#ifdef CPLEX

#include <ilcplex/cplex.h>

// CPLEX does zero indexed rows/cols, unlike GLPK
#define PLACEMENT_LP_E_IJ_COL \
    ((flex_prob->num_servers)*i+j)
#define PLACEMENT_LP_Y_IJ_COL \
    ((flex_prob->num_servers)*((flex_prob->num_services)+i)+j)
#define PLACEMENT_LP_OBJ_COL  \
    (2*(flex_prob->num_services)*(flex_prob->num_servers))

typedef struct linear_program_s {
    CPXENVptr *env;
    CPXLPptr *lp;
} *linear_program_t;

// for now lps always maximize...
linear_program_t new_linear_program(void)
{
    int status = 0;
    linear_program_t lp = malloc(sizeof(linear_program_s));
    lp->env = CPXopenCPLEX(&status);
    lp->lp = CPXcreateprob(lp->env, &status, "")
    CPXchgobjsen(lp->env, lp->lp, CPX_MAX);
}

void set_column_params(linear_program_t lp, int col, 
    int rational, double lb, double ub, double obj)
{
    char type = rational ? (CPX_CONTINUOUS) : (CPX_BINARY);
    ub -= lb; // ub is actually range, which is the amount allowed over the lb
    CPXnewcols(lp->env, lp->lp, 1, &obj, &lb, &ub, &type, NULL);
    return;
}

void set_row_params(
    linear_program_t lp, int row, double bound, int xxx)
{
    type = 'E'; bound = 1.0;
    CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
}

void load_matrix(linear_program_t lp, int elts, int ia[], int ja[], float ra[])
{
    CPXchgcoeflist(lp->env, lp->lp, elts, ia, ja, ra);
}

void free_linear_program(linear_program_t lp)
{
}

int solve_linear_program(linear_program_t lp, int rational) {
{
    int status;

    if (rational) {
        status = CPXlpopt(lp->env, lp->lp);
    } else {
        status = CPXmipopt(lp->env, lp->lp);
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

typedef glp_prob *linear_program_t;

// for now lps always maximize...
// fixme... parames for rows/cols?
linear_program_t new_linear_program(void)
{
    linear_program_t lp = glp_create_prob();
    glp_add_cols(lp, PLACEMENT_LP_NUM_COLS);
    glp_add_rows(lp, PLACEMENT_LP_NUM_ROWS);
    glp_set_obj_dir(lp, GLP_MAX);
    return lp;
}

void set_column_params(linear_program_t lp, int col, 
    int rational, double lb, double ub, double obj)
{
    glp_set_col_kind(lp, col, rational ? (GLP_CV) : (GLP_BV));
    glp_set_col_bnds(lp, col, GLP_DB, lb, ub);
    glp_set_obj_coef(lp, col, obj);
}

void set_row_params(
    linear_program_t lp, int row, double bound, int xxx)
{
}

void load_matrix(linear_program_t lp, int elts, int ia[], int ja[], float ra[])
{
}

void free_linear_program(linear_program_t lp)
{
}

int glp_stderr_null_out(void *info, const char *s) 
{
    return 1;
}

int solve_linear_program(linear_program_t lp, int rational)
{
    int solver_status, solution_status;

    // kill lp output
    glp_term_hook(&glp_stderr_null_out, NULL);

    if (rational) {
        solver_status = glp_simplex(prob, NULL);
        solution_status = (glp_get_status(prob) != GLP_OPT);
    } else {
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        parm.tm_lim = GLPK_TIME_LIMIT;
        solver_status = glp_intopt(prob, &parm);
        solution_status = (glp_mip_status(prob) != GLP_OPT);
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

#endif

/***********************************************************/
/* Schedulers based on the new linear program              */
/***********************************************************/

#define RATIONAL 1
#define INTEGER 0

#define LP_IGNORE 666

// Number of variables (called "columns"): 
//  e_ij = H*(i-1)+j: 1 .. NH
//  y_ih = N*H+(i-1)+j: NH+1 .. 2NH
//  Y: 2NH+1
#define PLACEMENT_LP_NUM_COLS \
    (2*(flex_prob->num_services)*(flex_prob->num_servers)+1)

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
    int i, j, k;

    /********** START CPLEX ***************/
    int row = 0, elt = 0;

    int ia[PLACEMENT_LP_NUM_ELTS];\
    int ja[PLACEMENT_LP_NUM_ELTS];\
    double ra[PLACEMENT_LP_NUM_ELTS];

    /********** END CPLEX ***************/

    // Add e_ij columns
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            type = rational ? (CPX_CONTINUOUS) : (CPX_BINARY);
            CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);
        }
    }

    // Add y_ij columns
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            type = CPX_CONTINUOUS;
            CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);
        }
    }

    // objective column
    obj = 1.0;
    type = rational ? (CPX_CONTINUOUS) : (CPX_BINARY);
    CPXnewcols(env, lp, 1, &obj, &bound, &range, &type, NULL);

    // (A) for all J: sum_N e_{j,n} = 1                
    for (i = 0; i < flex_prob->num_services; i++) {
        type = 'E'; bound = 1.0;
        CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = 1.0; elt++;
        }
        row++;
    }

    // (B) for all J: sum_N y_{j,n} >= Y
    for (i = 0; i < flex_prob->num_services; i++) {
        type = 'L'; bound = 0.0;
        CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
        for (j = 0; j < flex_prob->num_servers; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = -1.0; elt++;
        }
        ia[elt] = row; ja[elt] = PLACEMENT_LP_OBJ_COL;
        ra[elt] = flex_prob->slas[i] 1.0;
        elt++; row++;
    }

    // (C) for all J,N: y_{j,n} <= e_{j,n}
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            type = 'L'; bound = 0.0;
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = 1.0; elt++;
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    // (D) for all J,N,R: e_{j,n}*a_{j,r}+y_{j,n}*b_{j,r} <= u_{n,r}
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            for (k = 0; k < flex_prob->num_resources; k++) {
                type = 'L';
                bound = flex_prob->servers[j]->unit_capacities[k];
                CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
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
            type = 'L';
            bound = flex_prob->servers[j]->total_capacities[k];
            CPXnewrows(env, lp, 1, &bound, &type, NULL, NULL);
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
    CPXchgcoeflist (env, lp, PLACEMENT_LP_NUM_ELTS, ia, ja, ra);

    return blah;
}

linear_program_t create_placement_lp(
    flexsched_problem_t flex_prob, int rational) 
{
    int i, j, k;

    /********** START GLPK ***************/
    int row = 1, elt = 1;
    int ia[PLACEMENT_LP_NUM_ELTS+1];    
    int ja[PLACEMENT_LP_NUM_ELTS+1];   
    double ra[PLACEMENT_LP_NUM_ELTS+1];

    /*********** END GLPK *****************/

    // Set bounds for the min yield: between 0 and 1.0
    glp_set_col_bnds(prob, LP_OBJ_COL, GLP_DB, 0.0, 1.0);
    glp_set_obj_coef(prob, LP_OBJ_COL, 1.0);

    // Set bounds for the e_ij: binary if !rational
    for (i = 0; i < flex_prob->num_services; i++) {
        for (j = 0; j < flex_prob->num_servers; j++) {
            if (rational) {
                glp_set_col_kind(prob, LP_E_IJ_COL, GLP_CV);
            } else {
                glp_set_col_kind(prob, LP_E_IJ_COL, GLP_BV);
            }
            glp_set_col_bnds(prob, LP_E_IJ_COL, GLP_DB, 0.0, 1.0);
            glp_set_col_bnds(prob, LP_Y_IJ_COL, GLP_DB, 0.0, 1.0);
        }
    }

    glp_set_row_bnds(prob, row, GLP_FX, 1.0, 1.0);
    // upper bound 0...
    glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 0.0);
    // lower bound 0...
    glp_set_row_bnds(prob, row, GLP_UP, LP_IGNORE, 0.0);

    // load matrix...
    glp_load_matrix(prob, LP_NUM_ELTS, ia, ja, ra);

}
flexsched_solution_t LPBOUND_scheduler(char *name, char **options)
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
