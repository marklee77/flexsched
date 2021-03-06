#include "flexsched.h"
#include "linearprog.h"

// Number of variables (called "columns"): 
//  e_ij = H*(i-1)+j: 1 .. NH
//  y_ih = N*H+(i-1)+j: NH+1 .. 2NH
//  Y: 2NH+1
#define PLACEMENT_LP_NUM_COLS \
    (2*(flex_prob->num_jobs)*(flex_prob->num_nodes)+1)

#define PLACEMENT_LP_COL_E_IJ \
    ((flex_prob->num_nodes)*i+j)
#define PLACEMENT_LP_COL_Y_IJ \
    ((flex_prob->num_nodes)*((flex_prob->num_jobs)+i)+j)
#define PLACEMENT_LP_COL_OBJ \
    (2*(flex_prob->num_jobs)*(flex_prob->num_nodes))

// Number of rows  
// (A) for all J: sum_N e_{j,n} = 1                
// (B) for all J: sum_N y_{j,n} >= Y
// (C) for all J,N: y_{j,n} <= e_{j,n}
// (D) for all J,N,R: e_{j,n}*a_{j,r}+y_{j,n}*b_{j,r} <= u_{n,r}
// (E) for all N,R: sum_J e_{j,n}*c_{j,r} + y_{j,n}*b_{j,r} <= t_{n,r}
#define PLACEMENT_LP_NUM_ROWS (2*flex_prob->num_jobs+\
    flex_prob->num_jobs*flex_prob->num_nodes+\
    flex_prob->num_jobs*flex_prob->num_nodes*flex_prob->num_resources+\
    flex_prob->num_nodes*flex_prob->num_resources)

// non-zero matrix elements
// (A) J rows with N elements     (J*N)
// (B) J rows with N+1 elements   (J*(N+1))
// (C) J*N rows with 2 elements   (2*J*N)
// (D) J*N*R rows with 2 elements (2*J*N*R)
// (E) N*R rows of 2*J elements   (2*J*N*R)
// total: J*(1+4*N*(R+1))
#define PLACEMENT_LP_NUM_ELTS (flex_prob->num_jobs*\
    (1+4*flex_prob->num_nodes*(1+flex_prob->num_resources)))
        
linear_program_t create_placement_lp(
    flexsched_problem_t flex_prob, int rational)
{
    int *ia, *ja;
    double *ra;
        
    ia = calloc(PLACEMENT_LP_NUM_ELTS, sizeof(int));
    ja = calloc(PLACEMENT_LP_NUM_ELTS, sizeof(int));
    ra = calloc(PLACEMENT_LP_NUM_ELTS, sizeof(double));

    linear_program_t lp = new_linear_program(PLACEMENT_LP_NUM_COLS,
        PLACEMENT_LP_NUM_ROWS);

    int i, j, k;
    int row, elt;

    // Add set column types and bounds for all e_ij and y_ij
    for (i = 0; i < flex_prob->num_jobs; i++) {
        for (j = 0; j < flex_prob->num_nodes; j++) {
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
    for (i = 0; i < flex_prob->num_jobs; i++) {
        set_row_params(lp, row, 1.0, 1.0);
        for (j = 0; j < flex_prob->num_nodes; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = 1.0; elt++;
        }
        row++;
    }

    // (B) for all J: sum_N y_{j,n} >= Y <==> sum_N y_{j,n} - Y >= 0
    for (i = 0; i < flex_prob->num_jobs; i++) {
        set_row_params(lp, row, 0.0, LP_IGNORE);
        for (j = 0; j < flex_prob->num_nodes; j++) {
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = 1.0; elt++;
        }
        ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_OBJ;
        ra[elt] = -1.0;
        elt++; row++;
    }

    // (C) for all J,N: y_{j,n} <= e_{j,n} <==> e_{j,n} - y_{j,n} >= 0
    for (i = 0; i < flex_prob->num_jobs; i++) {
        for (j = 0; j < flex_prob->num_nodes; j++) {
            set_row_params(lp, row, 0.0, LP_IGNORE);
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
            ra[elt] = 1.0; elt++;
            ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
            ra[elt] = -1.0; elt++;
            row++;
        }
    }
   
    // (D) for all J,N,R: e_{j,n}*a_{j,r}+y_{j,n}*b_{j,r} <= u_{n,r}
    for (i = 0; i < flex_prob->num_jobs; i++) {
        for (j = 0; j < flex_prob->num_nodes; j++) {
            for (k = 0; k < flex_prob->num_resources; k++) {
                set_row_params(lp, row, LP_IGNORE,
                    flex_prob->nodes[j]->unit_capacities[k] - EPSILON);
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
                ra[elt] = flex_prob->jobs[i]->unit_rigid_requirements[k]; 
                elt++;
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
                ra[elt] = flex_prob->jobs[i]->unit_fluid_needs[k]; elt++;
                row++;
            }
        }
    }

    // (E) for all N,R: sum_J e_{j,n}*c_{j,r} + y_{j,n}*b_{j,r} <= t_{n,r}
    for (j = 0; j < flex_prob->num_nodes; j++) {
        for (k = 0; k < flex_prob->num_resources; k++) {
            set_row_params(lp, row, LP_IGNORE, 
                flex_prob->nodes[j]->total_capacities[k] - EPSILON);
            for (i = 0; i < flex_prob->num_jobs; i++) {
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_E_IJ; 
                ra[elt] = flex_prob->jobs[i]->total_rigid_requirements[k]; 
                elt++;
                ia[elt] = row; ja[elt] = PLACEMENT_LP_COL_Y_IJ; 
                ra[elt] = flex_prob->jobs[i]->total_fluid_needs[k]; elt++;
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
    for (i = 0; i < flex_prob->num_jobs; i++) {
        flex_soln->yields[i] = objval - EPSILON;
    }
    free_linear_program(lp);
    return flex_soln;
}

flexsched_solution_t MILP_scheduler(
    flexsched_problem_t flex_prob, char *name, char **options)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    linear_program_t lp = create_placement_lp(flex_prob, INTEGER);
    int i, j;
    int status = 1;
    double val;

#ifdef CPLEX
    strcpy(flex_soln->misc_output, "CPLEX");
#else
    strcpy(flex_soln->misc_output, "GLPK");
#endif

    status = solve_linear_program(lp, INTEGER);

    if (status == 2) {
        strcat(flex_soln->misc_output, " T");
    }

    if (!status) {
        flex_soln->success = 1;
        // Retrieve the mapping
        for (i = 0; i < flex_prob->num_jobs; i++) {
            for (j = 0; j < flex_prob->num_nodes; j++) {
                if(get_mip_col_val(lp, PLACEMENT_LP_COL_E_IJ)) {
                    flex_soln->mapping[i] = j;
                    flex_soln->yields[i] = 
                        get_col_val(lp, PLACEMENT_LP_COL_Y_IJ) - EPSILON;
                    break;
                }
            }
        }
    }

    free_linear_program(lp);
    return flex_soln;
}

/*
 *      - RRND: LPROUNDING(0.0)
 *      - RRNZ: LPROUNDING(xxx)
 *
 * FIXME: bad flow, kind of spaghetti-like
 */
flexsched_solution_t LPROUNDING_solver(
    flexsched_problem_t flex_prob, double min_weight)
{
    flexsched_solution_t flex_soln = new_flexsched_solution(flex_prob);
    int i, j, k;
    double weights[flex_prob->num_nodes];
    double x, total_weight, select_over_weight;
    int feasible_nodes[flex_prob->num_nodes];
    double yields[flex_prob->num_nodes];
    int num_feasible_nodes;
    int repeat, couldreduce;

    // Create the placement problem in rational mode
    // and solve it to compute rational mappings,
    // adding an epsilon to zero values
    linear_program_t lp = create_placement_lp(flex_prob, RATIONAL);

    if (solve_linear_program(lp, RATIONAL)) {
        free_linear_program(lp);
        return flex_soln;
    }

    srand(RANDOM_SEED);
    initialize_global_resource_availabilities_and_loads(flex_prob);

    // For each job, pick on which node it lands
    for (i = 0; i < flex_prob->num_jobs; i++) {
        total_weight = 0.0;
        num_feasible_nodes = 0;
        for (j = 0; j < flex_prob->num_nodes; j++) {
            x = MAX(get_col_val(lp, PLACEMENT_LP_COL_E_IJ), min_weight);
            if (job_can_fit_on_node_fast(flex_soln, i, j) && x > 0.0) {
                feasible_nodes[num_feasible_nodes] = j;
                weights[num_feasible_nodes] = x;
                yields[num_feasible_nodes] = MAX(EPSILON, 
                    get_col_val(lp, PLACEMENT_LP_COL_Y_IJ) - EPSILON);
                total_weight += x;
                num_feasible_nodes++;
            }
        }

        // If nobody works, forget it
        if (!num_feasible_nodes) {
            free_linear_program(lp);
            return flex_soln;
        }

        // Pick a probability
        select_over_weight = total_weight * (rand() / (RAND_MAX + 1.0));
        x = 0.0;
        for (j = 0; j < num_feasible_nodes; j++) {
            x += weights[j];
            // need to break to skip unfeasible nodes...
            if (x >= select_over_weight) break;
        }

        if (j >= num_feasible_nodes) {
            free_linear_program(lp);
            return flex_soln;
        }

        // set the mappings appropriately
        put_job_on_node_fast(flex_soln, i, feasible_nodes[j]);
        flex_soln->yields[i] = MAX(EPSILON, yields[j]);
    }
    free_linear_program(lp);
    free_global_resource_availabilities_and_loads(flex_prob);

    // a little hacky, but adjust the final solution...
    do {
        repeat = 0;
        for (j = 0; j < flex_prob->num_nodes; j++) {
            for (k = 0; k < flex_prob->num_resources; k++) {
                if (compute_allocated_resource(flex_soln, j, k) >
                        flex_prob->nodes[j]->total_capacities[k]) {
                    repeat = 1;
                    couldreduce = 0;
                    for (i = 0; i < flex_prob->num_jobs; i++) {
                        if (flex_soln->mapping[i] != j || flex_soln->yields[i]
                                <= EPSILON) continue;
                        couldreduce = 1;
                        flex_soln->yields[i] = MAX(EPSILON, flex_soln->yields[i] - EPSILON);
                    } 
                    if (!couldreduce) return flex_soln; // can't fix problem
                }
            }
        }
    } while (repeat);

    flex_soln->success = 1;
    return flex_soln;
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

#define ALLOC_LP_NUM_COLS (flex_soln->prob->num_jobs)
#define ALLOC_LP_COL_Y_I i

// (A) for all J,R: y_{j}*b_{j,r} <= u_{n,r} - a_{j,r}
// (B) for all N,R: sum_J y_{j}*b_{j,r} <= t_{n,r} - sum_J c_{j,r}
// total rows: J*R+N*R = (J+N)*R, total elts: J*R+J*R = 2*J*R
#define ALLOC_LP_NUM_ROWS ((flex_soln->prob->num_jobs+\
    flex_soln->prob->num_nodes)*flex_soln->prob->num_resources)
#define ALLOC_LP_NUM_ELTS \
    (2*flex_soln->prob->num_jobs*flex_soln->prob->num_resources)
#define ALLOC_LP_ROW_B_JK \
    (flex_soln->prob->num_resources*(flex_soln->prob->num_jobs+j)+k)

linear_program_t create_allocation_lp(
    flexsched_solution_t flex_soln, double minyield)
{
    int ia[ALLOC_LP_NUM_ELTS];
    int ja[ALLOC_LP_NUM_ELTS];
    double ra[ALLOC_LP_NUM_ELTS];

    int i, j, k;
    int row, elt;

    linear_program_t lp = new_linear_program(ALLOC_LP_NUM_COLS, 
        ALLOC_LP_NUM_ROWS);

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        set_col_params(lp, i, RATIONAL, minyield, 1.0, 1.0);
    }

    // initialize matrix row and element
    elt = 0;
    row = 0;

    // (A) for all J,R: y_{j}*b_{j,r} <= u_{n,r} - a_{j,r}
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        for (k = 0; k < flex_soln->prob->num_resources; k++) {
            set_row_params(lp, row, LP_IGNORE, flex_soln->prob->nodes[
                flex_soln->mapping[i]]->unit_capacities[k] -
                flex_soln->prob->jobs[i]->unit_rigid_requirements[k] - 
                EPSILON);
            ia[elt] = row; ja[elt] = ALLOC_LP_COL_Y_I; 
            ra[elt] = flex_soln->prob->jobs[i]->unit_fluid_needs[k];
            elt++; row++;
        }
    }

    // last constraint is kind of funny...
    for (j = 0; j < flex_soln->prob->num_nodes; j++) {
        for (k = 0; k < flex_soln->prob->num_resources; k++) {
            set_row_params(lp, ALLOC_LP_ROW_B_JK, LP_IGNORE, 
                compute_available_resource(flex_soln, j, k) - EPSILON);
        }
    }

    // (B) for all N,R: sum_J y_{j}*b_{j,r} <= t_{n,r} - sum_J c_{j,r}
    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        j = flex_soln->mapping[i];
        for (k = 0; k < flex_soln->prob->num_resources; k++) {
            ia[elt] = ALLOC_LP_ROW_B_JK; ja[elt] = ALLOC_LP_COL_Y_I; 
            ra[elt] = flex_soln->prob->jobs[i]->total_fluid_needs[k];
            elt++;
        }
    }

    // now load matrix...
    load_matrix(lp, ALLOC_LP_NUM_ELTS, ia, ja, ra);

    return lp;
}

void maximize_average_yield_given_minimum(
    flexsched_solution_t flex_soln, double minyield)
{
    int i;
    linear_program_t lp = create_allocation_lp(flex_soln, minyield);

    if (solve_linear_program(lp, RATIONAL)) return;

    for (i = 0; i < flex_soln->prob->num_jobs; i++) {
        flex_soln->yields[i] = MAX(EPSILON, get_col_val(lp, ALLOC_LP_COL_Y_I) - EPSILON);
    }
    free_linear_program(lp);

    return;
}
