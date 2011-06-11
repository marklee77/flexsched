#define RATIONAL 1
#define INTEGER 0

#define LP_IGNORE 666

#ifdef CPLEX

#include <ilcplex/cplex.h>

typedef struct linear_program_s {
    CPXENVptr *env;
    CPXLPptr *lp;
} *linear_program_t;

#else

#define GLPK_TIME_LIMIT 10*60*1000

#include <glpk.h>

typedef glp_prob *linear_program_t;

#endif

linear_program_t new_linear_program(int num_cols, int num_rows);
void free_linear_program(linear_program_t lp);
void set_col_params(linear_program_t lp, int col, 
    int rational, double lb, double ub, double obj);
void set_row_params(linear_program_t lp, int row, double lb, double ub);
void load_matrix(linear_program_t lp, int elts, int ia[], int ja[], float ra[]);
int solve_linear_program(linear_program_t lp, int rational);
double get_obj_val(linear_program_t lp);
int get_mip_col_val(linear_program_t lp, int col);
double get_col_val(linear_program_t lp, int col);
