#include "flexsched.h"
#include "linearprog.h"

#ifdef CPLEX

// for now lps always maximize...
linear_program_t new_linear_program(int num_cols, int num_rows)
{
    int status = 0;
    linear_program_t lp = malloc(sizeof(struct linear_program_s));
    lp->env = CPXopenCPLEX(&status);
    lp->lp = CPXcreateprob(lp->env, &status, "");
    CPXnewcols(lp->env, lp->lp, num_cols, NULL, NULL, NULL, NULL, NULL);
    CPXnewrows(lp->env, lp->lp, num_rows, NULL, NULL, NULL, NULL);
    CPXchgobjsen(lp->env, lp->lp, CPX_MAX);
    return lp;
}

void free_linear_program(linear_program_t lp)
{
    CPXfreeprob(lp->env, &(lp->lp));
    CPXcloseCPLEX(&(lp->env));
    return;
}

void set_col_params(linear_program_t lp, int col, 
    int rational, double lb, double ub, double obj)
{
    char type;
    type = 'L';
    CPXchgbds(lp->env, lp->lp, 1, &col, &type, &lb);
    type = 'U';
    CPXchgbds(lp->env, lp->lp, 1, &col, &type, &ub);
    type = rational ? (CPX_CONTINUOUS) : (CPX_BINARY);
    CPXchgctype(lp->env, lp->lp, 1, &col, &type);
    CPXchgobj(lp->env, lp->lp, 1, &col, &obj);
    return;
}

void set_row_params(
    linear_program_t lp, int row, double lb, double ub)
{
    char sense;
    double bound, range;
    if (LP_IGNORE == lb) { 
        sense = 'L';
        bound = ub;
        range = 0.0;
    } else if (LP_IGNORE == ub) {
        sense = 'G';
        bound = lb;
        range = 0.0;
    } else if (lb == ub) {
        sense = 'E';
        bound = lb;
        range = 0.0;
    } else {
        sense = 'R';
        bound = lb;
        range = ub - lb;
    }
    CPXchgsense(lp->env, lp->lp, 1, &row, &sense);
    CPXchgrhs(lp->env, lp->lp, 1, &row, &bound);
    CPXchgrngval(lp->env, lp->lp, 1, &row, &range);
    return;
}

void load_matrix(linear_program_t lp, int elts, int *ia, int *ja, double *ra)
{
    CPXchgcoeflist(lp->env, lp->lp, elts, ia, ja, ra);
    return;
}

int solve_linear_program(linear_program_t lp, int rational)
{
    int status;
    double val;

    if (rational) {
        CPXchgprobtype(lp->env, lp->lp, CPXPROB_FIXEDMILP);
        status = CPXlpopt(lp->env, lp->lp);
    } else {
        CPXchgprobtype(lp->env, lp->lp, CPXPROB_MILP);
        CPXsetintparam(lp->env, CPX_PARAM_TILIM, LP_MILP_MAX_SECONDS);
        CPXmipopt(lp->env, lp->lp);
        status = CPXgetobjval(lp->env, lp->lp, &val);
    }

    return status;
}

double get_obj_val(linear_program_t lp) {
    double retval;
    CPXgetobjval(lp->env, lp->lp, &retval);
    return retval;
}

int get_mip_col_val(linear_program_t lp, int col) {
    double val;
    CPXgetx(lp->env, lp->lp, &val, col, col);
    return (int)(val + EPSILON);
}

double get_col_val(linear_program_t lp, int col) {
    double val;
    CPXgetx(lp->env, lp->lp, &val, col, col);
    return val;
}

#else

linear_program_t new_linear_program(int num_cols, int num_rows)
{
    linear_program_t lp = glp_create_prob();
    glp_add_cols(lp, num_cols);
    glp_add_rows(lp, num_rows);
    glp_set_obj_dir(lp, GLP_MAX);
    return lp;
}

void free_linear_program(linear_program_t lp)
{
    glp_delete_prob(lp);
    return;
}

void set_col_params(linear_program_t lp, int col, 
    int rational, double lb, double ub, double obj)
{
    glp_set_col_kind(lp, col+1, rational ? (GLP_CV) : (GLP_BV));
    glp_set_col_bnds(lp, col+1, GLP_DB, lb, ub);
    glp_set_obj_coef(lp, col+1, obj);
}

void set_row_params(
    linear_program_t lp, int row, double lb, double ub)
{
    int type = GLP_FR;
    if (LP_IGNORE == lb) { 
        if (LP_IGNORE == ub) {
            type = GLP_FR;
        } else {
            type = GLP_UP;
        }
    } else if (LP_IGNORE == ub) {
        type = GLP_LO;
    } else if (lb == ub) {
        type = GLP_FX;
    } else {
        type = GLP_DB;
    }
    glp_set_row_bnds(lp, row+1, type, lb, ub);
    return;
}

void load_matrix(linear_program_t lp, int elts, int *ia, int *ja, double *ra)
{
    int i;
    int *ib, *jb;
    double *rb;
    ib = calloc(elts+1, sizeof(int));
    jb = calloc(elts+1, sizeof(int));
    rb = calloc(elts+1, sizeof(double));

    for (i = 0; i < elts; i++) {
        ib[i+1] = ia[i]+1;
        jb[i+1] = ja[i]+1;
        rb[i+1] = ra[i];
    }
    glp_load_matrix(lp, elts, ib, jb, rb);
    return;
}

static int glp_stderr_null_out(void *info, const char *s) 
{
    return 1;
}

int solve_linear_program(linear_program_t lp, int rational)
{
    int solver_status, solution_status;

    // kill lp output
    glp_term_hook(&glp_stderr_null_out, NULL);

    if (rational) {
        solver_status = glp_simplex(lp, NULL);
        solution_status = (glp_get_status(lp) != GLP_OPT);
    } else {
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.msg_lev = GLP_MSG_ERR;
        parm.tm_lim = GLPK_TIME_LIMIT;
        parm.presolve = GLP_ON;
        parm.binarize = GLP_ON;
        glp_intopt(lp, &parm);
        solution_status = (glp_mip_status(lp) != GLP_OPT);
    }

    // Reached the time limit
    if (solver_status == GLP_ETMLIM) {
        return 2;
    }

    return (solver_status || solution_status);
}

double get_obj_val(linear_program_t lp) {
    if (GLP_OPT == glp_mip_status(lp)) return glp_mip_obj_val(lp);
    return glp_get_obj_val(lp);
}

int get_mip_col_val(linear_program_t lp, int col) {
    return (int)(glp_mip_col_val(lp, col+1) + EPSILON);
}

double get_col_val(linear_program_t lp, int col) {
    if (GLP_OPT == glp_mip_status(lp)) return glp_mip_col_val(lp, col+1);
    return glp_get_col_prim(lp, col+1);
}

#endif
