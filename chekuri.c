#include "flexsched.h"

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
  //  (A) for all i        sum_j e_ij = 1:                n rows
  //  (B) for all k      sum_i p_ik x_ij <= 1                 d*m rows
  //fprintf(stderr,"adding %d rows\n",n+n*m);
  glp_add_rows(prob, n + n*m);

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
  // vp->num_bins = last_successful_num_bins;

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

