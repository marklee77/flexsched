#include "flexsched.h"

/* Global for some VP stuff acceleration*/ 
float **global_vp_bin_loads; 

/* generate_vp_instance() */
struct vp_instance *generate_vp_instance(float yield) 
{
  int i, j, dim_counter;
  struct vp_instance *instance;

  // INPUT
  instance = (struct vp_instance *)calloc(1, sizeof(struct vp_instance));
  instance->num_dims = INS.numrigid + INS.numfluid;
  instance->num_vectors = INS.numservices;
  instance->num_bins = INS.numservers;
  instance->vectors = (struct vp_vector *)
    calloc(instance->num_vectors, sizeof(struct vp_vector));
  instance->bins = (struct vp_vector *)
    calloc(instance->num_bins, sizeof(struct vp_vector));
  instance->mapping = (int *)calloc(instance->num_vectors,sizeof(int));
  for (i=0; i<instance->num_vectors; i++) {
    instance->vectors[i].service = i;
    instance->vectors[i].num_dims = INS.numrigid + INS.numfluid;
    instance->vectors[i].x = (float *)calloc(instance->num_dims, sizeof(float));
    dim_counter = 0;
    for (j=0; j < INS.numrigid; dim_counter++, j++) {
      instance->vectors[i].x[dim_counter] = INS.rigidneeds[i][j]; 
    }
    for (j=0; j < INS.numfluid; dim_counter++, j++) {
      instance->vectors[i].x[dim_counter] = 
              ((1.0 - INS.slas[i])*yield + INS.slas[i]) * INS.fluidneeds[i][j]; 
    }
  }
  for (i = 0; i < instance->num_bins; i++) {
    instance->bins[i].num_dims = INS.numrigid + INS.numfluid;
    instance->bins[i].x = (float *)calloc(instance->num_dims, sizeof(float));
    dim_counter = 0;
    for (j = 0; j < INS.numrigid; dim_counter++, j++) {
        instance->bins[i].x[dim_counter] = INS.rigidcapacities[i][j];
    }
    for (j = 0; j < INS.numrigid; dim_counter++, j++) {
        instance->bins[i].x[dim_counter] = INS.fluidcapacities[i][j];
    }
  }

  // OUTPUT
  for (i=0; i < instance->num_vectors; i++)
    instance->mapping[i] = -1; // initialized to "not mapped"

#if 0
  fprintf(stderr,"VP instance:\n");
  for (i=0; i<instance->num_vectors; i++) {
    int j;
    for(j=0; j<instance->num_dims; j++) 
      fprintf(stderr,"%.2f ",instance->vectors[i].x[j]);
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"----------------\n");
#endif

  return instance;
}

/* free_vp_instance */
void free_vp_instance(struct vp_instance *instance)
{
  int i;
  for (i=0; i<instance->num_vectors; i++) {
    free(instance->vectors[i].x);
  }
  for (i=0; i<instance->num_bins; i++) {
    free(instance->bins[i].x);
  }
  free(instance->vectors);
  free(instance->bins);
  free(instance->mapping);
  free(instance);
  return;
}

/* Lexicographical comparison of vectors */
int VP_sort_LEX(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;

  for (i=0; i < vx->num_dims; i++) {
    if (vx->x[i] < vy->x[i])
      return 1;
    if (vx->x[i] > vy->x[i])
      return -1;
  }
  return 0;
}

x comparison of vectors */
int VP_sort_MAX(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
  float max_x, max_y;

  max_x = vx->x[0];
  max_y = vy->x[0];

  for (i=0; i < vx->num_dims; i++) {
    if (vx->x[i] > max_x)
      max_x = vx->x[i];
    if (vy->x[i] > max_y)
      max_y = vy->x[i];
  }
  if (max_x > max_y)
    return -1;
  if (max_x < max_y)
    return 1;
  return 0;
}

/* Sum comparison of vectors */
int VP_sort_SUM(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;
  float sum_x, sum_y;

  sum_x = 0.0;
  sum_y = 0.0;

  for (i=0; i < vx->num_dims; i++) {
    sum_x += vx->x[i];
    sum_y += vy->x[i];
  }

  if (sum_x > sum_y)
    return -1;
  if (sum_x < sum_y)
    return 1;
  return 0;
}

/* Misc comparison of vectors  */
/* By decreasing order of MISC */
int VP_sort_MISC(const void *x, const void *y)
{
  int i;
  struct vp_vector *vx = (struct vp_vector *)x;
  struct vp_vector *vy = (struct vp_vector *)y;

  if (vx->misc > vy->misc)
    return -1;
  if (vx->misc < vy->misc)
    return 1;
  return 0;
}

int vp_vector_can_fit_in_bin(struct vp_instance *vp, int v, int b)
{
  int i,j;
  float load[vp->num_dims];

  for (j=0; j < vp->num_dims; j++)
    load[j] = 0.0;

  // Compute the load
  for (i=0; i < vp->num_vectors; i++) {
    if (vp->mapping[i] != b)
      continue;
    for (j=0; j < vp->num_dims; j++) {
      load[j] += vp->vectors[i].x[j];
    }
  }

  // Check that the load can accomodate the new vector
  for (j=0; j < vp->num_dims; j++) {
    if (load[j] + vp->vectors[v].x[j] > vp->bins[b].x[j])
      return 0;
  }
  return 1;
}

// load computed as the sum of all dims of all objects
int vp_vector_can_fit_in_bin_fast(struct vp_instance *vp, int v, int b)
{
  int i,j;

  for (j=0; j < vp->num_dims; j++)  {
    if (global_vp_bin_loads[b][j] + vp->vectors[v].x[j] > vp->bins[b].x[j]) {
      return 0;
    }
  }
  return 1;
}

// load computed as the sum of all dims of all objects
float vp_compute_bin_load(struct vp_instance *vp, int b)
{
  float load = 0.0;
  int i,j;

  for (i=0; i < vp->num_vectors; i++) {
    if (vp->mapping[i] != b)
      continue;
    for (j=0; j < vp->num_dims; j++)
      load += vp->vectors[i].x[j];
  }
  return load;
}

// load computed as the sum of all dims of all objects
float vp_compute_bin_load_fast(struct vp_instance *vp, int b)
{
  float load = 0.0;
  int j;

  for (j=0; j < vp->num_dims; j++) {
    load += global_vp_bin_loads[b][j];
  }
  return load;
}

/* A helper function to compute the thingies from
 * the Maruyama article, and puts them into the misc 
 * field of the vp->vectors
 */
void vp_compute_degrees_of_dominance(struct vp_instance *vp)
{
  int k,j,i, status;
  float *sums, *degrees;

  degrees = (float *)calloc(vp->num_dims, sizeof(float));
  sums = (float *)calloc(vp->num_bins, sizeof(float));

  // For each dimension compute its degree of dominance
  for (j=0; j < vp->num_dims; j++) {
    // sets the sums to zero
    for (k=0; k < vp->num_bins; k++)
      sums[k] = 0.0;
    // compute the sums
    for (i=0; i < vp->num_vectors; i++) {
      if (vp->mapping[i] < 0) {
        free(degrees);
        free(sums);
        fprintf(stderr,"Warning: can't compute degrees of dominance due to unmapped vector!\n");
        return;
      }
      sums[vp->mapping[i]] += vp->vectors[i].x[j];
    }
    // compute the ceil
    for (k=0; k < vp->num_bins; k++)
      sums[k] = ceilf(sums[k] - 1.0);
    // compute the degree
    degrees[j] = 0.0;
    for (k=0; k < vp->num_bins; k++)
      degrees[j] += sums[k];
    degrees[j] /= vp->num_bins;
  }

  for (i=0; i < vp->num_vectors; i++) {
    vp->vectors[i].misc = 0.0;
    for (j=0; j < vp->num_dims; j++)
      vp->vectors[i].misc += degrees[j] * vp->vectors[i].x[j];
  }

  free(sums);
  free(degrees);
  return;
}

/* Implementation of The standard "Fit" algorithms for
 * solving binpacking:
 *    "FIRST"   or "BEST":  fitting policy
 *    "LEX", "MAX", or "SUM": sorting policy
 *    "MARUYAMA" sorting: the Maruyama optimization
 *
 * Returns the number of used bins
 */
//#define FIT_DEBUG 0
int solve_vp_instance_FITD(const char *fit_type,
                           const char *sort_type,
                           struct vp_instance *vp)
{
  int i,j,k;
  int (*compar)(const void *, const void *);
  int max_bin;

  // Sort the vectors in the instance according to the sort type
  if (!strcmp(sort_type,"LEX")) {
    compar = VP_sort_LEX;
  } else if (!strcmp(sort_type,"MAX")) {
    compar = VP_sort_MAX;
  } else if (!strcmp(sort_type,"SUM")) {
    compar = VP_sort_SUM;
  } else if (!strcmp(sort_type,"MARUYAMA")) {
    compar = VP_sort_MISC;
  } else {
    fprintf(stderr,"Invalid VP sort type '%s'\n",sort_type);
    exit(1);
  }
  qsort(vp->vectors, vp->num_vectors,
        sizeof(struct vp_vector), compar);

 // Allocate and initialize vp_bin_loads 
  global_vp_bin_loads = (float **)calloc(vp->num_bins,sizeof(float*));
  for (k=0; k < vp->num_bins; k++) {
    global_vp_bin_loads[k] = (float *)calloc(vp->num_dims, sizeof(float));
  }

  max_bin = 0;

  // Place vectors into bins
  for (i=0; i < vp->num_vectors; i++) {
    if (!strcmp(fit_type, "FIRST")) { // First Fit
      for (j=0; j < vp->num_bins; j++) {
        if (vp_vector_can_fit_in_bin_fast(vp, i, j)) {
          vp->mapping[i] = j;
          for (k=0; k < vp->num_dims; k++) {
            global_vp_bin_loads[j][k] += vp->vectors[i].x[k];
          }
      if (j > max_bin)
            max_bin = j;
          break;
        }
      }
    } else if (!strcmp(fit_type, "BEST")) { // Best Fit
      float load, max_load = -1.0;
#ifdef FIT_DEBUG
      fprintf(stderr,"Trying to fit vector ");
      for (j=0; j<vp->num_dims; j++)
        fprintf(stderr,"%.2f ",vp->vectors[i].x[j]);
      fprintf(stderr,"(svc=%d) into a bin\n",i);
#endif
      for (j=0; j<vp->num_bins; j++) {
        if (!vp_vector_can_fit_in_bin_fast(vp, i, j)) {
#ifdef FIT_DEBUG
          fprintf(stderr,"  can't fit in bin %d\n",j);
#endif
          continue;
        }
        load = vp_compute_bin_load_fast(vp, j);
#ifdef FIT_DEBUG
        fprintf(stderr,"  bin %d has load %.3f\n",j,load);
#endif
        if ((max_load == -1) || (load > max_load)) {
          max_load = load;
          vp->mapping[i] = j;
          for (k=0; k < vp->num_dims; k++) {
            global_vp_bin_loads[j][k] += vp->vectors[i].x[k];
          }
      if (j > max_bin)
            max_bin = j;
        }
      }
    } else {
      fprintf(stderr,"Invalid VP fit type '%s'\n",fit_type);
      exit(1);
    }
#ifdef FIT_DEBUG
    fprintf(stderr,"vector %d (%.3f %.3f) was put in bin %d\n",
           i,
           vp->vectors[i].x[0],
           vp->vectors[i].x[1],
           vp->mapping[i]);
#endif
  }
  //vp->num_bins = (max_bin+1);

  for (k=0; k < vp->num_vectors; k++) {
    free(global_vp_bin_loads[k]);
  }
  free(global_vp_bin_loads);
 // Apply the Marumaya optimization in case it helps, if
  // we didn't just do it
  // FIXME: num_bins should never be bigger than num_servers these days
  if (vp->num_bins > INS.numservers) {
    if (strcmp(sort_type,"MARUYAMA")) {
      int status;
      int old_num_bins = vp->num_bins;
      vp_compute_degrees_of_dominance(vp);
      status = solve_vp_instance_FITD(fit_type,"MARUYAMA",vp);
      if (vp->num_bins <= INS.numservers) {
        strcat(INS.misc_output,"M");
        // fprintf(stderr,"MARUYAMA HELPED! (%d -> %d)\n", old_num_bins, vp->num_bins);
      }
      return status;
    }
  }

  return VP_SOLVED;
}

void find_top_two_dims(float *array, int n,
                       int *d1, int *d2)
{
  int i;
  float max2;

  *d1 = array_argmax(array,n);
  max2 = -1.0;
  for (i=0; i < n; i++) {
    if (i == *d1)
      continue;
    if ((max2 == -1.0) || (array[i] > max2)) {
      max2 = array[i];
      *d2 = i;
    }
  }

  return;
}

/* MCB_compar_MAX */
int MCB_compar_MAX(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);

  if (max1 > max2) {
    return -1;
  } else if (max2 > max1) {
    return 1;
  } else {
    return 0;
  }
}

/* MCB_compar_SUM */
int MCB_compar_SUM(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float sum1, sum2;

  sum1 = array_sum(v1->x, v1->num_dims);
  sum2 = array_sum(v2->x, v2->num_dims);

  if (sum1 > sum2)
    return -1;
  else if (sum2 > sum1)
    return 1;
  else
    return 0;
}

/* MCB_compar_MAXRATIO */
int MCB_compar_MAXRATIO(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;
  float min1, min2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);
  min1 = array_min(v1->x, v1->num_dims);
  min2 = array_min(v2->x, v2->num_dims);

  if (max1/min1 > max2/min2)
    return -1;
  else if (max2/min2 > max1/min1)
    return 1;
  else
    return 0;
}

/* MCB_compar_MAXDIFF */
int MCB_compar_MAXDIFF(const void *p1, const void *p2)
{
  struct vp_vector *v1 = *(struct vp_vector **)p1;
  struct vp_vector *v2 = *(struct vp_vector **)p2;
  float max1, max2;
  float min1, min2;

  max1 = array_max(v1->x, v1->num_dims);
  max2 = array_max(v2->x, v2->num_dims);
  min1 = array_min(v1->x, v1->num_dims);
  min2 = array_min(v2->x, v2->num_dims);

  if (max1-min1 > max2-min2)
    return -1;
  else if (max2-min2 > max1-min1)
    return 1;
  else
    return 0;
}

int int_float_compar(const void *x, const void *y)
{
  struct int_float ifx = *(struct int_float *)x;
  struct int_float ify = *(struct int_float *)y;

  if (ifx.f > ify.f)
    return 1;
  else if (ifx.f < ify.f)
    return -1;
  else
    return 0;
}

/* helper function 
 *   Takes a bin_load and returns an array of dimension indices
 *   sorted from the lowest loaded dim to the highest loaded dim 
 */
void sort_dims_by_increasing_load(int num_dims, int *sorted_dims, float *bin_load)
{
  struct int_float *sorted;
  int i;

  sorted = (struct int_float*)calloc(num_dims,sizeof(struct int_float));
  for (i=0; i<num_dims; i++) {
    sorted[i].i = i;
    sorted[i].f = bin_load[i];
  }

  qsort(sorted,num_dims,sizeof(struct int_float),int_float_compar);
  for (i=0; i<num_dims; i++)
    sorted_dims[i] = sorted[i].i;
  free(sorted);
}

void add_to_list_of_lists(int first_dim, int second_dim,
                         int *num_lists, int **first_dims, int **second_dims)
{
  int i;

  // If already in the lists, foret it
  for (i=0; i < *num_lists; i++) {
    if ((first_dim == (*first_dims)[i]) &&
        (second_dim == (*second_dims)[i])) {
        return;
    }
  }

  // Add to the list 
  (*num_lists)++;
  *first_dims = (int *)REALLOC(*first_dims, (*num_lists) * sizeof(int));
  *second_dims = (int *)REALLOC(*second_dims, (*num_lists) * sizeof(int));
  (*first_dims)[(*num_lists)-1] = first_dim;
  (*second_dims)[(*num_lists)-1] = second_dim;

  return;
}

/* MCB_find_and_extract_feasible_vector_in_list() 
 *  Given a vector list, and a bin load, find the first vector
 *  that fits, and remove it from the list (return it)
 *
 *  The extraction just erases the end of the list (with shifting)
 *  just to be a bit more efficient
 */
struct vp_vector *MCB_find_and_extract_feasible_vector_in_list(
                  int num_dims, float *bin_load, float *bin_capacity,
                  struct vp_vector **list, int *list_size)
{
  struct vp_vector *vector;
  int i,j;

  for (i=0; i < *list_size; i++) {
    // fprintf(stderr,"  Looking at: ");
    // for (j=0; j < num_dims; j++) 
    //   fprintf(stderr,"%.2f ",list[i]->x[j]);
    // fprintf(stderr,"\n");
    for (j=0; j < num_dims; j++) {
      // FIXME: should probably be using vector fits in bin function here...
      if (bin_capacity[j] - bin_load[j] < list[i]->x[j])
        break;
    }
    if (j >= num_dims)
      break;
  }

  if (i >= *list_size)
    return NULL;

  // extract it
  vector = list[i];
  for (j=i; j < *list_size-1; j++)
    list[j] = list[j+1];
  (*list_size)--;

  return vector;
}

/* MCB_pick_vector(): The core of the BCB implementation 
 *   right now implements CP only
 */
struct vp_vector *MCB_pick_vector(int isPP,
                                  int num_dims, float *bin_load,
                                  float *bin_capacity,
                                  struct vp_vector ****lists,
                                  int **list_sizes)
{
  int i,j;
  int *sorted_dims;
  int bottom_d, top_d;

  //fprintf(stderr,"IN MCB_PICK_VECTOR\n");

  // sort the bin dims by increasing load
  sorted_dims = (int *)calloc(num_dims,sizeof(int));
  sort_dims_by_increasing_load(num_dims, sorted_dims, bin_load);

  // Go through the grinds...
  //   **.....  (best case)
  //   *.*....  (d1 ok, and other)
  //   .**....  (d2 ok, and other)
  //   *..*...  etc..
  //   .*.*...
  //   *...*..
  //   .*..*..

  // Construct an array of list indices to look at, in the
  // right order 
  int num_lists = 0;
  int *first_dims=NULL, *second_dims=NULL;

 for (bottom_d = 1; bottom_d < num_dims; bottom_d++) {
    int first_dim, second_dim;
    for (top_d = 0; top_d < bottom_d; top_d++) {
      // Top dim and d
      first_dim  = sorted_dims[bottom_d];
      second_dim = sorted_dims[top_d];
      if (isPP) {
        if (first_dim < second_dim) {
         add_to_list_of_lists(first_dim, second_dim,
                         &num_lists, &first_dims, &second_dims);
         add_to_list_of_lists(second_dim, first_dim,
                         &num_lists, &first_dims, &second_dims);
        } else {
         add_to_list_of_lists(second_dim, first_dim,
                         &num_lists, &first_dims, &second_dims);
         add_to_list_of_lists(first_dim, second_dim,
                         &num_lists, &first_dims, &second_dims);
        }
      } else {
        first_dim  = MAX(sorted_dims[bottom_d],sorted_dims[top_d]);
        second_dim = MIN(sorted_dims[bottom_d],sorted_dims[top_d]);
        add_to_list_of_lists(first_dim, second_dim,
                         &num_lists, &first_dims, &second_dims);
      }
    }
  }

  free(sorted_dims);

#ifdef MCB_DEBUG
  fprintf(stderr,"Lists order: ");
  for (i=0; i< num_lists; i++)
    fprintf(stderr,"[%d][%d],",first_dims[i],second_dims[i]);
  fprintf(stderr,"\n");
#endif
 // Go through the lists in order and look for a feasible vector to pick
  struct vp_vector *vector = NULL;
  for (i=0; i < num_lists; i++) {
    struct vp_vector **list = lists[first_dims[i]][second_dims[i]];
    int *list_size = &(list_sizes[first_dims[i]][second_dims[i]]);
#ifdef MCB_DEBUG
    fprintf(stderr,"  Looking for a feasible vector in LIST[%d][%d] (size=%d):",
                    first_dims[i],second_dims[i],*list_size);
    for (j=0; j < *list_size; j++)
      fprintf(stderr,"%d ",list[j]->service);
    fprintf(stderr,"\n");
#endif
    vector = MCB_find_and_extract_feasible_vector_in_list(
                  num_dims, bin_load, bin_capacity, list, list_size);
    if (vector) {// We found a vector in the list
//      fprintf(stderr,"    Picked vector corr to svc. %d in LIST[%d][%d]\n",
//                vector->service, first_dims[i],second_dims[i]);
      // fprintf(stderr,"  Found vector: ");
      // int j;
      // for (j=0; j < vector->num_dims; j++)
      //   fprintf(stderr,"%.2f ",vector->x[j]);
      // fprintf(stderr,"\n");
      break;
    }
  }

  free(first_dims); free(second_dims);

  // If we found a vector, then great
  if (vector)
    return vector;
  else
    return NULL;
}
/*
 * The MCB implementation
 *   "PP" (permutation pack) or "CP" (choose pack)
 *   sorting type: "MAX", "SUM", "MAXDIFF", "MAXRATIO"
 * 
 *  Everything implemented for w=2
 */
int solve_vp_instance_MCB(char *pack, char *sort_type,
                          struct vp_instance *vp)
{
  int num_lists,i,j,d1,d2;
  int (*compar)(const void *, const void *);
  char isPP,isCP;

  if (!strcmp(pack,"PP")) {
    isPP = 1; isCP = 0;
  } else if (!strcmp(pack,"CP")) {
    isPP = 0; isCP = 1;
  } else {
    fprintf(stderr,"MCB: unknown pack type '%s'\n",pack);
    exit(1);
  }

#ifdef MCB_DEBUG
  fprintf(stderr,"#############################################\n");
  fprintf(stderr,"########## SOLVING AN INSTANCE ##############\n");
  fprintf(stderr,"#############################################\n");
#endif

  // Initialize lists
  struct vp_vector ****lists;
  int **list_sizes;

  lists = (struct vp_vector ****)calloc(vp->num_dims,
             sizeof(struct vp_vector ***));

  list_sizes = (int **)calloc(vp->num_dims,
             sizeof(int *));
  for (d1=0; d1 < vp->num_dims; d1++) {
    lists[d1] = (struct vp_vector ***)calloc(vp->num_dims,
             sizeof(struct vp_vector **));
    list_sizes[d1] = (int *)calloc(vp->num_dims,
             sizeof(int));
  }

  // For CP we need only to initialize the triangular
  // matrix but what the heck.. this doesn't hurt
  for (d1=0; d1 < vp->num_dims; d1++) {
    for (d2=0; d2 < vp->num_dims; d2++) {
      lists[d1][d2] = NULL;
      list_sizes[d1][d2] = 0;
    }
  }

  // Put vectors in lists
  for (i=0; i < vp->num_vectors; i++) {
    int d1, d2;

    // Find the indices of the two largest coordinates of the vector
    // with d1 being the index of the dimension with the 
    // largest value and d2 the index of the dimension with
    // second largest value
    find_top_two_dims(vp->vectors[i].x,
                      vp->vectors[i].num_dims,&d1,&d2);

    // If CP, then sort the dims indices, since we only
    // need a triangular matrix
    if (isCP) {
      if (d1 < d2) {
        int tmp = d1;
        d1 = d2;
        d2 = tmp;
      }
    }

    // Put the vector in the list
#ifdef MCB_DEBUG
    fprintf(stderr,"Putting vector ");
    for (j=0; j < vp->num_dims; j++)
      fprintf(stderr,"%.2f ",vp->vectors[i].x[j]);
    fprintf(stderr," (svc %d) in LIST[%d][%d]\n",
                  vp->vectors[i].service,
                  d1,d2);
#endif
    list_sizes[d1][d2]++;
    lists[d1][d2] = (struct vp_vector **)REALLOC(lists[d1][d2],
                         list_sizes[d1][d2] * sizeof (struct vp_vector *));
    lists[d1][d2][list_sizes[d1][d2] - 1] = &(vp->vectors[i]);
  }

  // Sort the lists according to the criteria
  // Select the sorting function
  if (!strcmp(sort_type,"MAX")) {
    compar = MCB_compar_MAX;
  } else if (!strcmp(sort_type,"SUM")) {
    compar = MCB_compar_SUM;
  } else if (!strcmp(sort_type,"MAXDIFF")) {
    compar = MCB_compar_MAXDIFF;
  } else if (!strcmp(sort_type,"MAXRATIO")) {
    compar = MCB_compar_MAXRATIO;
  } else {
    fprintf(stderr,"Invalid MCB compare type '%s'\n",sort_type);
    exit(1);
  }

  // Sort all the lists
  for (d1=0; d1 < vp->num_dims; d1++) {
    int d2_loop_bound;
    if (isCP)
      d2_loop_bound = d1;  // only triangular
    else
      d2_loop_bound = vp->num_dims; // full matrix
    for (d2=0; d2 < d2_loop_bound; d2++) {
      qsort(lists[d1][d2],list_sizes[d1][d2], sizeof(struct vp_vector *), compar);
#ifdef MCB_DEBUG
      fprintf(stderr,"  Sorted LIST[%d][%d] (size=%d):",
                      d1,d2,list_sizes[d1][d2]);
      for (j=0; j < list_sizes[d1][d2]; j++)
        fprintf(stderr,"%d ",lists[d1][d2][j]->service);
      fprintf(stderr,"\n");
#endif
    }
  }

  // Go through bins
  int bin = 0;
  int num_vectors_picked = 0;

  float load[vp->num_bins][vp->num_dims];

  // Set the loads to zero
  for (i=0; i < vp->num_bins; i++)
    for (j=0; j < vp->num_dims; j++)
      load[i][j] = 0.0;

  while (bin < vp->num_bins && num_vectors_picked < vp->num_vectors) {
    struct vp_vector *vector;
    // Pick the vector to put in the bin
    // This call removes the vector from lists
#ifdef MCB_DEBUG
    fprintf(stderr,"** BIN %d load= ",bin);
    for (j=0; j < vp->num_dims; j++)
      fprintf(stderr,"%.2f ",load[bin][j]);
    fprintf(stderr,"\n");
#endif

    vector = MCB_pick_vector(isPP, vp->num_dims, load[bin], vp->bins[bin].x,  lists, list_sizes);
    if (!vector) { // Couldn't find a vector, go to next bin
      bin++;
      continue;
    }

    // Update the mapping
    vp->mapping[vector->service] = bin;
    // Update the load of the bin
    for (j=0; j < vp->num_dims; j++)
      load[bin][j] += vector->x[j];
    // Are we done?
    num_vectors_picked++;
  }

  // Free some memory
  for (d1=0; d1 < vp->num_dims; d1++) {
    free(list_sizes[d1]);
    for (d2=0; d2 < vp->num_dims; d2++)  {
      if (lists[d1][d2])
        free(lists[d1][d2]);
    }
    free(lists[d1]);
  }
  free(lists);
  free(list_sizes);

  return (num_vectors_picked < vp->num_vectors) ? VP_NOT_SOLVED: VP_SOLVED;
}

/* solve_vp_instance() 
 *  Returns the number of bins used
 */
int solve_vp_instance(struct vp_instance *vp,
                      char *vp_algorithm)
{
  int status;

  if (!strcmp(vp_algorithm,"FFDLEX")) {
    status = solve_vp_instance_FITD("FIRST","LEX",vp);
  } else if (!strcmp(vp_algorithm,"FFDMAX")) {
    status = solve_vp_instance_FITD("FIRST","MAX",vp);
  } else if (!strcmp(vp_algorithm,"FFDSUM")) {
    status = solve_vp_instance_FITD("FIRST","SUM",vp);
  } else if (!strcmp(vp_algorithm,"BFDLEX")) {
    status = solve_vp_instance_FITD("BEST","LEX",vp);
  } else if (!strcmp(vp_algorithm,"BFDMAX")) {
    status = solve_vp_instance_FITD("BEST","MAX",vp);
  } else if (!strcmp(vp_algorithm,"BFDSUM")) {
    status = solve_vp_instance_FITD("BEST","SUM",vp);
  } else if (!strcmp(vp_algorithm,"CPMAX")) {
    status = solve_vp_instance_MCB("CP","MAX",vp);
  } else if (!strcmp(vp_algorithm,"CPSUM")) {
    status = solve_vp_instance_MCB("CP","SUM",vp);
  } else if (!strcmp(vp_algorithm,"CPMAXRATIO")) {
    status = solve_vp_instance_MCB("CP","MAXRATIO",vp);
  } else if (!strcmp(vp_algorithm,"CPMAXDIFF")) {
    status = solve_vp_instance_MCB("CP","MAXDIFF",vp);
  } else if (!strcmp(vp_algorithm,"PPMAX")) {
    status = solve_vp_instance_MCB("PP","MAX",vp);
  } else if (!strcmp(vp_algorithm,"PPSUM")) {
    status = solve_vp_instance_MCB("PP","SUM",vp);
  } else if (!strcmp(vp_algorithm,"PPMAXRATIO")) {
    status = solve_vp_instance_MCB("PP","MAXRATIO",vp);
  } else if (!strcmp(vp_algorithm,"PPMAXDIFF")) {
    status = solve_vp_instance_MCB("PP","MAXDIFF",vp);
  } else if (!strcmp(vp_algorithm,"CHEKURI")) {
    status = solve_vp_instance_CHEKURI(vp);
  } else {
    fprintf(stderr,"Unknown vp_algorithm '%s'\n",vp_algorithm);
    exit(1);
  }
  return status;
}

/* 
 * VP_scheduler() 
 */
int VP_scheduler(char *vp_algorithm, char *ignore2, char *ignore3)
{
  double yield, yieldlb, yieldub, best_yield;
  struct vp_instance *vp = NULL;
  int i, status;
  int num_bins;

  yieldlb = 0.0;
  yieldub = compute_LP_bound();

  best_yield = -1.0;
  yield = 0.0;


  while((yield < 1.0 - EPSILON) && (yieldub - yieldlb > 0.001)) {
    yield = (yieldub + yieldlb) / 2.0;
#if 0
    fprintf(stderr,"Trying yield=%.6f (up=%.4f down=%.4f)\n",yield,yieldub,yieldlb);
#endif

    // Generate the VP instance
    if (vp)
      free_vp_instance(vp);
    vp = generate_vp_instance(yield);

    // Solve the VP instance
    if ((status = solve_vp_instance(vp, vp_algorithm)) == VP_SOLVED) {
      yieldlb = yield;
      best_yield = yield;
#if 0
      fprintf(stderr,"Instance was solved with %d bins\n",vp->num_bins);
#endif

    // Save the computed mapping
      for (i=0; i < INS.numservices; i++) {
        INS.mapping[vp->vectors[i].service] = vp->mapping[i];
        INS.allocation[vp->vectors[i].service] = best_yield;
#if 0
        fprintf(stderr,"SERVICE %d : on server %d with yield %.2f\n",
                vp->vectors[i].service, 
                mapping[vp->vectors[i].service], 
                allocation[vp->vectors[i].service]);
#endif
      }
    } else {
      yieldub = yield;
#if 0
      fprintf(stderr,"Instance was not solved (needed %d bins, and we have %d servers)\n",
              vp->num_bins,INS.numservers);
#endif
    }
  }
  if (vp)
    free_vp_instance(vp); // free the last instance,vp->num_bins

  // We never solved the VP instance
  if (best_yield <= 0)
    return RESOURCE_ALLOCATION_FAILURE;

  // We set the yield and mapping of all services

  return RESOURCE_ALLOCATION_SUCCESS;
}
