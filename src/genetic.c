#include "flexsched.h"

int GA_scheduler(char *arg1, char *arg2, char *arg3)
{
  int i;
  int initial = !strcmp(arg1,"F");
  int mutated = !strcmp(arg2,"F");
  int crossed = !strcmp(arg3,"F");

  // Map each service to a particular host
  if (GA_compute_mapping(initial,mutated,crossed) == RESOURCE_ALLOCATION_FAILURE) {
    return RESOURCE_ALLOCATION_FAILURE;
  }
#if 0
  fprintf(stderr,"GA MAPPING:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d on server %d\n",i,INS.mapping[i]);
  }
#endif

  // For each host, compute its services' allocations
  for (i=0; i<INS.numservers; i++) {
    compute_allocations_given_mapping(i);
  }
#if 0
  fprintf(stderr,"GA ALLOCATIONS:\n");
  for (i=0; i<INS.numservices; i++) {
    fprintf(stderr,"Service %d has an allocation (yield) of %f\n",i,INS.allocation[i]);
  }

#endif

  return RESOURCE_ALLOCATION_SUCCESS;
}

