#include "flexsched.h"

#include <ga/ga.h>
#include <ga/GA1DArrayGenome.h>
#include <ga/std_stream.h>

#define cout STD_COUT
#define endl STD_ENDL
#define ofstream STD_OFSTREAM

#define OBJECTIVE_FEASIBLE_WEIGHT 2.0

/** C/C++ interface **/
extern "C" void initialize_global_server_loads();
extern "C" int GA_compute_mapping(int,int,int);
extern "C" int service_can_fit_on_server_fast(int, int);
extern "C" float compute_minimum_yield();
extern "C" void compute_allocations_given_mapping(int);

/* Prototypes */
float Objective(GAGenome&);
static void myInitializer(GAGenome&);
static int myMutator(GAGenome&, float);
static int myMutator2(GAGenome&, float);
static int makeGenomeFeasible(GAGenome&g);
static int myCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
static int pick_service_to_evict(GAGenome &g, int server, char type, int j); 
static int evict_service(GAGenome &g, int to_evict, int from_server);
static void print_genome(GAGenome &g);

/* Global variable */
static int makeInitialGenomesFeasible;
static int makeMutatedGenomesFeasible;
static int makeCrossedGenomesFeasible;

/* GA_compute_mapping() */
int GA_compute_mapping(int initial, int mutated, int crossed)
{
  int i,j;
  unsigned int seed=0;

  makeInitialGenomesFeasible = initial;
  makeMutatedGenomesFeasible = mutated;
  makeCrossedGenomesFeasible = crossed;

  /* Seed the Random Number Generator */
  for (i=0; i<INS.numservices; i++) {
    seed += (int)(100.0 * INS.fluidneeds[i][0]);
  }
  GARandomSeed(seed);
  
  /* Initialize the server load (for speed of computation) */
  initialize_global_server_loads();

  /* Create a genome: used to clone a population of genomes */
  /* Use my own objective */
  GA1DArrayGenome<int> genome(INS.numservices, Objective);

  /* Use my own initializer */
  genome.initializer(myInitializer);
  /* Use my own mutator */
  genome.mutator(myMutator);
  /* Use my own crossover */
  genome.crossover(myCrossover);

  /* Instantiate the genetic algorithm */
  GASimpleGA ga(genome);

  /* Set the GA parameters */
#if 1
  int popsize = 100;
  int ngen = 2000;
  float pmut = 0.1;
  float pcross = 0.25;
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
#endif

//  ga.parameters("./settings.txt");
  
  /* Solve */
  ga.evolve();

  /* Print out statistics */
//  cout << "The GA did " << ga.statistics().crossovers() << " crossovers\n";
//  cout << "The GA did " << ga.statistics().mutations() << " mutations\n";
//  cout << "The GA did " << ga.statistics().generation() << " generations\n";
//  cout << "The GA found: " << ga.statistics().bestIndividual() << "\n";

  GA1DArrayGenome<int>& best_genome =
            (GA1DArrayGenome<int>&)ga.statistics().bestIndividual();

  /* Retrieve the mapping */
  for (int service=0; service < INS.numservices; service++) {
    INS.mapping[service] = best_genome.gene(service);
  }
  
  /* Retrieve its fitness */
  float best_fitness = Objective(best_genome);

  if (best_fitness >= INS.numservers * OBJECTIVE_FEASIBLE_WEIGHT) {
    return RESOURCE_ALLOCATION_SUCCESS;
  } else {
    return RESOURCE_ALLOCATION_FAILURE;
  }

}

/* Objective function 
 *  - If the allocation is not feasible: return 2*number of overloaded servers
 *  - Otherwise, return the maximum minimum yield
 */
float Objective(GAGenome & g)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  int size = genome.length();
  float fitness = OBJECTIVE_FEASIBLE_WEIGHT * INS.numservers;

  /* Compute the server loads */
  initialize_global_server_loads();

  for (int service=0; service < genome.length(); service++) {
    int server = genome.gene(service);

    for (int j=0; j < INS.numrigid; j++)  {
      global_server_rigid_loads[server][j] += INS.rigidneeds[service][j];
    }
    for (int j=0; j < INS.numfluid; j++)  {
      global_server_fluid_loads[server][j] += INS.fluidneeds[service][j];
    }
    for (int j=0; j < INS.numfluid; j++)  {
      global_server_fluidmin_loads[server][j] += 
                          INS.slas[service] * INS.fluidneeds[service][j];
    }
  }

  /* For each overloaded server, reduce the fitness */
  int overloaded;
  for (int server=0; server < INS.numservers; server++) {
    overloaded = 0;

    for (int j=0; j < INS.numrigid; j++)  {
      if (global_server_rigid_loads[server][j] > 1.0 + EPSILON)  {
        overloaded = 1;
        break;
      }
    }

    if (!overloaded) {
      for (int j=0; j < INS.numfluid; j++)  {
        if (global_server_fluidmin_loads[server][j] > 1.0 + EPSILON)  {
          overloaded = 1;
          break;
        }
      }
    }
    if (overloaded) {
      fitness -= OBJECTIVE_FEASIBLE_WEIGHT;
    }
  }

  /* If there is but one overloaded server, we're done */
  if (fitness < OBJECTIVE_FEASIBLE_WEIGHT * INS.numservers) {
    return fitness;
  }
  
  /* We have a valid allocation, let's compute the minimum yield */

  /* Compute the mapping */
  for (int service=0; service < INS.numservices; service++) {
    int server = genome.gene(service);
    INS.mapping[service] = server;
  }
  
  /* Compute the allocation */
  for (int server=0; server<INS.numservers; server++)
    compute_allocations_given_mapping(server);

  /* Compute the minimum yield */
  float min_yield = compute_minimum_yield();
  
  /* Add the minimym yield to the fitness, scaled by the # of servers */
  fitness += min_yield * INS.numservers;

  return fitness;
}

/* Initializer function 
 *  Does not make any attempt to ensure the solution is feasible.
 */
static void myInitializer(GAGenome & g)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  int size = genome.length();

  //cout << "Creating an initial random genome\n";

  /* Come up with a random assignment */
  for (int i=0; i<size; i++) {
    genome.gene(i,GARandomInt(0,INS.numservers-1));
    //cout << "  Service #" << i << "[ ";
    //cout << "sla=" << INS.slas[i] << ", ";
    //for (int j=0; j<INS.numrigid; j++)
    //  cout << "r" << j << "=" << INS.rigidneeds[i][j] << ", ";
    //for (int j=0; j<INS.numfluid; j++)
    //  cout << "f" << j << "=" << INS.fluidneeds[i][j] << ", ";
    //
    //cout << " ] on server " << genome.gene(i) << "\n";
  }

  if (makeInitialGenomesFeasible)
    makeGenomeFeasible(genome);

  return;
}

/* Mutator function 
 *  Randomly swap services
 */
static int myMutator(GAGenome & g, float pmut)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;

  /* Perhaps we don't mutate at all */
  if ((pmut <= 0.0) || (!GAFlipCoin(pmut))) {
    return 0;
  }

  // cout << "Mutating "; print_genome(genome);

  /* Pick random services and swap them */
  int service1, service2;
  
  service1  = GARandomInt(0,INS.numservices-1);
  while ( (service2 = GARandomInt(0,INS.numservices-1)) == service1) ;

  // genome.write(&cout);
  
  int tmp = genome.gene(service1);
  genome.gene(service1, genome.gene(service2));
  genome.gene(service2, tmp);

  if (makeMutatedGenomesFeasible)
    makeGenomeFeasible(genome);

  // cout << "   into "; print_genome(genome); cout << "\n";
  return 1;
}

/* Mutator function 
 *  Randomly swap services
 */
static int myMutator2(GAGenome & g, float pmut)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;

  /* Perhaps we don't mutate at all */
  if ((pmut <= 0.0) || (!GAFlipCoin(pmut))) {
    return 0;
  }

  // cout << "Mutating "; print_genome(genome);

  /* Pick a random service and a random server */
  int service;
  int server;
  
  service  = GARandomInt(0,INS.numservices-1);
  while ( (server = GARandomInt(0,INS.numservers-1)) == genome.gene(service) );

  // genome.write(&cout);
  genome.gene(service, server);

  if (makeMutatedGenomesFeasible)
    makeGenomeFeasible(genome);

  // cout << "   into "; print_genome(genome); cout << "\n";
  return 1;
}

/* Crossover function 
 *  Randomly selects from one parent or the other
 */
int static myCrossover(const GAGenome& g1, const GAGenome& g2, 
                       GAGenome* c1, GAGenome* c2)
{

  /* Do the uniform crossover */
  GA1DArrayGenome<int>::OnePointCrossover(g1,g2,c1,c2);

  if (makeCrossedGenomesFeasible) {
    makeGenomeFeasible(*c1);
    makeGenomeFeasible(*c2);
  }
  return 1;
}

/* makeGenomeFeasible function */
static int makeGenomeFeasible(GAGenome & g)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  int j;

  // cout << "Making the genome feasible...\n";
  
  /* Compute the server loads */
  initialize_global_server_loads();
  for (int service=0; service < genome.length(); service++) {
    int server = genome.gene(service);

    for (j=0; j < INS.numrigid; j++)  {
      global_server_rigid_loads[server][j] += INS.rigidneeds[service][j];
    }

    for (j=0; j < INS.numfluid; j++)  {
      global_server_fluid_loads[server][j] += INS.fluidneeds[service][j];
    }
    for (j=0; j < INS.numfluid; j++)  {
      global_server_fluidmin_loads[server][j] += 
                          INS.slas[service] * INS.fluidneeds[service][j];
    }
  }

  /* Go through the server and make sure that all constraints are respected */
  for (int server=0; server < INS.numservers; server++ ) {
    for (j=0; j < INS.numrigid; j++) {
      //cout << "dealing with rigid dim #" << j << "\n";
      while (global_server_rigid_loads[server][j] > 1.0 + EPSILON) {
        // pick a service to evict
        //cout << "Picking a serve to evict from server " << server << "\n";
        int to_evict = pick_service_to_evict(genome, server, 'R', j); 
        //cout << "  Picked service" << to_evict << "\n";


        // adjust the server's load 
        for (int j1=0; j1 < INS.numrigid; j1++) {
          global_server_rigid_loads[server][j1] -= 
                        INS.rigidneeds[to_evict][j1];
        } 
        for (int j1=0; j1 < INS.numfluid; j1++) {
          global_server_fluid_loads[server][j1] -= 
                        INS.fluidneeds[to_evict][j1];
        }
        for (int j1=0; j1 < INS.numfluid; j1++) {
          global_server_fluidmin_loads[server][j1] -= 
                        INS.slas[to_evict] * INS.fluidneeds[to_evict][j1];
        }
       
        // evict the service
        if (evict_service(genome, to_evict, server) == -1) {
          //cout << "Cannot find server for evicted service!";
          return -1;
        }
      }
    }    
    for (j=0; j < INS.numfluid; j++) {
      //cout << "dealing with fluid dim #" << j << "\n";
      while (global_server_fluidmin_loads[server][j] > 1.0 + EPSILON) {
        // pick a service to evict
        int to_evict = pick_service_to_evict(genome, server, 'F', j); 

        // adjust the server's load 
        for (int j1=0; j1 < INS.numrigid; j1++) {
          global_server_rigid_loads[server][j1] -= 
                        INS.rigidneeds[to_evict][j1];
        } 
        for (int j1=0; j1 < INS.numfluid; j1++) {
          global_server_fluid_loads[server][j1] -= 
                        INS.fluidneeds[to_evict][j1];
        }
        for (int j1=0; j1 < INS.numfluid; j1++) {
          global_server_fluidmin_loads[server][j1] -= 
                        INS.slas[to_evict] * INS.fluidneeds[to_evict][j1];
        }
       
        // evict the service
        if (evict_service(genome, to_evict, server) == -1) {
          //cout << "Cannot find server for evicted service!";
          return -1;
        }
      }
    }    
  }

  return 0;

}

/* Pick which service to evict from an overloaded server
 * in dimension (type,j)
 */
static int pick_service_to_evict(GAGenome &g, int server, char type, int j)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  int largest = -1;
  int smallest_that_fixes_it = -1;
  int i;
  
  /* Try to find the smallest service that can fix this */
  // cout << "  Tring to find the smallest service that can fix the overload\n";
  for (i = 0; i < INS.numservices; i++) {
    if (genome.gene(i) != server)
      continue;
    if (type == 'R') {
      if (INS.rigidneeds[i][j] > 
             EPSILON + global_server_rigid_loads[server][j] - 1.0) {
        if (smallest_that_fixes_it == -1) {
          smallest_that_fixes_it = i;
        } else {
          if (INS.rigidneeds[i][j] < INS.rigidneeds[smallest_that_fixes_it][j] ) {
            smallest_that_fixes_it = i;
          }
        }
      }
    } else {
      if (INS.slas[i]*INS.fluidneeds[i][j] >  
             EPSILON + global_server_fluidmin_loads[server][j] - 1.0) {
        if (smallest_that_fixes_it == -1) {
          smallest_that_fixes_it = i;
        } else {
          if (INS.slas[i] * INS.fluidneeds[i][j] < 
              INS.slas[smallest_that_fixes_it] * 
                    INS.fluidneeds[smallest_that_fixes_it][j] ) {
            smallest_that_fixes_it = i;
          }
        }
      }
    }
  }

  /* If we found it, great */
  if (smallest_that_fixes_it != -1) {
    // cout << "    Yes! service " << smallest_that_fixes_it << "\n";
    return smallest_that_fixes_it;
  }

  /* If we didn't find it, then just get rid of the biggest job */
  // cout << "  Tring to find the largest service to evict\n";
  for (i = 0; i < INS.numservices; i++) {
    if (genome.gene(i) != server)
      continue;
    if (type == 'R') {
      if (largest == -1) {
        largest = i;
      } else {
        if (INS.rigidneeds[i][j] > INS.rigidneeds[largest][j] ) {
          largest = i;
        }
      }
    } else {
      if (largest == -1) {
        largest = i;
      } else {
        if (INS.slas[i] * INS.fluidneeds[i][j] > 
            INS.slas[largest] * INS.fluidneeds[largest][j] ) {
          largest = i;
        }
      }
    }
  }

  // cout << "    Found service " << largest << "\n";
  return largest;
}


static int evict_service(GAGenome &g, int to_evict, int from_server)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  int i;

  for (i=0; i < INS.numservers; i++) {
    if (i == from_server)
      continue;
    if (!service_can_fit_on_server_fast(to_evict, i))
      continue;

    /* found one! */

    // cout << "Eviction to server " << i << "\n";

    /* update the server loads */
    for (int j=0; j < INS.numrigid; j++) {
      global_server_rigid_loads[i][j] += 
                    INS.rigidneeds[to_evict][j];
    } 
    for (int j=0; j < INS.numfluid; j++) {
      global_server_fluid_loads[i][j] += 
                    INS.fluidneeds[to_evict][j];
    }
    for (int j=0; j < INS.numfluid; j++) {
      global_server_fluidmin_loads[i][j] += 
                    INS.slas[to_evict] * INS.fluidneeds[to_evict][j];
    }

    /* update the gene */
    genome.gene(to_evict, i);

    break;
  }

  /* Did we find one? */
  if (i == INS.numservers) {
    return -1;
  } else {
    return 0;
  }
}

static void print_genome(GAGenome &g)
{
  GA1DArrayGenome<int>& genome = (GA1DArrayGenome<int>&)g;
  cout << "[";
  for (int i=0; i < genome.length(); i++) {
    cout << " " << genome.gene(i);
  }
  cout << " ]";
  return;
}
