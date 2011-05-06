#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

double gsl_ran_2phase_uniform(
        gsl_rng *rng, double l, double m, double h, double p) {
    if (gsl_rng_uniform(rng) < p) {
        return gsl_rng_uniform(rng)*(m-l);
    } else {
        return m + gsl_rng_uniform(rng)*(h-m);
    }
}

double gsl_ran_normal(gsl_rng *rng, double median, double sigma) {
    return median + gsl_ran_gaussian(rng, sigma);
}

int round_to_nearest_int(double x) {
    double ceil_val, floor_val;
    ceil_val = ceil(x);
    floor_val = floor(x);
    if (x - floor_val > ceil_val - x)
        return (int)ceil_val;
    else
    return (int)floor_val;
}

int generate_job_size(gsl_rng *rng, int maxhosts) {
    double log_minhosts = 1;
    double log_maxhosts = log(maxhosts)/log(2);
    double p1 = 0.24;
    double p2 = 0.75;
    double l  = (double)log_minhosts - 0.2;
    double m  = (double)log_maxhosts - 2.5;
    double h  = (double)log_maxhosts;
    double p  = 0.86;
    double log_of_size;
    int size;

    /* Sequential Job? */
    if (gsl_rng_uniform(rng) < p1)
        return 1;

    log_of_size = gsl_ran_2phase_uniform(rng, l, m, h, p);

    /* Round to nearest integer? */
    if (gsl_rng_uniform(rng) < p2)
        log_of_size = (double)round_to_nearest_int(log_of_size);

    size = round_to_nearest_int(pow(2,log_of_size));

    if (size < 1)
        return 1;

    if (size > maxhosts)
        return maxhosts;

    return size;

}

int min(int a, int b) {

    if (a < b)
        return a;

    return b;

}


int main(int argc, char *argv[]) {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    int seed  = atoi(argv[1]);
    int hosts = atoi(argv[2]);
    int tasks = atoi(argv[3]);
    float slack, memvar, cpumed, cpuvar;
    float cpusigma, memsigma;
    float memmed = 0.5, memtotal = 0.0, memscale;
    int jobs = 0;
    float *jobcpu, *jobmem;
    int *jobsize;

    sscanf(argv[4], "%f", &slack);
    sscanf(argv[5], "%f", &memvar);
    sscanf(argv[6], "%f", &cpumed);
    sscanf(argv[7], "%f", &cpuvar);

    cpusigma = cpumed * cpuvar;
    memsigma = memmed * memvar;

    gsl_rng_set(rng, seed);

    jobcpu  = malloc(tasks * sizeof(float));
    jobmem  = malloc(tasks * sizeof(float));
    jobsize = malloc(tasks * sizeof(int));

    printf("%d\n", hosts);
    printf("%d\n", tasks);

    while (tasks > 0) {
        jobsize[jobs] = generate_job_size(rng, min(hosts, tasks));

        jobcpu[jobs] = 0.0;
        while (jobcpu[jobs] < 0.01 || jobcpu[jobs] > 1.0) {
            jobcpu[jobs] = gsl_ran_normal(rng, cpumed, cpusigma);
        }

        jobmem[jobs] = 0.0;
        while (jobmem[jobs] < 0.01 || jobmem[jobs] > 1.0) {
            jobmem[jobs] = gsl_ran_normal(rng, memmed, memsigma);
        }

        memtotal += (float)jobsize[jobs] * jobmem[jobs];
        
        tasks -= jobsize[jobs];
        jobs++;

    }

    memscale = (1.0 - slack) * (float)hosts / memtotal; 

    printf("%d\n", jobs);

    while (jobs > 0) {
       jobs--;
       printf("%d %f %f\n", jobsize[jobs], jobcpu[jobs], jobmem[jobs] * memscale);
    }

    return 0;
}
