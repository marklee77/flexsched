#!/usr/bin/python

def main(argv=None):
    from argparse import ArgumentParser
    import random
    parser = ArgumentParser(description="Generat a random problem instance.")
    parser.add_argument('-s', '--seed', type=int, help='seed for random number generator')
    parser.add_argument('-r', '--resources', type=int, help='number of resources')
    parser.add_argument('-n', '--nodes', type=int, help='number of nodes')
    parser.add_argument('-t', '--tasks', type=int, help='number of tasks')
    parser.add_argument('-m', '--mean', type=float, help='node mean resource')
    parser.add_argument('-c', '--cov', type=float, help='node resource cov')

    args = parser.parse_args()

    random.seed(args.seed) 

    print args.resources
    print args.nodes
    print args.tasks

    for node in range(args.nodes):
        # truncate between 0.001 and 1.000
        resources = [ max(0.001, min(random.normalvariate(args.mean, args.mean * args.cov), 1.0))
                      for x in range(args.resources)]

        # unit resources equal aggregate for now...
        for resource in resources:
            print "%0.3f" % (resource),
        for resource in resources:
            print "%0.3f" % (resource),

        print

    for task in range(args.tasks):

        # ignoring units for this...
        for resource in range(args.resources):
            print "0.000",
        for resource in range(args.resources):
            print "0.000",

        # do we get viable problems with just uniform random?
        for resource in resources:
            print "%0.3f" % (random.random()),
        for resource in resources:
            print "%0.3f" % (random.random()),
        
        print

if __name__ == "__main__":
    main()
