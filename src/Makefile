PROG = flexsched
OBJS = main.o util.o greedy.o linutil.o linearprog.o vecutil.o vector.o hetero.o

CC = gcc
#CFLAGS = -I/opt/ilog/cplex112/include -O3 -DNO_QSORT_R -DCPLEX
#LDFLAGS = -L/opt/ilog/cplex112/bin/x86-64_debian4.0_4.1 -lcplex112
CFLAGS = -I$(HOME)/progs/base/include -g
LDFLAGS = -L$(HOME)/progs/base/lib -lglpk

default: $(PROG) 

$(PROG): $(OBJS) 
	$(CC) $(LDFLAGS) $(OBJS) -o $(PROG)

linutil.o: linutil.c linearprog.h flexsched.h Makefile
%.o: %.c flexsched.h Makefile
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	/bin/rm -f *.o $(PROG)
