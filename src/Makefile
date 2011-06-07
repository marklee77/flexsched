PROG = flexsched
OBJS = main.o util.o #greedy.o linearprog.o vecutil.o vector.o hetero.o

CC = gcc
CFLAGS = -I$(HOME)/progs/base/include -g #-DDEBUG
LDFLAGS = -L$(HOME)/progs/base/lib -lglpk

default: $(PROG) 

$(PROG): $(OBJS) 
	$(CC) $(LDFLAGS) $(OBJS) -o $(PROG)

%.o: %.c flexsched.h Makefile
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	/bin/rm -f *.o $(PROG)
