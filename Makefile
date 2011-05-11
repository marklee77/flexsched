PROG = flexsched

CC = gcc
CFLAGS = -I $(HOME)/progs/base/include -L $(HOME)/progs/base/lib -g

default: $(PROG) 

$(PROG): main.o util.o greedy.o linearprog.o
	$(CC) $(CFLAGS) main.o util.o greedy.o linearprog.o -o $(PROG)

%.o: %.c flexsched.h
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	/bin/rm -f *.o $(PROG)
