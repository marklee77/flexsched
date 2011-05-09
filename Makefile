PROG1 = flexsched
PROG2 = flexsched2


# CFLAGS = -Wall -Wstrict-prototypes -Wcast-align -Wmissing-prototypes -Wmissing-prototypes -Wmissing-declarations -g

GLPK_INSTALLPATH=$(HOME)/progs/base
LIBGA_PATH=$(HOME)/progs/base

CFLAGS = -g

default: $(PROG1) 

$(PROG1): $(PROG1).o ga_compute_mapping.o
	g++ $(CFLAGS) $(PROG1).o ga_compute_mapping.o -o $(PROG1) -lc -lm -L$(LIBGA_PATH)/lib -lga -L$(GLPK_INSTALLPATH)/lib -lglpk

$(PROG1).o: $(PROG1).c
	gcc $(CFLAGS) -c $(PROG1).c -I$(GLPK_INSTALLPATH)/include -o $(PROG1).o

ga_compute_mapping.o: ga_compute_mapping.C
	g++ $(CFLAGS) -I$(LIBGA_PATH)/include -c ga_compute_mapping.C -o ga_compute_mapping.o


$(PROG2): $(PROG2).o
	gcc $(CFLAGS) $(PROG2).o -o $(PROG2) -lm -L$(GLPK_INSTALLPATH)/lib -lglpk

$(PROG2).o: $(PROG2).c
	gcc $(CFLAGS) -c $(PROG2).c -I$(GLPK_INSTALLPATH)/include/ -o $(PROG2).o

clean:
	/bin/rm -f *.o $(PROG1) $(PROG2)
