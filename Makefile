PROG1 = flexsched
PROG2 = flexsched2


# CFLAGS = -Wall -Wstrict-prototypes -Wcast-align -Wmissing-prototypes -Wmissing-prototypes -Wmissing-declarations -g

GLPK_INSTALLPATH=$(HOME)/PACKAGES/
LIBGA_PATH=$(HOME)/PACKAGES/galib247/ga/

CFLAGS = -O3

default: $(PROG1) 

$(PROG1): $(PROG1).o ga_compute_mapping.o
	g++ $(CFLAGS) $(PROG1).o ga_compute_mapping.o $(GLPK_INSTALLPATH)/lib/libglpk.a -o $(PROG1) -lm -L$(LIBGA_PATH) -lga

$(PROG1).o: $(PROG1).c
	gcc $(CFLAGS) -c $(PROG1).c -I$(GLPK_INSTALLPATH)/include/ -o $(PROG1).o

ga_compute_mapping.o: ga_compute_mapping.C
	g++ $(CFLAGS) -I$(LIBGA_PATH)/.. -c ga_compute_mapping.C -o ga_compute_mapping.o


$(PROG2): $(PROG2).o
	gcc $(CFLAGS) $(PROG2).o $(GLPK_INSTALLPATH)/lib/libglpk.a -o $(PROG2) -lm

$(PROG2).o: $(PROG2).c
	gcc $(CFLAGS) -c $(PROG2).c -I$(GLPK_INSTALLPATH)/include/ -o $(PROG2).o

clean:
	/bin/rm -f *.o $(PROG1) $(PROG2)
