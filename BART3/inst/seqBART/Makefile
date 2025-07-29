MPIDIR=.
LIBS=-lm 
CC=g++
MPICC=mpic++
MPICFLAGS=-DMPIBART
CCOPS=-Wall -g -O2 

OBS=tree.o funs.o bd.o rng.o
MIOBS=mi_tree.o mi_funs.o mi_bd.o mi_rng.o

#--------------------------------------------------
pbartall: mi_bart mi_bart_y mi_bart_y1 load-pbart

#--------------------------------------------------
#serial compilation
mi_bart: mi_bart.o $(MIOBS) 
	$(CC) -o mi_bart mi_bart.o $(MIOBS) $(LIBS)

mi_bart.o: mi-bart.cpp
	$(CC) $(CCOPS) -c mi-bart.cpp -o mi_bart.o

mi_bart_y: mi_bart_y.o $(MIOBS) 
	$(CC) -o mi_bart_y mi_bart_y.o $(MIOBS) $(LIBS)

mi_bart_y.o: mi-bart-y.cpp
	$(CC) $(CCOPS) -I $(MPIDIR) -c mi-bart-y.cpp -o mi_bart_y.o


mi_bart_y1: mi_bart_y1.o $(MIOBS) 
	$(CC) -o mi_bart_y1 mi_bart_y1.o $(MIOBS) $(LIBS)

mi_bart_y1.o: mi-bart-y1.cpp
	$(CC) $(CCOPS) -I $(MPIDIR) -c mi-bart-y1.cpp -o mi_bart_y1.o

mi_tree.o: tree.h tree.cpp
	$(CC) $(CCOPS) -c tree.cpp -o mi_tree.o

mi_bd.o: bd.cpp bd.h
	$(CC) $(CCOPS) -c bd.cpp -o mi_bd.o

mi_funs.o: funs.cpp funs.h
	$(CC) $(CCOPS) -c funs.cpp -o mi_funs.o

mi_rng.o: rng.C
	$(CC) $(CCOPS) -c rng.C -o mi_rng.o 

load-pbart:
	rm -f load-pbart.R
	echo "source(\"$(MPIDIR)/pbart.R\")" > load-pbart.R
	echo "bartmpicstem=\"$(MPIDIR)/\"" >> load-pbart.R

#--------------------------------------------------
clean:
	rm -f ./mi_bart
	rm -f ./mi_bart_y
	rm -f ./mi_bart_y1
	rm -f *.o
	rm -f load-pbart.R


