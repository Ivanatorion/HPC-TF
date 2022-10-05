OPTIMIZATION_FLAG = -O2

all: gnn gnnmpi gnnmpi-nopartition AuxPrograms/mconvert AuxPrograms/paranalyser AuxPrograms/mpitimecalc

include/pcheaders.h.gch: include/pcheaders.h
	gcc $(OPTIMIZATION_FLAG) -o include/pcheaders.h.gch -c include/pcheaders.h

AuxPrograms/mconvert: src/AuxPrograms/mconvert.c
	gcc $(OPTIMIZATION_FLAG) -o AuxPrograms/mconvert src/AuxPrograms/mconvert.c -lm -Wall

AuxPrograms/paranalyser: bin/paranalyser.o bin/Graph.o
	gcc $(OPTIMIZATION_FLAG) -o AuxPrograms/paranalyser bin/Graph.o bin/paranalyser.o -lm -lmetis -Wall

AuxPrograms/mpitimecalc: src/AuxPrograms/mpitimecalc.c
	gcc $(OPTIMIZATION_FLAG) -o AuxPrograms/mpitimecalc src/AuxPrograms/mpitimecalc.c -lm -Wall

bin/paranalyser.o: include/pcheaders.h src/AuxPrograms/paranalyser.c include/Graph.h
	gcc $(OPTIMIZATION_FLAG) -o bin/paranalyser.o -c src/AuxPrograms/paranalyser.c -lm -Wall

gnn: bin/main.o bin/MLPFeedF.o bin/GNN.o bin/Graph.o
	gcc $(OPTIMIZATION_FLAG) -o gnn bin/main.o bin/MLPFeedF.o bin/GNN.o bin/Graph.o -static -lm -Wall

gnnmpi: bin/main_mpi.o bin/MLPFeedF_mpi.o bin/GNN_mpi.o bin/Graph_mpi.o
	mpicc $(OPTIMIZATION_FLAG) -o gnnmpi bin/main_mpi.o bin/MLPFeedF_mpi.o bin/GNN_mpi.o bin/Graph_mpi.o -laky -lrastro -lm -lmetis -Wall

gnnmpi-nopartition: bin/main_mpi.o bin/MLPFeedF_mpi.o bin/GNN_mpi-nopartition.o bin/Graph_mpi.o
	mpicc $(OPTIMIZATION_FLAG) -o gnnmpi-nopartition bin/main_mpi.o bin/MLPFeedF_mpi.o bin/GNN_mpi-nopartition.o bin/Graph_mpi.o -lm -Wall

bin/main.o: include/pcheaders.h include/Graph.h include/MLPFeedF.h include/GNN.h src/main.c
	gcc $(OPTIMIZATION_FLAG) -o bin/main.o -c src/main.c -Wall

bin/MLPFeedF.o: include/pcheaders.h include/MLPFeedF.h src/MLPFeedF.c
	gcc $(OPTIMIZATION_FLAG) -o bin/MLPFeedF.o -c src/MLPFeedF.c -Wall

bin/GNN.o: include/pcheaders.h include/Graph.h include/MLPFeedF.h include/GNN.h src/GNN.c
	gcc $(OPTIMIZATION_FLAG) -o bin/GNN.o -c src/GNN.c -Wall

bin/Graph.o: include/pcheaders.h include/Graph.h src/Graph.c
	gcc $(OPTIMIZATION_FLAG) -o bin/Graph.o -c src/Graph.c -Wall

bin/main_mpi.o: include/pcheaders.h include/Graph.h include/MLPFeedF.h include/GNN_MPI.h src/main_mpi.c
	mpicc $(OPTIMIZATION_FLAG) -o bin/main_mpi.o -c src/main_mpi.c -Wall

bin/MLPFeedF_mpi.o: include/pcheaders.h include/MLPFeedF.h src/MLPFeedF.c
	mpicc $(OPTIMIZATION_FLAG) -o bin/MLPFeedF_mpi.o -c src/MLPFeedF.c -Wall

bin/GNN_mpi.o: include/pcheaders.h include/Graph.h include/MLPFeedF.h include/GNN_MPI.h src/GNN_mpi.c
	mpicc $(OPTIMIZATION_FLAG) -o bin/GNN_mpi.o -c src/GNN_mpi.c -Wall

bin/GNN_mpi-nopartition.o: include/pcheaders.h include/Graph.h include/MLPFeedF.h include/GNN_MPI.h src/GNN_mpi-nopartition.c
	mpicc $(OPTIMIZATION_FLAG) -o bin/GNN_mpi-nopartition.o -c src/GNN_mpi-nopartition.c -Wall

bin/Graph_mpi.o: include/pcheaders.h include/Graph.h src/Graph.c
	mpicc $(OPTIMIZATION_FLAG) -o bin/Graph_mpi.o -c src/Graph.c -Wall

clean:
	rm -f include/pcheaders.h.gch
	rm -f gnn gnnmpi gnnmpi-nopartition AuxPrograms/mconvert AuxPrograms/paranalyser AuxPrograms/mpitimecalc
	rm -rf bin/*.o
