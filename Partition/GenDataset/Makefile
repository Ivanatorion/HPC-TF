all: prog createDataset

createDataset: src/createDataset.c
	gcc -O3 -o createDataset src/createDataset.c

paranalyser: bin/Grafo.o bin/Fila.o bin/paranalyser.o
	g++ -o paranalyser bin/Grafo.o bin/Fila.o bin/paranalyser.o -lm

bin/paranalyser.o: src/paranalyser.cpp include/Fila.h include/Grafo.h
	g++ -O3 -o bin/paranalyser.o -c src/paranalyser.cpp

prog: bin/main.o bin/Grafo.o bin/Fila.o
	g++ -std=c++11 -O3 -o prog bin/*.o -lm

bin/main.o: src/main.cpp include/Grafo.h include/Fila.h
	g++ -std=c++11 -O3 -o bin/main.o -c src/main.cpp

bin/Grafo.o: src/Grafo.cpp include/Grafo.h include/Fila.h
	g++ -std=c++11 -O3 -o bin/Grafo.o -c src/Grafo.cpp

bin/Fila.o: src/Fila.cpp include/Fila.h
	g++ -std=c++11 -O3 -o bin/Fila.o -c src/Fila.cpp

clean:
	rm -f prog createDataset paranalyser Graphs/mconvert
	rm -rf bin/*.o
