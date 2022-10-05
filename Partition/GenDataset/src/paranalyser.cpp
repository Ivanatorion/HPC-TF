#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include "../include/Fila.h"
#include "../include/Grafo.h"

#define N_GRAPHS 200

int main(int argc, char *argv[]){
    char buffer[256];
    char buffer2[256];

    long long int metisCE = 0;
    long long int searchCE = 0;

    for(int i = 0; i < N_GRAPHS; i++){
        sprintf(buffer, "Graphs/ConvertedGraphs/G%d.txt", i);
        sprintf(buffer2, "Graphs/Partitions/G%d.2", i);
        Grafo g(buffer, buffer2);

        sprintf(buffer2, "Results/G%d.out", i);
        Grafo g2(buffer, buffer2);

        metisCE = metisCE + g.getValorParticao();
        searchCE = searchCE + g2.getValorParticao();
    }

    printf("%0.2f\n", ((double) searchCE) / metisCE);
    return 0;
}