#include <metis.h>
#include "../../include/pcheaders.h"
#include "../../include/Graph.h"

int intFromStr(char str[]){
    int result = 0;
    int i = 0;
    while(str[i] != '\0'){
        if(str[i] >= '0' && str[i] <= '9')
            result = result * 10 + (str[i] - '0');
        else
            return -1;
        i++;
    }
    return result;
}

void getPartitions(char file[], int *partitions, int nNodes){
    char buffer[256];
    
    FILE *fp = fopen(file, "r");

    if(!fp){
        fprintf(stderr, "Could not open file: %s\n", file);
        exit(1);
    }

    char c;
    int i, auxI;
    for(i = 0; i < nNodes; i++){
        auxI = 0;
        c = fgetc(fp);
        while(!feof(fp) && c != '\n'){
            buffer[auxI] = c;
            auxI++;
            c = fgetc(fp);
        }
        buffer[auxI] = '\0';
        partitions[i] = intFromStr(buffer);
    }

    fclose(fp);
}

int getCutEdges(Graph *g, int *partitions){
    int result = 0;

    int i;
    for(i = 0; i < g->nEdges; i++)
        if(partitions[g->edges[i].n1] != partitions[g->edges[i].n2])
            result++;

    return result / 2;
}

int getSimplePartitionEdges(Graph *g, int nParts){
    int result = 0;

    int i;
    for(i = 0; i < g->nEdges; i++)
        if(g->edges[i].n1 % nParts != g->edges[i].n2 % nParts)
            result++;

    return result / 2;
}

int Partition(Graph *graph, int *graphPartition, int nParts){
    int i;

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    int edgeCut = 0;

    int *xadj = (int *) malloc(sizeof(int) * (graph->nNodes + 1));
    int *adjncy = (int *) malloc(sizeof(int) * graph->nEdges * 2);

    int cv = 0;
    xadj[0] = 0;
    for(i = 0; i < graph->nEdges; i++){
        if(graph->edges[i].n1 != cv){
            cv++;
            xadj[cv] = i;
        }
        adjncy[i] = graph->edges[i].n2;
    }
    xadj[graph->nNodes] = graph->nEdges * 2;

    int nCon = 1;
    METIS_PartGraphRecursive(&graph->nNodes, &nCon, xadj, adjncy, NULL, NULL, NULL, &nParts, NULL, NULL, options, &edgeCut, graphPartition);

    free(xadj);
    free(adjncy);

    return edgeCut;
}

int main(int argc, char *argv[]){
    if(argc != 4){
        fprintf(stderr, "Usage: %s <Nparts> <Input Graph> <Partition File>\n", argv[0]);
        return 1;
    }

    int nParts = intFromStr(argv[1]);
    if(nParts < 1){
        fprintf(stderr, "Invalid number of partitions: %s\n", argv[1]);
        return 1;
    }

    Graph *g = Graph_Create();
    Graph_LoadFromFile(g, argv[2]);

    int *nodePartitions = (int *) malloc(sizeof(int) * g->nNodes);
    getPartitions(argv[3], nodePartitions, g->nNodes);

    int nCutEdges = getCutEdges(g, nodePartitions);
    int nSimplePartitionEdges = getSimplePartitionEdges(g, nParts);
    Partition(g, nodePartitions, nParts);
    int metisCEdgeCut = getCutEdges(g, nodePartitions);
    printf("Cut Edges: %d\nSimple Partition: %d\nMETIS_PartGraphRecursive: %d\n", nCutEdges, nSimplePartitionEdges, metisCEdgeCut);

    Graph_Destroy(g);
    free(nodePartitions);
    return 0;
}