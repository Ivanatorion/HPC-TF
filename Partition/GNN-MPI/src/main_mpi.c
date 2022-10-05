#include "../include/pcheaders.h"
#include <mpi.h>
#include "../include/Graph.h"
#include "../include/MLPFeedF.h"
#include "../include/GNN.h"

double cpu_time(){
   return ( double ) clock() / ( double ) CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]){
    int rank, size;
    double time1 = cpu_time();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc != 4){
        if(rank == 0)
            fprintf(stderr, "Usage: %s <InputGraph> <InputPartition> <UseMetisPartition? (0 / 1)>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    Graph *g = Graph_Create();
    Graph_LoadFromFile(g, argv[1], argv[2]);

    GNN* gnn = GNN_Create();
    GNN_LoadFromFile(gnn, "GNN/model_weights.dat");

    double *outputPerNode = malloc(sizeof(double) * g->nNodes);

    GNN_Apply(gnn, g, outputPerNode, (argv[3][0] == '1') ? 1 : 0);

    /*
    int i;
    if(rank == 0){
        printf("Output:\n");
        for(i = 0; i < g->nNodes; i++)
            printf("%d: %.4f\n", i+1, outputPerNode[i]);
    }
    */

    free(outputPerNode);

    GNN_Destroy(gnn);
    Graph_Destroy(g);

    double time2 = cpu_time();
    if(rank == 0)
        printf("%.4f\n", time2 - time1);

    MPI_Finalize();
    return 0;
}
