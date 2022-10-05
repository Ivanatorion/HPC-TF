#include "../include/pcheaders.h"
#include "../include/Graph.h"
#include "../include/MLPFeedF.h"
#include "../include/GNN.h"

int main(int argc, char *argv[]){
    int i;

    if(argc != 2){
        fprintf(stderr, "Usage: %s <InputGraph>\n", argv[0]);
        return 1;
    }

    Graph *g = Graph_Create();
    Graph_LoadFromFile(g, argv[1]);

    GNN* gnn = GNN_Create();
    GNN_LoadFromFile(gnn, "GNN/model_weights.dat");

    double *outputPerNode = malloc(sizeof(double) * g->nNodes);

    GNN_Apply(gnn, g, outputPerNode);
    
    printf("Output:\n");
    for(i = 0; i < g->nNodes; i++)
        printf("%d: %.4f\n", i, outputPerNode[i]);

    printf("SP:\n");
    for(i = 0; i < g->nNodes; i++)
        if(outputPerNode[i] > 0.5)
            printf("%d: %.4f\n", i, outputPerNode[i]);

    free(outputPerNode);

    GNN_Destroy(gnn);
    Graph_Destroy(g);
    return 0;
}
