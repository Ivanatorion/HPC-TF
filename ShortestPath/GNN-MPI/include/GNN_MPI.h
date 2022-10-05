#ifndef GNN_MPI_H
#define GNN_MPI_H

enum GNN_AGGREGATION {GNN_AGGREGATION_NONE = 0, GNN_AGGREGATION_MIN = 1, GNN_AGGREGATION_MAX = 2, GNN_AGGREGATION_SUM = 3, GNN_AGGREGATION_AVG = 4};

typedef struct gnnS{
    MLPFeedF *msg;
    MLPFeedF *agg;
    MLPFeedF *upd;
    MLPFeedF *ro;

    int nIterations;
    int gnn_aggregation;
    int hiddenStateDim;
} GNN;

GNN* GNN_Create();
void GNN_Destroy(GNN *gnn);
void GNN_LoadFromFile(GNN* gnn, const char *file);

void GNN_Apply(GNN* gnn, Graph *graph, double *outputPerNode, int metisPartition);

void GNN_FixGraphHSDim(GNN* gnn, Graph *graph);

#endif