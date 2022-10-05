#ifndef GRAPH_H
#define GRAPH_H

typedef struct gnn_graph_node{
    double *hiddenState; //Dimension determined in the Graph struct

    int nInEdges, nOutEdges;
    int *inEdgesIndexes, *outEdgesIndexes;
} GraphNode;

typedef struct gnn_graph_edge{
    int n1, n2;
} GraphEdge;

typedef struct gnn_graph{
    int hiddenStateDim;

    int nNodes, nEdges;
    GraphNode* nodes;
    GraphEdge* edges;

    //For Shortest Path
    int src, tgt;
} Graph;

Graph* Graph_Create();
void Graph_Destroy(Graph* graph);
void Graph_LoadFromFile(Graph* graph, const char *file, const char *partitionFile);

void Graph_Print(Graph* graph);

#endif
