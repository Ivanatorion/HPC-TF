#include <mpi.h>
#include <metis.h>
#include "../include/pcheaders.h"
#include "../include/MLPFeedF.h"
#include "../include/Graph.h"
#include "../include/GNN_MPI.h"

GNN* GNN_Create(){
    GNN *gnn = (GNN *) malloc(sizeof(GNN));

    gnn->upd = NULL;
    gnn->agg = NULL;
    gnn->msg = NULL;
    gnn->ro = NULL;

    return gnn;
}

void GNN_Destroy(GNN *gnn){
    if(gnn->upd != NULL)
        MLPFeedF_Destroy(gnn->upd);
    if(gnn->agg != NULL)
        MLPFeedF_Destroy(gnn->agg);
    if(gnn->msg != NULL)
        MLPFeedF_Destroy(gnn->msg);
    if(gnn->ro != NULL)
        MLPFeedF_Destroy(gnn->ro);
    
    free(gnn);
}

void GNN_LoadFromFile(GNN* gnn, const char *file){
    char buffer[256];
    int rb = 0;
    int i, j, k;

    FILE *fp = fopen(file, "r");

    if(!fp){
        fprintf(stderr, "ERROR: Could not open file: %s\n", file);
        exit(1);
    }

    int msgLayers, aggLayers, uptLayers, readoutLayers;
    int *msgLayersDims = NULL, *aggLayersDims = NULL, *uptLayersDims = NULL, *readoutLayersDims = NULL;

    double **msgWeights = NULL, **aggWeights = NULL, **uptWeights = NULL, **readoutWeights = NULL;
    double **msgBiases = NULL, **aggBiases = NULL, **uptBiases = NULL, **readoutBiases = NULL;

    rb += fscanf(fp, "%d %d %d %d\n", &uptLayers, &aggLayers, &msgLayers, &readoutLayers);

    if(msgLayers > 0){
        msgLayersDims = (int *) malloc(sizeof(int) * msgLayers * 2); //Line - col
        msgWeights = (double **) malloc(sizeof(double *) * msgLayers);
        msgBiases = (double **) malloc(sizeof(double *) * msgLayers);
    }
    if(aggLayers > 0){
        aggLayersDims = (int *) malloc(sizeof(int) * aggLayers * 2); //Line - col
        aggWeights = (double **) malloc(sizeof(double *) * aggLayers);
        aggBiases = (double **) malloc(sizeof(double *) * aggLayers);
    }
    if(uptLayers > 0){
        uptLayersDims = (int *) malloc(sizeof(int) * uptLayers * 2); //Line - col
        uptWeights = (double **) malloc(sizeof(double *) * uptLayers);
        uptBiases = (double **) malloc(sizeof(double *) * uptLayers);
    }
    if(readoutLayers > 0){
        readoutLayersDims = (int *) malloc(sizeof(int) * readoutLayers * 2); //Line - col
        readoutWeights = (double **) malloc(sizeof(double *) * readoutLayers);
        readoutBiases = (double **) malloc(sizeof(double *) * readoutLayers);
    }

    rb += fscanf(fp, "%d\n%s\n", &gnn->nIterations, buffer);

    if(!strcmp(buffer, "min"))
        gnn->gnn_aggregation = GNN_AGGREGATION_MIN;
    else if(!strcmp(buffer, "max"))
        gnn->gnn_aggregation = GNN_AGGREGATION_MAX;
    else if(!strcmp(buffer, "mean"))
        gnn->gnn_aggregation = GNN_AGGREGATION_AVG;
    else if(!strcmp(buffer, "sum"))
        gnn->gnn_aggregation = GNN_AGGREGATION_SUM;
    else{
        fprintf(stderr, "ERROR: Unknown aggregation \"%s\"\n", buffer);
        exit(1);
    }

    rb += fscanf(fp, "%d\n", &gnn->hiddenStateDim);

    int *allActivations = malloc(sizeof(int) * (msgLayers + aggLayers + uptLayers + readoutLayers));
    for(i = 0; i < msgLayers + aggLayers + uptLayers + readoutLayers; i++){
        rb += fscanf(fp, "%s\n", buffer);
        if(!strcmp(buffer, "None"))
            allActivations[i] = MLP_ACTIVATION_NONE;
        else if(!strcmp(buffer, "sigmoid"))
            allActivations[i] = MLP_ACTIVATION_SIGMOID;
        else if(!strcmp(buffer, "relu"))
            allActivations[i] = MLP_ACTIVATION_RELU;
        else{
            fprintf(stderr, "ERROR: Unknown layer activation \"%s\"\n", buffer);
            exit(1);
        }
    }

    for(i = 0; i < uptLayers; i++){
        rb += fscanf(fp, "%d %d\n", &uptLayersDims[i * 2], &uptLayersDims[i * 2 + 1]); //Line - col

        uptWeights[i] = (double *) malloc(sizeof(double) * uptLayersDims[i * 2] * uptLayersDims[i * 2 + 1]);
        uptBiases[i] = (double *) malloc(sizeof(double) * uptLayersDims[i * 2 + 1]);

        for(j = 0; j < uptLayersDims[i * 2]; j++){
            for(k = 0; k < uptLayersDims[i * 2 + 1]; k++){
                rb += fscanf(fp, "%lf ", &uptWeights[i][j * uptLayersDims[i * 2 + 1] + k]);
            }
        }

        rb += fscanf(fp, "%d\n", &uptLayersDims[i * 2 + 1]);
        for(k = 0; k < uptLayersDims[i * 2 + 1]; k++){
            rb += fscanf(fp, "%lf ", &uptBiases[i][k]);
        }
    }

    for(i = 0; i < aggLayers; i++){
        rb += fscanf(fp, "%d %d\n", &aggLayersDims[i * 2], &aggLayersDims[i * 2 + 1]); //Line - col

        aggWeights[i] = (double *) malloc(sizeof(double) * aggLayersDims[i * 2] * aggLayersDims[i * 2 + 1]);
        aggBiases[i] = (double *) malloc(sizeof(double) * aggLayersDims[i * 2 + 1]);

        for(j = 0; j < aggLayersDims[i * 2]; j++){
            for(k = 0; k < aggLayersDims[i * 2 + 1]; k++){
                rb += fscanf(fp, "%lf ", &aggWeights[i][j * aggLayersDims[i * 2 + 1] + k]);
            }
        }

        rb += fscanf(fp, "%d\n", &aggLayersDims[i * 2 + 1]);
        for(k = 0; k < aggLayersDims[i * 2 + 1]; k++){
            rb += fscanf(fp, "%lf ", &aggBiases[i][k]);
        }
    }

    for(i = 0; i < msgLayers; i++){
        rb += fscanf(fp, "%d %d\n", &msgLayersDims[i * 2], &msgLayersDims[i * 2 + 1]); //Line - col

        msgWeights[i] = (double *) malloc(sizeof(double) * msgLayersDims[i * 2] * msgLayersDims[i * 2 + 1]);
        msgBiases[i] = (double *) malloc(sizeof(double) * msgLayersDims[i * 2 + 1]);

        for(j = 0; j < msgLayersDims[i * 2]; j++){
            for(k = 0; k < msgLayersDims[i * 2 + 1]; k++){
                rb += fscanf(fp, "%lf ", &msgWeights[i][j * msgLayersDims[i * 2 + 1] + k]);
            }
        }

        rb += fscanf(fp, "%d\n", &msgLayersDims[i * 2 + 1]);
        for(k = 0; k < msgLayersDims[i * 2 + 1]; k++){
            rb += fscanf(fp, "%lf ", &msgBiases[i][k]);
        }
    }

    for(i = 0; i < readoutLayers; i++){
        rb += fscanf(fp, "%d %d\n", &readoutLayersDims[i * 2], &readoutLayersDims[i * 2 + 1]); //Line - col

        readoutWeights[i] = (double *) malloc(sizeof(double) * readoutLayersDims[i * 2] * readoutLayersDims[i * 2 + 1]);
        readoutBiases[i] = (double *) malloc(sizeof(double) * readoutLayersDims[i * 2 + 1]);

        for(j = 0; j < readoutLayersDims[i * 2]; j++){
            for(k = 0; k < readoutLayersDims[i * 2 + 1]; k++){
                rb += fscanf(fp, "%lf ", &readoutWeights[i][j * readoutLayersDims[i * 2 + 1] + k]);
            }
        }

        rb += fscanf(fp, "%d\n", &readoutLayersDims[i * 2 + 1]);
        for(k = 0; k < readoutLayersDims[i * 2 + 1]; k++){
            rb += fscanf(fp, "%lf ", &readoutBiases[i][k]);
        }
    }

    fclose(fp);

    int atvCounter = 0;
    if(uptLayers > 0){
        int *neuronsPL = malloc(sizeof(int) * (uptLayers + 1));
        for(i = 0; i < uptLayers; i++){
            neuronsPL[i] = uptLayersDims[i * 2];
        }
        neuronsPL[uptLayers] = uptLayersDims[uptLayers * 2 - 1];

        gnn->upd = MLPFeedF_Create(uptLayers+1, neuronsPL);
        for(i = 0; i < uptLayers; i++){
            for(j = 0; j < uptLayersDims[i * 2 + 1]; j++){ //Col -> Line
                for(k = 0; k < uptLayersDims[i * 2]; k++){
                    gnn->upd->weights[i][j * uptLayersDims[i * 2] + k] = uptWeights[i][k * uptLayersDims[i * 2 + 1] + j];
                }
                gnn->upd->biases[i][j] = uptBiases[i][j];
            }
            gnn->upd->activations[i] = allActivations[atvCounter];
            atvCounter++;
            free(uptWeights[i]);
            free(uptBiases[i]);
        }

        free(uptWeights);
        free(uptBiases);
        free(uptLayersDims);
        free(neuronsPL);
    }

    if(aggLayers > 0){
        int *neuronsPL = malloc(sizeof(int) * (aggLayers + 1));
        for(i = 0; i < aggLayers; i++){
            neuronsPL[i] = aggLayersDims[i * 2];
        }
        neuronsPL[aggLayers] = aggLayersDims[aggLayers * 2 - 1];

        gnn->agg = MLPFeedF_Create(aggLayers+1, neuronsPL);
        for(i = 0; i < aggLayers; i++){
            for(j = 0; j < aggLayersDims[i * 2 + 1]; j++){ //Col -> Line
                for(k = 0; k < aggLayersDims[i * 2]; k++){
                    gnn->agg->weights[i][j * aggLayersDims[i * 2] + k] = aggWeights[i][k * aggLayersDims[i * 2 + 1] + j];
                }
                gnn->agg->biases[i][j] = aggBiases[i][j];
            }
            gnn->agg->activations[i] = allActivations[atvCounter];
            atvCounter++;
            free(aggWeights[i]);
            free(aggBiases[i]);
        }

        free(aggWeights);
        free(aggBiases);
        free(aggLayersDims);
        free(neuronsPL);
    }

    if(msgLayers > 0){
        int *neuronsPL = malloc(sizeof(int) * (msgLayers + 1));
        for(i = 0; i < msgLayers; i++)
            neuronsPL[i] = msgLayersDims[i * 2];
        neuronsPL[msgLayers] = msgLayersDims[msgLayers * 2 - 1];

        gnn->msg = MLPFeedF_Create(msgLayers+1, neuronsPL);
        for(i = 0; i < msgLayers; i++){
            for(j = 0; j < msgLayersDims[i * 2 + 1]; j++){ //Col -> Line
                for(k = 0; k < msgLayersDims[i * 2]; k++){
                    gnn->msg->weights[i][j * msgLayersDims[i * 2] + k] = msgWeights[i][k * msgLayersDims[i * 2 + 1] + j];
                }
            }
            gnn->msg->activations[i] = allActivations[atvCounter];
            atvCounter++;
            free(msgWeights[i]);
            free(msgBiases[i]);
        }

        free(msgWeights);
        free(msgBiases);
        free(msgLayersDims);
        free(neuronsPL);
    }

    if(readoutLayers > 0){
        int *neuronsPL = malloc(sizeof(int) * (readoutLayers + 1));
        for(i = 0; i < readoutLayers; i++){
            neuronsPL[i] = readoutLayersDims[i * 2];
        }
        neuronsPL[readoutLayers] = readoutLayersDims[readoutLayers * 2 - 1];

        gnn->ro = MLPFeedF_Create(readoutLayers+1, neuronsPL);
        for(i = 0; i < readoutLayers; i++){
            for(j = 0; j < readoutLayersDims[i * 2 + 1]; j++){ //Col -> Line
                for(k = 0; k < readoutLayersDims[i * 2]; k++){
                    gnn->ro->weights[i][j * readoutLayersDims[i * 2] + k] = readoutWeights[i][k * readoutLayersDims[i * 2 + 1] + j];
                }
                gnn->ro->biases[i][j] = readoutBiases[i][j];
            }
            gnn->ro->activations[i] = allActivations[atvCounter];
            atvCounter++;
            free(readoutWeights[i]);
            free(readoutBiases[i]);
        }

        free(readoutWeights);
        free(readoutBiases);
        free(readoutLayersDims);
        free(neuronsPL);
    }
    free(allActivations);
}

void GNN_FixGraphHSDim(GNN* gnn, Graph *graph){
    if(graph->hiddenStateDim == gnn->hiddenStateDim)
        return;

    int i;

    const int elemsToCopy = (graph->hiddenStateDim < gnn->hiddenStateDim) ? graph->hiddenStateDim : gnn->hiddenStateDim;
    graph->hiddenStateDim = gnn->hiddenStateDim;
    double *hsAux = (double *) malloc(sizeof(double) * elemsToCopy);

    for(i = 0; i < graph->nNodes; i++){
        memcpy(hsAux, graph->nodes[i].hiddenState, elemsToCopy * sizeof(double));
        free(graph->nodes[i].hiddenState);
        graph->nodes[i].hiddenState = (double *) malloc(sizeof(double) * gnn->hiddenStateDim);
        memset(graph->nodes[i].hiddenState, 0, sizeof(double) * gnn->hiddenStateDim);
        memcpy(graph->nodes[i].hiddenState, hsAux, sizeof(double) * elemsToCopy);
    }

    free(hsAux);
}

int getCutEdges(Graph *g, int *partitions){
    int result = 0;

    int i;
    for(i = 0; i < g->nEdges; i++)
        if(partitions[g->edges[i].n1] != partitions[g->edges[i].n2])
            result++;

    return result / 2;
}

void Partition(Graph *graph, int *graphPartition, int nParts, int metisPartition){
    int i;

    //Simple Partitioning
    for(i = 0; i < graph->nNodes; i++)
        graphPartition[i] = i % nParts;
    
    if(metisPartition == 0)
        return;

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    int edgeCut = 0;

    int *xadj = (int *) malloc(sizeof(int) * (graph->nNodes + 1) * 2);
    int *adjncy = (int *) malloc(sizeof(int) * graph->nEdges * 4);

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
}

void GNN_Apply(GNN* gnn, Graph *graph, double *outputPerNode, int metisPartition){
    GNN_FixGraphHSDim(gnn, graph);

    MPI_Status status;
    int i, j, k, l;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *graphPartition = (int *) malloc(sizeof(int) * graph->nNodes * 2);
    if(rank == 0){
        double timeStart = MPI_Wtime();
        if(size == 1)
            memset(graphPartition, 0, sizeof(int) * graph->nNodes);
        else
            Partition(graph, graphPartition, size, metisPartition);
        double timeEnd = MPI_Wtime();
        printf("Partition Time: %.4fs\n", timeEnd - timeStart);
        printf("Cut Edges: %d\n", getCutEdges(graph, graphPartition));
    }
    
    MPI_Bcast(graphPartition, graph->nNodes, MPI_INTEGER, 0, MPI_COMM_WORLD);

    const int HIDDEN_STATE_BYTES = sizeof(double) * graph->hiddenStateDim;

    //Allocate memory for the hidden node states
    double *nextHiddenStates = (double *) malloc(HIDDEN_STATE_BYTES * graph->nNodes);
    
    double *msgInput = (double *) malloc(HIDDEN_STATE_BYTES + HIDDEN_STATE_BYTES); //[source, destination]
    double *msgOutput = (double *) malloc(HIDDEN_STATE_BYTES);

    double *updateInput = (double *) malloc(HIDDEN_STATE_BYTES + HIDDEN_STATE_BYTES);
    double *updateOutput = (double *) malloc(HIDDEN_STATE_BYTES);
    
    //Message passing between processes
    int messagesRequired = 0;
    for(i = 0; i < graph->nNodes; i++){
        if(graphPartition[i] == rank){
            for(k = 0; k < graph->nodes[i].nInEdges; k++){
                if(graphPartition[graph->edges[graph->nodes[i].inEdgesIndexes[k]].n1] != rank){
                    messagesRequired++;
                }
            }
        }
    }

    double *messagesMPI = (double *) malloc(HIDDEN_STATE_BYTES * messagesRequired);
    MPI_Request *sendRequests = (MPI_Request *) malloc(sizeof(MPI_Request) * messagesRequired);
    MPI_Request *recvRequests = (MPI_Request *) malloc(sizeof(MPI_Request) * messagesRequired);

    int *messagesReceived = (int *) malloc(sizeof(int) * graph->nNodes);

    int waitingEdges, messageCounter;
    //Assuming hidden states already initialized...
    for(i = 0; i < gnn->nIterations; i++){
        //Message and Aggregation
        messageCounter = 0;
        for(j = 0; j < graph->nNodes; j++){
            messagesReceived[j] = 0;
            waitingEdges = 0;
            if(graphPartition[j] == rank){
                for(k = 0; k < graph->nodes[j].nInEdges; k++){
                    int sourceOwner = graphPartition[graph->edges[graph->nodes[j].inEdgesIndexes[k]].n1];
                    if(sourceOwner != rank){
                        MPI_Isend(graph->nodes[j].hiddenState, graph->hiddenStateDim, MPI_DOUBLE, sourceOwner, j, MPI_COMM_WORLD, &sendRequests[messageCounter]);
                        MPI_Irecv(&messagesMPI[messageCounter * graph->hiddenStateDim], graph->hiddenStateDim, MPI_DOUBLE, sourceOwner, graph->edges[graph->nodes[j].inEdgesIndexes[k]].n1, MPI_COMM_WORLD, &recvRequests[messageCounter]);
                        waitingEdges++;
                        messageCounter++;
                    }
                    else{
                        memcpy(msgInput, graph->nodes[graph->edges[graph->nodes[j].inEdgesIndexes[k]].n1].hiddenState, HIDDEN_STATE_BYTES);
                        memcpy(msgInput + graph->hiddenStateDim, graph->nodes[j].hiddenState, HIDDEN_STATE_BYTES);

                        MLPFeedF_Forward(gnn->msg, msgInput, msgOutput);

                        if(messagesReceived[j] == 0)
                            memcpy(&nextHiddenStates[graph->hiddenStateDim * j], msgOutput, HIDDEN_STATE_BYTES);
                        else{
                            switch(gnn->gnn_aggregation){
                                case GNN_AGGREGATION_MIN:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        if(nextHiddenStates[graph->hiddenStateDim * j + l] > msgOutput[l])
                                            nextHiddenStates[graph->hiddenStateDim * j + l] = msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_MAX:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        if(nextHiddenStates[graph->hiddenStateDim * j + l] < msgOutput[l])
                                            nextHiddenStates[graph->hiddenStateDim * j + l] = msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_AVG:
                                case GNN_AGGREGATION_SUM:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        nextHiddenStates[graph->hiddenStateDim * j + l] = nextHiddenStates[graph->hiddenStateDim * j + l] + msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_NONE: //TODO: Error...
                                    break;
                            }
                        }
                        messagesReceived[j]++;
                    }
                }

                if(waitingEdges == 0){
                    if(gnn->gnn_aggregation == GNN_AGGREGATION_AVG){
                        for(l = 0; l < graph->hiddenStateDim; l++)
                            nextHiddenStates[graph->hiddenStateDim * j + l] = nextHiddenStates[graph->hiddenStateDim * j + l] / graph->nodes[j].nInEdges;
                    }
                }
            }
        }
        messageCounter = 0;
        for(j = 0; j < graph->nNodes; j++){
            if(graphPartition[j] == rank){
                for(k = 0; k < graph->nodes[j].nInEdges; k++){
                    int sourceOwner = graphPartition[graph->edges[graph->nodes[j].inEdgesIndexes[k]].n1];
                    if(sourceOwner != rank){
                        MPI_Wait(&recvRequests[messageCounter], &status);

                        memcpy(msgInput, &messagesMPI[messageCounter * graph->hiddenStateDim], HIDDEN_STATE_BYTES);
                        memcpy(msgInput + graph->hiddenStateDim, graph->nodes[j].hiddenState, HIDDEN_STATE_BYTES);

                        MLPFeedF_Forward(gnn->msg, msgInput, msgOutput);

                        if(messagesReceived[j] == 0)
                            memcpy(&nextHiddenStates[graph->hiddenStateDim * j], msgOutput, HIDDEN_STATE_BYTES);
                        else{
                            switch(gnn->gnn_aggregation){
                                case GNN_AGGREGATION_MIN:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        if(nextHiddenStates[graph->hiddenStateDim * j + l] > msgOutput[l])
                                            nextHiddenStates[graph->hiddenStateDim * j + l] = msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_MAX:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        if(nextHiddenStates[graph->hiddenStateDim * j + l] < msgOutput[l])
                                            nextHiddenStates[graph->hiddenStateDim * j + l] = msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_AVG:
                                case GNN_AGGREGATION_SUM:
                                    for(l = 0; l < graph->hiddenStateDim; l++)
                                        nextHiddenStates[graph->hiddenStateDim * j + l] = nextHiddenStates[graph->hiddenStateDim * j + l] + msgOutput[l];
                                    break;
                                case GNN_AGGREGATION_NONE: //TODO: Error...
                                    break;
                            }
                        }

                        messagesReceived[j]++;
                        messageCounter++;

                        if(messagesReceived[j] == graph->nodes[j].nInEdges){
                            if(gnn->gnn_aggregation == GNN_AGGREGATION_AVG){
                                for(l = 0; l < graph->hiddenStateDim; l++)
                                    nextHiddenStates[graph->hiddenStateDim * j + l] = nextHiddenStates[graph->hiddenStateDim * j + l] / graph->nodes[j].nInEdges;
                            }
                        }
                    }
                }
            }
        }

        for(messageCounter = 0; messageCounter < messagesRequired; messageCounter++)
            MPI_Wait(&sendRequests[messageCounter], &status);

        //Update
        for(j = 0; j < graph->nNodes; j++){
            if(graphPartition[j] == rank){
                memcpy(updateInput, &nextHiddenStates[graph->hiddenStateDim * j], HIDDEN_STATE_BYTES);
                memcpy(updateInput + graph->hiddenStateDim, graph->nodes[j].hiddenState, HIDDEN_STATE_BYTES);
                MLPFeedF_Forward(gnn->upd, updateInput, updateOutput);
                memcpy(graph->nodes[j].hiddenState, updateOutput, HIDDEN_STATE_BYTES);
            }
        }
    }

    //Readout
    //Assumes 1 double per output...
    for(j = 0; j < graph->nNodes; j++){
        if(graphPartition[j] == rank){
            MLPFeedF_Forward(gnn->ro, graph->nodes[j].hiddenState, outputPerNode + j);
        }
    }

    for(j = 0; j < graph->nNodes; j++){
        if(graphPartition[j] != rank){
            if(rank == 0)
                MPI_Recv(outputPerNode + j, 1, MPI_DOUBLE, graphPartition[j], 0, MPI_COMM_WORLD, &status);     
        }
        else
            MPI_Send(outputPerNode + j, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    free(nextHiddenStates);
    free(msgInput);
    free(msgOutput);
    free(updateInput);
    free(updateOutput);

    free(graphPartition);
    free(messagesMPI);
    free(sendRequests);
    free(recvRequests);

    free(messagesReceived);
}