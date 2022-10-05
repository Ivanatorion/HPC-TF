#include "../include/pcheaders.h"
#include "../include/MLPFeedF.h"

double _mlp_sigmoid(const double n){
    return 1.0 / (1.0 + exp(-n));
}

double _mlp_random_weight(){
    return 1.0 - ((double) rand() / ((double) RAND_MAX / 2.0));
}

MLPFeedF* MLPFeedF_Create(int nLayers, int* neuronsPerLayer){
    if(nLayers < 2){
        fprintf(stderr, "ERROR: MPLFeedF must have at least 2 layers (in / out)\n");
        exit(1);
    }

    MLPFeedF *mlpff = (MLPFeedF *) malloc(sizeof(MLPFeedF));
    
    mlpff->nLayers = nLayers;
    mlpff->neuronsPerLayer = (int *) malloc(sizeof(int) * nLayers);
    memcpy(mlpff->neuronsPerLayer, neuronsPerLayer, sizeof(int) * nLayers);

    mlpff->activations = (int *) malloc(sizeof(int) * (nLayers - 1));
    memset(mlpff->activations, 0, sizeof(int) * (nLayers - 1)); //Activations set to NONE

    mlpff->weights = (double **) malloc(sizeof(double *) * (nLayers - 1));
    mlpff->biases = (double **) malloc(sizeof(double *) * (nLayers - 1));

    int i;
    for(i = 0; i < mlpff->nLayers - 1; i++){
        mlpff->weights[i] = (double *) malloc(sizeof(double) * mlpff->neuronsPerLayer[i] * mlpff->neuronsPerLayer[i+1]);
        mlpff->biases[i] = (double *) malloc(sizeof(double) * mlpff->neuronsPerLayer[i+1]);
        memset(mlpff->weights[i], 0, sizeof(double) * mlpff->neuronsPerLayer[i] * mlpff->neuronsPerLayer[i+1]);
        memset(mlpff->biases[i], 0, sizeof(double) * mlpff->neuronsPerLayer[i+1]);
    }

    return mlpff;
}

void MLPFeedF_Destroy(MLPFeedF *mlpff){
    int i;
    for(i = 0; i < mlpff->nLayers - 1; i++){
        free(mlpff->weights[i]);
        free(mlpff->biases[i]);
    }

    free(mlpff->weights);
    free(mlpff->biases);
    free(mlpff->neuronsPerLayer);
    free(mlpff->activations);
    free(mlpff);
}

void MLPFeedF_RandomizeWeights(MLPFeedF *mlpff){
    int i, j, l;
    for(l = 0; l < mlpff->nLayers - 1; l++){
        for(i = 0; i < mlpff->neuronsPerLayer[l + 1]; i++){
            mlpff->biases[l][i] = _mlp_random_weight();
            for(j = 0; j < mlpff->neuronsPerLayer[l]; j++)
                mlpff->weights[l][i * mlpff->neuronsPerLayer[l] + j] = _mlp_random_weight();
        }
    }
}

void MLPFeedF_Forward(MLPFeedF *mlpff, double *input, double *output){
    int i, j, l;
    int maxNPL = 0;
    for(i = 0; i < mlpff->nLayers; i++)
        if(mlpff->neuronsPerLayer[i] > maxNPL)
            maxNPL = mlpff->neuronsPerLayer[i];

    double *auxI = malloc(maxNPL * sizeof(double));
    double *auxO = malloc(maxNPL * sizeof(double));

    memcpy(auxI, input, mlpff->neuronsPerLayer[0] * sizeof(double));

    for(l = 0; l < mlpff->nLayers - 1; l++){
        for(i = 0; i < mlpff->neuronsPerLayer[l+1]; i++){
            auxO[i] = mlpff->biases[l][i];

            for(j = 0; j < mlpff->neuronsPerLayer[l]; j++)
                auxO[i] = auxO[i] + auxI[j] * mlpff->weights[l][i * mlpff->neuronsPerLayer[l] + j];
            
            switch(mlpff->activations[l]){
                case MLP_ACTIVATION_RELU:
                    auxO[i] = (auxO[i] < 0.0) ? 0.0 : auxO[i];
                    break;
                case MLP_ACTIVATION_SIGMOID:
                    auxO[i] = _mlp_sigmoid(auxO[i]);
                    break;
                case MLP_ACTIVATION_NONE:
                    break;
            }
        }
        memcpy(auxI, auxO, mlpff->neuronsPerLayer[l+1] * sizeof(double));
    }

    memcpy(output, auxO, mlpff->neuronsPerLayer[mlpff->nLayers - 1] * sizeof(double));

    free(auxI);
    free(auxO);
}

void MLPFeedF_Print(MLPFeedF *mlpff){
    int i, j, l;
    for(l = 0; l < mlpff->nLayers - 1; l++){
        for(i = 0; i < mlpff->neuronsPerLayer[l + 1]; i++){
            for(j = 0; j < mlpff->neuronsPerLayer[l]; j++)
                printf("%.4f ", mlpff->weights[l][i * mlpff->neuronsPerLayer[l] + j]);
            printf("\n");
        }
        printf("\n");
        for(i = 0; i < mlpff->neuronsPerLayer[l + 1]; i++){
            printf("%.4f ", mlpff->biases[l][i]);
        }
        printf("\n\n");
    }
}