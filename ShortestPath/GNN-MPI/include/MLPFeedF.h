#ifndef MLP_FEED_FORWARD_H
#define MLP_FEED_FORWARD_H

enum MLP_ACTIVATION {MLP_ACTIVATION_NONE = 0, MLP_ACTIVATION_RELU = 1, MLP_ACTIVATION_SIGMOID = 2};

typedef struct mlpff{
    int nLayers;
    int *neuronsPerLayer;

    double **weights; //Each points to 2d matrix
    double **biases;  //Each points to 1d array

    int *activations;
} MLPFeedF;

MLPFeedF* MLPFeedF_Create(int nLayers, int* neuronsPerLayer);
void MLPFeedF_Destroy(MLPFeedF *mlpff);
void MLPFeedF_RandomizeWeights(MLPFeedF *mlpff);
void MLPFeedF_Forward(MLPFeedF *mlpff, double *input, double *output);
void MLPFeedF_Print(MLPFeedF *mlpff);

#endif