#include <stdio.h>

int main(int argc, char *argv[]){
    char buffer[512];

    FILE *fp = fopen("dataset.json", "w");
    fprintf(fp, "[");

    for(int i = 0; i < 200; i++){
        sprintf(buffer, "Graphs/GraphSamples/Graph_%d.txt", i);
        FILE *graphFile = fopen(buffer, "r");
        if(!graphFile){
            fprintf(stderr, "ERROR: Could not open file: %s\n", buffer);
            return 1;
        }
        sprintf(buffer, "Results/G%d.out", i);
        FILE *partitionFile = fopen(buffer, "r");
        if(!partitionFile){
            fprintf(stderr, "ERROR: Could not open file: %s\n", buffer);
            return 1;
        }
        sprintf(buffer, "Graphs/Partitions/G%d.2", i);
        FILE *partitionFileStart = fopen(buffer, "r");
        if(!partitionFileStart){
            fprintf(stderr, "ERROR: Could not open file: %s\n", buffer);
            return 1;
        }

        int rb = 0;
        int n, m;
        rb += fscanf(graphFile, "%d\n%d\n", &n, &m);

        fprintf(fp, "{\"directed\" : true, \"multigraph\" : false, \"graph\": {}, \"nodes\": [");
        for(int j = 0; j < n; j++){
            int partIn, partOut;
            rb += fscanf(partitionFileStart, "%d\n", &partIn);
            rb += fscanf(partitionFile, "%d\n", &partOut);
            fprintf(fp, "{\"entity\": \"node\", \"id\": %d, \"pi\": %d, \"po\": %d}", j, partIn, partOut);

            if(j < n - 1)
                fprintf(fp, ", ");
        }
        fprintf(fp, "], \"links\": [");
        for(int j = 0; j < m; j++){
            int src, tgt;
            rb += fscanf(graphFile, "%d %d\n", &src, &tgt);
            fprintf(fp, "{\"source\": %d, \"target\": %d}", src, tgt);

            if(j < m - 1)
                fprintf(fp, ", ");
        }
        fprintf(fp, "]}");

        if(i < 200 - 1)
            fprintf(fp, ", ");

        fclose(graphFile);
        fclose(partitionFile);
    }

    fprintf(fp, "]");
    fclose(fp);

    return 0;
}
