#include <stdio.h>

int main(int argc, char *argv[]){
    if(argc != 3){
        fprintf(stderr, "Usage: %s <Input File> <Output File>\n", argv[0]);
        return 1;
    }

    FILE *fin = fopen(argv[1], "r");
    if(!fin){
        fprintf(stderr, "Could not open file %s\n", argv[1]);
        return 1; 
    }

    FILE *fout = fopen(argv[2], "w");
    if(!fout){
        fclose(fin);
        fprintf(stderr, "Could not open file %s\n", argv[2]);
        return 1; 
    }

    int n, m, rb;
    rb += fscanf(fin, "%d\n%d\n", &n, &m);
    fprintf(fout, "%d %d\n", n, m / 2);

    int src, tgt;
    rb += fscanf(fin, "%d\n%d\n", &src, &tgt);

    int i, curV = 0;
    int eSrc, eTgt, eWgt;
    for(i = 0; i < m; i++){
        rb += fscanf(fin, "%d %d %d\n", &eSrc, &eTgt, &eWgt);
        if(eSrc != curV){
            fprintf(fout, "\n%d ", eTgt+1);
            curV++;
        }
        else{
            fprintf(fout, "%d ", eTgt+1);
        }
    }
    fprintf(fout, "\n");

    fclose(fin);
    fclose(fout);

    return 0;
}