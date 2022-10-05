#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct lline{
    int rank;
    double start;
    double end;
    double duration;
} LogLine;

typedef struct outLine{
    int rank;
    double mpiTime;
    double fullTime;
    double compTime;
    double commRatio;

    double maxEnd;
    double minStart;
} OutLine;

LogLine* readFile(char file[], int *nLines, int *nRanks){
    char ignoredStringField[3][256];
    char stringRank[256];
    double ignoredDoubleField;

    FILE *fp = fopen(file, "r");

    if(!fp){
        fprintf(stderr, "Could not open file: %s\n", file);
        exit(1);
    }

    *nLines = 0;
    while(!feof(fp)){
        if(fgetc(fp) == '\n')
            *nLines = *nLines + 1;
    }

    fclose(fp);

    LogLine *lines = (LogLine *) malloc(sizeof(LogLine) * (*nLines));

    fp = fopen(file, "r");
    int i, j, rb = 0;
    *nRanks = 0;
    for(i = 0; i < *nLines; i++){
        rb += fscanf(fp, "%s %s %s %lf, %lf, %lf, %lf, %s\n", ignoredStringField[0], stringRank, ignoredStringField[1], &lines[i].start, &lines[i].end, &lines[i].duration, &ignoredDoubleField, ignoredStringField[2]);
        lines[i].rank = 0;
        j = 4;
        while(stringRank[j] != ','){
            lines[i].rank = lines[i].rank * 10 + (stringRank[j] - '0');
            j++;
        }
        if(lines[i].rank > *nRanks)
            *nRanks = lines[i].rank;
    }
    fclose(fp);

    *nRanks = *nRanks + 1;

    return lines;
}

int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Usage: %s <InputFile.csv> <OutputFile.csv>\n", argv[0]);
        return 0;
    }

    int nLines, nRanks;
    LogLine *lines = readFile(argv[1], &nLines, &nRanks);

    int i = 0;

    OutLine *outLines = (OutLine *) malloc(sizeof(OutLine) * nRanks);
    for(i = 0; i < nRanks; i++){
        outLines[i].rank = i;
        outLines[i].fullTime = 0.0;
        outLines[i].mpiTime = 0.0;
        outLines[i].compTime = 0.0;
        outLines[i].commRatio = 0.0;

        outLines[i].maxEnd = -1.0;
        outLines[i].minStart = -1.0;
    }

    for(i = 0; i < nLines; i++){
        int rank = lines[i].rank;
        outLines[rank].mpiTime = outLines[rank].mpiTime + lines[i].duration;
        if(outLines[rank].minStart < 0.0 || lines[i].start < outLines[rank].minStart)
            outLines[rank].minStart = lines[i].start;
        if(outLines[rank].maxEnd < 0.0 || lines[i].end > outLines[rank].maxEnd)
            outLines[rank].maxEnd = lines[i].end;
    }

    for(i = 0; i < nRanks; i++){
        outLines[i].fullTime = outLines[i].maxEnd - outLines[i].minStart;
        outLines[i].compTime = outLines[i].fullTime - outLines[i].mpiTime;
        outLines[i].commRatio = (outLines[i].mpiTime * 100.0) / outLines[i].fullTime;
    }

    FILE *fp = fopen(argv[2], "w");
    fprintf(fp, "Rank,MPI.Time,Full.Time,Comp.Time,Comm.Ratio\n");

    for(i = 0; i < nRanks; i++){
        fprintf(fp, "%d,%.4f,%.4f,%.4f,%.2f\n", outLines[i].rank, outLines[i].mpiTime, outLines[i].fullTime, outLines[i].compTime, outLines[i].commRatio);
    }
    
    fclose(fp);

    free(lines);
    free(outLines);
    return 0;
}