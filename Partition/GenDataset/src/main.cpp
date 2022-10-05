#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>

#include "../include/Fila.h"
#include "../include/Grafo.h"

int strToInt(char *s){
  int result = 0;
  while(*s != '\0'){
    if(*s < '0' || *s > '9')
      return -1;

    result = result * 10 + (*s - '0');
    s++;
  }
  return result;
}

void parseArgs(int argc, char* argv[], bool *help, bool *verbose,  bool *randD, float *randomP, int *rSeed, int *timeout, char oFile[], char iFile[], char pFile[]){
  int i = 1;
  while(i < argc){
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
      *help = true;
    }
    else if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")){
      *verbose = true;
    }
    else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--seed")){
      i++;
      if(i < argc)
        *rSeed = strToInt(argv[i]);
    }
    else if(!strcmp(argv[i], "-rp") || !strcmp(argv[i], "--randomP")){
      i++;
      if(i < argc){
        *randomP = atof(argv[i]);
        if(*randomP < -0.00001 || *randomP > 1.00001){
          printf("Invalid p: %s\n", argv[i]);
          exit(0);
        }
      }
    }
    else if(!strcmp(argv[i], "--timeout")){
      i++;
      if(i < argc){
        *timeout = strToInt(argv[i]);
        if(*timeout < 1){
          printf("Invalid timeout: %s\n", argv[i]);
          exit(0);
        }
      }
    }
    else if(!strcmp(argv[i], "-if") || !strcmp(argv[i], "--inputFile")){
      i++;
      if(i < argc)
        strcpy(iFile, argv[i]);
    }
    else if(!strcmp(argv[i], "-of") || !strcmp(argv[i], "--outputFile")){
      i++;
      if(i < argc)
        strcpy(oFile, argv[i]);
    }
    else if(!strcmp(argv[i], "-pf") || !strcmp(argv[i], "--partitionFile")){
      i++;
      if(i < argc)
        strcpy(pFile, argv[i]);
    }
    i++;
  }
}

void showHelp(){
  printf("Use: ./prog [Args]\n\nArgs:\n");
  printf("(--help | -h)                : Displays help message\n");
  printf("(--verbose | -v)             : Print execution information\n");
  printf("(--partitionFile | -pf [F]   : Initial partition file\n");
  printf("(--seed [S] | -s [S])        : Use S as the seed for random numbers\n");
  printf("(--randomP [P] | -rp [P])    : Use random prob P (0.0 - 1.0)\n");
  printf("(--timeout [T])              : Timeout T seconds\n");
  printf("(--inputFile [F] | -if [F])  : Input file F\n");
  printf("(--outputFile [F] | -of [F]) : Output file F\n");
}

int main(int argc, char *argv[]){

  int rSeed = 1024;
  int repeat = 1;
  int timeout = 300;
  bool help = false;
  bool verbose = false;
  bool randD = false;

  float randomP = 0.05;

  char iFile[256] = "NO_FILE";
  char oFile[256] = "NO_FILE";
  char pFile[256] = "NO_FILE";

  parseArgs(argc, argv, &help, &verbose, &randD, &randomP, &rSeed, &timeout, oFile, iFile, pFile);

  if(help){
    showHelp();
    return 0;
  }

  if(!strcmp(iFile, "NO_FILE")){
    fprintf(stderr, "ERROR: No input file\n");
    return 1;
  }

  if(!strcmp(pFile, "NO_FILE")){
    fprintf(stderr, "ERROR: No partition file\n");
    return 1;
  }

  srand(rSeed);

  Grafo *g1 = new Grafo(iFile, pFile);

  int initialValue = g1->getValorParticao();

  g1->randomizada(verbose, randomP, timeout);

  FILE *fp = NULL;

  if(strcmp(oFile, "NO_FILE")){
    fp = fopen(oFile, "w");
    std::vector<bool> particao = g1->getParticao();
    for(int i = 0; i < (int) particao.size(); i++)
      fprintf(fp, "%d\n", particao[i] ? 1 : 0);
    fclose(fp);
  }
  
  printf("Iterations: %lld\nTime (us): %lld\nInitial Value: %d\nFinal Value: %d\n", g1->getIteracoes(), g1->getTempoMicro(), initialValue, g1->getValorParticao());

  delete g1;
  return 0;
}
