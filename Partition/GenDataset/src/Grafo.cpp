#include <vector>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstdint>
#include <queue>
#include <unordered_set>
#include <algorithm>

#include "../include/Fila.h"
#include "../include/Grafo.h"

#define INT_MAX 2000000000

Grafo::Grafo(const char *filePath, const char *partitionFile){
  FILE *fp = fopen(filePath, "r");

  if(!fp){
    printf("Erro - Nao pode abir o arquivo: %s\n", filePath);
    exit(1);
  }

  int n, m;
  char c;
  char aux[256];
  int i;

  c = fgetc(fp);
  while(c == ' '){
    c = fgetc(fp);
  }
  i = 0;
  while(c >= '0' && c <= '9'){
    aux[i] = c;
    i++;
    c = fgetc(fp);
  }
  aux[i] = '\0';
  n = atoi(aux);

  while(!(c >= '0' && c <= '9'))
    c = fgetc(fp);
  i = 0;
  while(c >= '0' && c <= '9'){
    aux[i] = c;
    i++;
    c = fgetc(fp);
  }
  aux[i] = '\0';
  m = atoi(aux);

  this->listaAdj.resize(n);

  for(int x = 0; x < n; x++){
    int v1, v2;
    Aresta aresta;

    bool endLine = false;

    v1 = x;
    while(!endLine){
      c = fgetc(fp);
      while(!(c >= '0' && c <= '9') && c != '\n' && !feof(fp))
        c = fgetc(fp);

      if(c == '\n' || feof(fp)){
        endLine = true;
      }
      else{
        i = 0;
        while(c >= '0' && c <= '9'){
          aux[i] = c;
          i++;
          c = fgetc(fp);
        }
        aux[i] = '\0';
        v2 = atoi(aux) - 1;

        aresta.origem = v1;
        aresta.destino = v2;
        this->listaAdj[aresta.origem].push_back(aresta);

        if(c == '\n')
          endLine = true;
      }
    }
  }

  fclose(fp);

  for(int i = 0; i < n; i++)
    std::sort(this->listaAdj[i].begin(), this->listaAdj[i].end(), compareAresta);

  this->particao.resize(n);

  /*
  for(int i = 0; i < n / 2; i++)
    this->particao[i] = false;
  for(int i = n / 2; i < n; i++)
    this->particao[i] = true;

  //Mistura a solucao inicial
  for(int i = 0; i < 10 * n; i++){
    int v1 = rand() % n;
    int v2 = rand() % n;
    while(this->particao[v1] == this->particao[v2])
      v2 = rand() % n;
    this->particao[v1] = !this->particao[v1];
    this->particao[v2] = !this->particao[v2];
  }
  */

  fp = fopen(partitionFile, "r");

  if(!fp){
    fprintf(stderr, "Error: Could not open file %s\n", partitionFile);
    exit(1);
  }

  int rb = 0;
  int readV;
  for(int i = 0; i < n; i++){
    rb += fscanf(fp, "%d\n", &readV);
    this->particao[i] = ((readV == 0) ? false : true);
  }

  fclose(fp);

  this->valorP = this->calcValorParticao();

  this->tempoMicro = 0;
  this->iteracoes = 0;
}

int Grafo::getValorParticao(){
  return this->valorP;
}

std::vector<bool> Grafo::getParticao(){
  return this->particao;
}

long long int Grafo::getTempoMicro(){
  return this->tempoMicro;
}

long long int Grafo::getIteracoes(){
  return this->iteracoes;
}

void Grafo::tabu(const bool verbose, const int dTabu, const bool randD){
  const int n = (int) this->listaAdj.size();

  auto tempoInicial = std::chrono::steady_clock::now();

  Fila fila(&this->listaAdj, &this->particao, dTabu, randD);
  std::pair<int, int> vizinho;

  int lastImprov = 0;

  std::vector<bool> incumbente = this->particao;
  int valorIncumbente = this->calcValorParticao();

  while(lastImprov < (n / 4)){
    this->iteracoes++;

    vizinho = fila.melhorMelhora();
    if(vizinho.first != -1){
      this->valorP = this->valorP + fila.swap(vizinho.first, vizinho.second);
      if(valorIncumbente > this->valorP){
        valorIncumbente = this->valorP;
        incumbente = this->particao;
        lastImprov = 0;
      }
      else{
        lastImprov++;
      }
    }
    else{
      lastImprov++;

      int v1 = rand() % n;
      while(fila.isTabu(v1))
        v1 = rand() % n;

      int v2 = rand() % n;
      while(fila.isTabu(v2) || this->particao[v1] == this->particao[v2])
        v2 = rand() % n;

      this->valorP = this->valorP + fila.swap(v1, v2);
    }

    if(verbose)
      printf("Valor: %d\n", this->valorP);
  }

  this->particao = incumbente;
  this->valorP = this->calcValorParticao();

  auto tempoFinal = std::chrono::steady_clock::now();
  long long int tempoTotal = (long long int) std::chrono::duration_cast<std::chrono::microseconds>(tempoFinal - tempoInicial).count();
  this->tempoMicro = this->tempoMicro + tempoTotal;
}

void Grafo::randomizada(const bool verbose, const float p, const int timeoutS){
  const int n = (int) this->listaAdj.size();

  auto tempoInicial = std::chrono::steady_clock::now();

  Fila fila(&this->listaAdj, &this->particao, 0, false);
  std::pair<int, int> vizinho;

  std::vector<bool> incumbente = this->particao;
  int valorIncumbente = this->calcValorParticao();

  bool stop = false;
  long long int timeout = ((long long int) timeoutS) * 1000000;
  while(!stop){
    this->iteracoes++;

    if((((float) rand()) / (RAND_MAX)) < p){
      int v1 = rand() % n;

      int v2 = rand() % n;
      while(this->particao[v1] == this->particao[v2])
        v2 = rand() % n;

      this->valorP = this->valorP + fila.swap(v1, v2);
    }
    else{
      vizinho = fila.melhorMelhora();
      if(vizinho.first != -1){
        this->valorP = this->valorP + fila.swap(vizinho.first, vizinho.second);
      }
      else{
        if(p < 0.00001)
          stop = true;
      }
    }

    if(valorIncumbente > this->valorP){
      valorIncumbente = this->valorP;
      incumbente = this->particao;
    }

    if(verbose)
      printf("Valor: %d\n", this->valorP);

    auto tempoAtual = std::chrono::steady_clock::now();
    long long int tempoDecorrido = (long long int) std::chrono::duration_cast<std::chrono::microseconds>(tempoAtual - tempoInicial).count();
    if(tempoDecorrido > timeout)
      stop = true;
  }

  this->particao = incumbente;
  this->valorP = this->calcValorParticao();

  auto tempoFinal = std::chrono::steady_clock::now();
  long long int tempoTotal = (long long int) std::chrono::duration_cast<std::chrono::microseconds>(tempoFinal - tempoInicial).count();
  this->tempoMicro = this->tempoMicro + tempoTotal;
}

int Grafo::calcValorParticao(){
  int valor = 0;

  for(int i = 0; i < (int) this->listaAdj.size(); i++){
    if(this->particao[i] == true){
      for(int j = 0; j < (int) this->listaAdj[i].size(); j++){
        if(this->particao[this->listaAdj[i][j].destino] == false){
          valor++;
        }
      }
    }
  }

  return valor;
}
