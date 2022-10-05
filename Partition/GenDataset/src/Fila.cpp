#include <vector>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstdint>
#include <queue>
#include <unordered_set>
#include <algorithm>

#include "../include/Fila.h"

#define INT_MAX 2000000000

bool compareAresta(const Aresta a, const Aresta b) {return (a.destino < b.destino); }

Fila::Fila(std::vector<std::vector<Aresta>> *listaAdj, std::vector<bool> *particao, int dTabuMax, bool randD){
  this->listaAdj = listaAdj;
  this->particao = particao;
  this->dTabuMax = dTabuMax;

  this->randD = randD;

  const int n = (int) (*this->listaAdj).size();
  this->swapDelta.resize(n);
  this->dTabu.resize(n);

  this->deg = 0;
  for(int i = 0; i < n; i++){
    if((int) (*this->listaAdj)[i].size() > this->deg)
      this->deg = (int) (*this->listaAdj)[i].size();
  }

  this->pbucks.resize(this->deg * 2 + 1);
  this->qbucks.resize(this->deg * 2 + 1);

  for(int i = 0; i < n; i++){
    this->swapDelta[i] = 0;
    this->dTabu[i] = 0;
    for(Aresta a : (*this->listaAdj)[i]){
      if((*this->particao)[a.destino] != (*this->particao)[a.origem]){
        this->swapDelta[i]--;
      }
      else{
        this->swapDelta[i]++;
      }
    }

    if((*this->particao)[i] == true)
      this->pbucks[this->swapDelta[i] + this->deg].insert(i);
    else
      this->qbucks[this->swapDelta[i] + this->deg].insert(i);
  }
}

std::pair<int, int> Fila::melhorMelhora(){
  const int n = (int) (*this->listaAdj).size();

  int pm, qm;
  this->sDist = 0;
  this->sY = 0;
  this->getSeq(&pm, &qm);

  int edgv1 = -1;
  int edgv2 = -1;

  while(pm + qm < 2*this->deg){
    for(int v1 : this->pbucks[pm]){
      if(dTabu[v1] == 0){
        for(int v2 : this->qbucks[qm]){
          if(dTabu[v2] == 0){
            if(edge(v1, v2) == 0){
              return std::make_pair(v1, v2);
            }
          }
          else{
            edgv1 = v1;
            edgv2 = v2;
          }
        }
      }
    }

    if(pm + qm + 1 < 2*this->deg){
      for(int v1 : this->pbucks[pm+1]){
        if(dTabu[v1] == 0){
          for(int v2 : this->qbucks[qm]){
            if(dTabu[v2] == 0){
              if(edge(v1, v2) == 0){
                return std::make_pair(v1, v2);
              }
            }
          }
        }
      }

      for(int v1 : this->pbucks[pm]){
        if(dTabu[v1] == 0){
          for(int v2 : this->qbucks[qm+1]){
            if(dTabu[v2] == 0){
              if(edge(v1, v2) == 0){
                return std::make_pair(v1, v2);
              }
            }
          }
        }
      }
    }

    int prev = pm + qm;
    this->getSeq(&pm, &qm);

    if(pm + qm > prev && edgv1 != -1)
      return std::make_pair(edgv1, edgv2);
  }

  return std::make_pair(-1, -1);
}

int Fila::swap(int i, int j){
  const int ij = 2 * this->edge(i,j);

  const int delta = this->swapDelta[i] + this->swapDelta[j] + ij;

  if((*this->particao)[i] != true){
    int aux = i;
    i = j;
    j = aux;
  }

  this->pbucks[this->swapDelta[i] + this->deg].erase(i);
  this->swapDelta[i] = -this->swapDelta[i] - ij;
  this->qbucks[this->swapDelta[i] + this->deg].insert(i);

  this->qbucks[this->swapDelta[j] + this->deg].erase(j);
  this->swapDelta[j] = -this->swapDelta[j] - ij;
  this->pbucks[this->swapDelta[j] + this->deg].insert(j);

  (*this->particao)[i] = false;
  (*this->particao)[j] = true;

  for(Aresta a: (*this->listaAdj)[i]){
    if(a.destino != j){
      if((*this->particao)[a.destino] == true){
        this->pbucks[this->swapDelta[a.destino] + this->deg].erase(a.destino);
        this->swapDelta[a.destino] = this->swapDelta[a.destino] - 2;
        this->pbucks[this->swapDelta[a.destino] + this->deg].insert(a.destino);
      }
      else{
        this->qbucks[this->swapDelta[a.destino] + this->deg].erase(a.destino);
        this->swapDelta[a.destino] = this->swapDelta[a.destino] + 2;
        this->qbucks[this->swapDelta[a.destino] + this->deg].insert(a.destino);
      }
    }
  }

  for(Aresta a: (*this->listaAdj)[j]){
    if(a.destino != i){
      if((*this->particao)[a.destino] == false){
        this->qbucks[this->swapDelta[a.destino] + this->deg].erase(a.destino);
        this->swapDelta[a.destino] = this->swapDelta[a.destino] - 2;
        this->qbucks[this->swapDelta[a.destino] + this->deg].insert(a.destino);
      }
      else{
        this->pbucks[this->swapDelta[a.destino] + this->deg].erase(a.destino);
        this->swapDelta[a.destino] = this->swapDelta[a.destino] + 2;
        this->pbucks[this->swapDelta[a.destino] + this->deg].insert(a.destino);
      }
    }
  }

  for(int v = 0; v < (*this->listaAdj).size(); v++)
    if(this->dTabu[v] > 0)
      this->dTabu[v]--;

  int tabuR = 0;

  if(this->randD)
    tabuR = 1 + rand() % 10;

  this->dTabu[i] = this->dTabuMax + tabuR;
  this->dTabu[j] = this->dTabuMax + tabuR;

  return delta;
}

bool Fila::isTabu(const int v){
  return (this->dTabu[v] > 0);
}

int Fila::edge(const int v1, const int v2){
  Aresta ar;
  ar.origem = v1;
  ar.destino = v2;
  return std::binary_search((*this->listaAdj)[v1].begin(), (*this->listaAdj)[v1].end(), ar, compareAresta) ? 1 : 0;
}

void Fila::getSeq(int *pm, int *qm){
  *pm = this->sY;

  *qm = this->sDist - this->sY;

  if(this->sY  == 0){
    this->sDist = this->sDist + 1;
    this->sY = this->sDist;
  }
  else{
    this->sY--;
  }
}
