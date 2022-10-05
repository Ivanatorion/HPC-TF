#ifndef FILA_H
#define FILA_H

typedef struct ars{
  int origem;
  int destino;
} Aresta;

bool compareAresta(const Aresta a, const Aresta b);

class Fila{
  public:
    Fila(std::vector<std::vector<Aresta>> *listaAdj, std::vector<bool> *particao, int dTabuMax, bool randD);

    std::pair<int, int> melhorMelhora();
    int swap(int i, int j);

    bool isTabu(const int v);

  private:

    bool randD;

    std::vector<std::vector<Aresta>> *listaAdj;
    std::vector<bool> *particao;
    int deg;

    std::vector<int> dTabu;
    int dTabuMax;

    int edge(const int v1, const int v2);
    std::vector<int> swapDelta;

    std::vector<std::unordered_set<int>> pbucks;
    std::vector<std::unordered_set<int>> qbucks;

    int sDist, sY;
    void getSeq(int *pm, int *qm);
};

#endif
