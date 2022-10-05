#ifndef GRAFO_H
#define GRAFO_H

class Grafo{
  public:
    Grafo(const char *filePath, const char *partitionFile);

    int getValorParticao();
    std::vector<bool> getParticao();

    long long int getTempoMicro();
    long long int getIteracoes();

    void tabu(const bool verbose, const int dTabu, const bool randD);
    void randomizada(const bool verbose, const float p, const int timeoutS);

  private:
    long long int tempoMicro;
    long long int iteracoes;

    std::vector<std::vector<Aresta>> listaAdj;

    //True e False indicam a particao do vertice
    std::vector<bool> particao;

    int valorP;
    int calcValorParticao();
};

#endif
