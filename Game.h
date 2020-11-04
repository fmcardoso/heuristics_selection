//
// Created by fmcardoso on 13/06/19.
//

#ifndef COOPGA_GAME_H
#define COOPGA_GAME_H

#include <vector>
#include <deque>
#include "Graph.h"

#define Cooperate true
#define Defect false
//#define nFParams 3

// Row player Payoff matrix - First letter = row player action, Second letter - column player action
// DO NOT USE NEGATIVE PAYOFFS
// TODO - normalize inside the game
#define PAYOFF_T 5.0
#define PAYOFF_R 3.0
#define PAYOFF_P 1.0
#define PAYOFF_S 0.0
//#define PAYOFF_R 5.0
//#define PAYOFF_S 0.0
//#define DC 3.0
//#define DD 4.0

//typedef struct individualF {double f[nFParams];} individualF;

class individualF{
public:
//    double f[nFParams];
    std::vector<int8_t> f;
};

typedef struct generationOutcome {double perCoop; double totalF;} generationOutcome;

class GenerationOutcome {
public:
    GenerationOutcome(const std::vector<double> &payoffs, const std::vector<int> &reputation,
                          double totalF, double perCoop);
    GenerationOutcome() = default;

    std::vector<double> payoffs;
    std::vector<int> reputation;

    double perCoop;
    double totalF;
};

class Game {
    int N, T, G, lags, nGenes, L;
    float pMutate;
    std::vector<std::vector<std::vector<int8_t >>> actions;

    float payoff_T, payoff_R, payoff_P, payoff_S;
    int nEdges;
    std::vector<int> degreeDist;
    std::string geneConf;

    GenerationOutcome performInteractions(const std::vector<individualF> &populationF, int lags);

public:
    Game(const int N, const int T, const int G, Graph *graph, int lags, const int L, const float pMutate,
         const std::string geneConf, const float payoffT);
    Graph* graph;
    void play(int nExec, bool saveMatrix, const std::string& agentsPathPrefix);
    std::vector<individualF> populationF;

    void buildX(std::vector<double> &xu, const std::vector<int> &rU, const std::vector<double> &pU,
                const std::vector<int> &rV, const std::vector<double> &pV, const std::vector<int8_t> &uC,
                const std::vector<int8_t> &vC, int vDegree, int u, int v);

    double getJaccardIndex(int u, int v);
};


#endif //COOPGA_GAME_H
