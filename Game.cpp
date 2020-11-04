//
// Created by fmcardoso on 13/06/19.
//

#include "Game.h"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>
#include <math.h>
#include <map>
#include "MyRandom.h"
#include <bitset>

#define K_SIGMOID 0.05
#define GENE_SIZE 8

static_assert(-1 == 0xffffffff, "not two's complement");
static_assert(127 == 0b01111111, "not two's complement");

static const std::string Cuv = "Cuv";
static const std::string Ruv = "Ruv";
static const std::string Puv = "Puv";
static const std::string NCuv = "NCuv";
static const std::string All = "All";
static const std::string Kin = "Kin"; // Add kin selection

Game::Game(const int N, const int T, const int G, Graph *graph, int lags, const int L, const float pMutate,
           const std::string geneConf, const float payoffT)
        : graph(graph) {
    this->N = N;
    this->T = T;
    this->G = G;
    this->lags = lags;
    populationF.resize(N);
    this->geneConf = geneConf;
    if (geneConf == All){
        this->nGenes = 1 + 4*lags;
    } else if (geneConf == Kin){
        this->nGenes = 2 + 4*lags;
    }
    else{
        this->nGenes = 1 + 3*lags;
    }

    this->L = L;
    this->pMutate = pMutate;

    degreeDist = graph->get_degreeDist();

    this->payoff_T = payoffT / payoffT;
    this->payoff_R = static_cast<float>(PAYOFF_R) / payoffT;
    this->payoff_P = static_cast<float>(PAYOFF_P) / payoffT;
    this->payoff_S = static_cast<float>(PAYOFF_S) / payoffT;

    this->nEdges = graph->get_edges().size();

    this->actions.resize(N);
    for (int i=0; i < N; i++){
        actions[i].resize(N);
        for (int j=0; j < N; j++){
            actions[i][j].resize(lags, 0);
        }
    }
}

void Game::play(int nExec, bool saveMatrix, const std::string& agentsPathPrefix) {
    GenerationOutcome outcome;

    MyRandom rng(rand());


    uniform_int8_distribution initialDistribution(-100, 100);
    uniform_int_distribution bitFlipDist(0, GENE_SIZE -1);
    uniform_real_distribution pMutateDist(0.0, 1.0);
    uniform_int_distribution nodesDist(0, N-1);
    uniform_int_distribution geneDist(0, nGenes - 1);

    static_assert(GENE_SIZE == 8, "GENE ERROR");

    std::vector<int> nodes_vector(N, 0); // Used for the death sample procedure

    for (int u =0; u < N; u++){
        nodes_vector[u] = u;
    }

    // Initialize first generation
    for (int u = 0; u < N; u++){
        this->populationF[u].f.resize(nGenes);
        for (int8_t & i : this->populationF[u].f){
            i = initialDistribution.Random(&rng);
        }
    }

    // Write headers
    std::cout << "Generation\tPerC\tTotalPayoff\tBeta0\tStdBeta0";

    if (geneConf == Kin){
        std::cout << "\tKin\tStdKin";
    }

    int nVars = 4;
    if ((geneConf != All) and (geneConf != Kin)){
        nVars = 3;
    }

    std::vector<std::string> labels(nVars);
    int lastVar = 0;

    if (this->geneConf != Cuv){
       labels[lastVar] = "Cother";
       lastVar++;
    }
    if (this->geneConf != Ruv){
        labels[lastVar] = "Rother";
        lastVar++;
    }
    if (this->geneConf != Puv){
        labels[lastVar] = "Pother";
        lastVar++;
    }
    if (this->geneConf != NCuv){
        labels[lastVar] = "NotCother";
    }


    for (const auto& lab : labels){
        for (int l = 0; l < lags; l++) {
            std::cout << "\t" << "Mean" << lab << l + 1;
            std::cout << "\t" << "Std" << lab << l + 1;
        }
    }
    std::cout << std::endl;

    for (int g = 0; g < G; g++) {
        outcome = performInteractions(this->populationF, this->lags);

        // Write current frac coop and total fitness
        std::cout << g << "\t" << outcome.perCoop << "\t" << outcome.totalF;

        // Write Mean and standard deviations of the distributions

        for (int var = 0; var < this->nGenes; var++) {
            float meanVar = 0.0;
            float stdVar = 0.0;

            for (int u = 0; u < N; u++) {
                meanVar += this->populationF[u].f.at(var);

            }
            meanVar = meanVar / (float) N;
            for (int u = 0; u < N; u++) {
                stdVar += pow(meanVar - this->populationF[u].f.at(var), 2);
            }
            stdVar = sqrt(stdVar / (float) (N - 1));
            std::cout << "\t" << meanVar << "\t" << stdVar;
        }
        std::cout << std::endl;

        if (saveMatrix and ((g % 10000 == 0) or g > (G-100))) {
            std::ofstream outfile;
            std::stringstream filePath;
            filePath << agentsPathPrefix << g << ".tsv";

            outfile.open(filePath.str().c_str(), std::ios_base::trunc);

            // Write headers
            outfile << "u\tBeta0";
            if (geneConf == Kin){
                outfile << "\tKin";
            }
            for (const auto &lab : labels) {
                for (int l = 0; l < lags; l++) {
                    outfile << "\t" << lab << l + 1;
                }
            }
            // Write last reputation and payoff headers
            outfile << "\tReputation" << "\tPayoff";
            outfile << std::endl;

            for (int u = 0; u < N; u++) {
                outfile << u;
                for (auto cU : this->populationF[u].f) {
                    outfile << "\t" << (int) cU;
                }
                // Write last reputation and payoff
                outfile << "\t" << outcome.reputation[u] << "\t" << outcome.payoffs[u];

                outfile << std::endl;
            }
        }

        // Save the old population vector for reproduction
        std::vector<individualF> currentPopulationF = this->populationF;

        for (int d = 0; d < N; d++){
//        for (auto d: deads){
            int chosen = -1;

            std::vector<int> ngb = graph->get_neighbors(d);
            ngb.push_back(d); // Add the node to the N2 neighborhood
            unsigned long nNgb = ngb.size();
            std::vector<double> pReproduce(ngb.size(), 0);

            // Get the fitness of each node
            double ngbTotalF = 0;
            for (unsigned int i = 0; i < nNgb; i++) {
                int u = ngb[i];
                pReproduce[i] = outcome.payoffs[u];

                ngbTotalF += pReproduce[i];
                if (i > 0) {
                    pReproduce[i] = pReproduce[i] + pReproduce[i - 1];
                }
            }

            // Choose one proportionally to its fitness
            uniform_real_distribution nodeSurvivalDistribution(0, ngbTotalF);
            double k = nodeSurvivalDistribution.Random(&rng);
            for (unsigned int i = 0; i < nNgb; i++) {
                if (pReproduce[i] > k) {
                    chosen = i;
                    break;
                }
            }
            chosen = ngb[chosen]; // Real index of neighbor

            // Change the behavior for the most suitable neighbor
            this->populationF[d] = currentPopulationF[chosen];
        }


        // Add mutations to pMutate*N
        for (int a = 0; a < N; a++) {
            if (pMutateDist.Random(&rng) <= pMutate) {
                int gene = geneDist.Random(&rng);
                individualF *u = &populationF[a];
                int bf = bitFlipDist.Random(&rng);
//            for (auto & j : u->f){
//                j += mutationDistribution(rng);
//            }
                u->f[gene] ^= 1UL << bf;
            }
        }
    }
}

double activationFunction(const std::vector<int8_t> &fu, const std::vector<double> &xu, const int nGenes) {
    double x = 0;

//    assert(fu.size() == xu.size()); // TODO - remove
    for (int i =0; i < nGenes; i++){
        x += fu[i] * xu[i];
    }

    double pCoop = 1.0 / (1.0 + exp(-K_SIGMOID * x));

    return pCoop;
}


void Game::buildX(std::vector<double> &xu, const std::vector<int> &rU, const std::vector<double> &pU,
                  const std::vector<int> &rV, const std::vector<double> &pV, const std::vector<int8_t> &uC,
                  const std::vector<int8_t> &vC, int vDegree, int u, int v) {

    xu[0] = 1.0;

    if (this->geneConf == Kin){
        xu[1] = getJaccardIndex(u, v);
    }

    int lastVar = 1;
    for (int l = 0; l < this->lags; l++){
        lastVar = 1 + (this->geneConf == Kin);
//        // Opposite player values

        if (this->geneConf != Cuv){
            xu[l + lastVar] = (vC[l] > 0);
            lastVar += this->lags;
        }
        if (this->geneConf != Ruv){
            int vCoop = vC[l] > 0;

            if (rV[l] - vCoop){
                assert(((rV[l] - vCoop) / (double) (vDegree - 1)) <= 1);
                xu[l + lastVar] = (rV[l] - vCoop) / (double) (vDegree - 1); // Remove the reputation with edge player to not count it twice
            }
            else{
                xu[l + lastVar] = 0;
            }

            lastVar += this->lags;
        }
        if (this->geneConf != Puv){
            xu[l + lastVar] = pV[l];
            lastVar += this->lags;
        }
        if (this->geneConf != NCuv){
            xu[l + lastVar] = (vC[l] < 0);
            lastVar += this->lags;
        }
         assert(lastVar == this->nGenes);

//        xu[l + 1 + 4*lags] = (rV[l] - vC) == 0; // Only applies agents who did not cooperated with anybody else
    }
//    static_assert(1 + 3 == N_GENES - 1, "Wrong lag numbers");
}

double Game::getJaccardIndex(int u, int v){
    double index = 0;

    // Compute how many bits differ
    for (int i = 0; i < this->nGenes; i++) {
        size_t diff = this->populationF[u].f[i] ^ this->populationF[v].f[i];
        std::bitset<8> diff_bits(diff);
        index += diff_bits.count();
    }

    index = 1 - (index / (nGenes * GENE_SIZE));

    return index;
}


/**
 * This function assume that the first interaction population is already initialized
 * @return
 */
GenerationOutcome Game::performInteractions(const std::vector<individualF> &populationF, int lags) {
    MyRandom rng(rand());
    uniform_real_distribution distribution(0.0, 1.0);

    std::vector<double> accPayoff (N, 0);

    std::vector<std::vector<int>> reps(N, std::vector<int>(lags, 0));
    std::vector<std::vector<double>> pffs(N, std::vector<double>(lags, 0.0));

    // TODO CHECK THIS
    // For each player and action
    std::vector<double> xu(this->nGenes, 1.0);
    std::vector<double> xv(this->nGenes, 1.0);

    std::vector<double> payoffs(static_cast<unsigned long>(this->N), 0.0);
    std::vector<int> reputation(static_cast<unsigned long>(this->N), 0.0);

    std::vector<bool> temp(lags, false);

    // Perform
    for (int t =0; t < T; t++){
        this->graph->wmUpdate(); // Update the network in the case of the well mixed
        std::vector<double> tPayoffs(static_cast<unsigned long>(this->N), 0.0);
        std::vector<int> tReputation(static_cast<unsigned long>(this->N), 0);

        // Perform pairwise play
        for (auto e: this->graph->get_edges()){
            bool uC;
            bool vC;

            // Nobody has previous actions in the beginning
            if (t==0){
                for (int l = 0; l < lags; l++) {
                    actions[e.v][e.u][l] = 0;
                    actions[e.u][e.v][l] = 0;
                }
            }

            buildX(xu, reps[e.u], pffs[e.u], reps[e.v], pffs[e.v], actions[e.u][e.v], actions[e.v][e.u],
                   this->degreeDist[e.v], e.u, e.v);
            buildX(xv, reps[e.v], pffs[e.v], reps[e.u], pffs[e.u], actions[e.v][e.u], actions[e.u][e.v],
                   this->degreeDist[e.u], e.v, e.u);

            uC = activationFunction(populationF[e.u].f, xu, this->nGenes) > distribution.Random(&rng);
            vC = activationFunction(populationF[e.v].f, xv, this->nGenes) > distribution.Random(&rng);


            // Update actions
            for (int l = lags -1; l > 0; l--){
                actions[e.u][e.v][l] = actions[e.u][e.v][l-1];
                actions[e.v][e.u][l] = actions[e.v][e.u][l-1];;
            }

            // Update local payoffs
            if (uC and vC){
                tPayoffs[e.u] += this->payoff_R / (float) this->degreeDist[e.u];
                tPayoffs[e.v] += this->payoff_R / (float) this->degreeDist[e.v];
                tReputation[e.u] += 1.0;
                tReputation[e.v] += 1.0;
                if (lags > 0){
                    actions[e.u][e.v][0] = 1;
                    actions[e.v][e.u][0] = 1;
                }
            } else if(uC){
                tPayoffs[e.u] += this->payoff_S / (float) this->degreeDist[e.u];
                tPayoffs[e.v] += this->payoff_T / (float) this->degreeDist[e.v];
                tReputation[e.u] += 1.0;
                if (lags > 0) {
                    actions[e.u][e.v][0] = 1;
                    actions[e.v][e.u][0] = -1;
                }
            } else if(vC){
                tPayoffs[e.u] += this->payoff_T  / (float) this->degreeDist[e.v];
                tPayoffs[e.v] += this->payoff_S / (float) this->degreeDist[e.v];
                tReputation[e.v] += 1.0;
                if (lags > 0) {
                    actions[e.u][e.v][0] = -1;
                    actions[e.v][e.u][0] = 1;
                }
            } else{
                tPayoffs[e.u] += this->payoff_P / (float) this->degreeDist[e.v];
                tPayoffs[e.v] += this->payoff_P / (float) this->degreeDist[e.v];
                if (lags > 0) {
                    actions[e.u][e.v][0] = -1;
                    actions[e.v][e.u][0] = -1;
                }
            }
        }

        // Update payoffs and reputations
       if(lags > 0){
           // TODO - probably this can be more efficient
           for (int u = 0 ; u < N; u++){
               reps[u].pop_back();
               reps[u].insert(reps[u].begin(), tReputation[u]);
               pffs[u].pop_back();
               pffs[u].insert(pffs[u].begin(), tPayoffs[u]);
           }
       }

        payoffs = tPayoffs;
        reputation = tReputation;
        // Save the accumulated payoff
        for (int u = 0 ; u < N; u++){
            accPayoff[u] += tPayoffs[u];
        }
    }

    // Payoff is the accumulated payoff
    payoffs = accPayoff;

    double perCoop = std::accumulate(reputation.begin(), reputation.end(), 0.0)/(2*this->nEdges);
    double totalF = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);

    GenerationOutcome outcome(payoffs, reputation, totalF, perCoop);

    return outcome;
}

GenerationOutcome::GenerationOutcome(const std::vector<double> &payoffs, const std::vector<int> &reputation,
                                     double totalF, double perCoop)
        : payoffs(payoffs), reputation(reputation), perCoop(perCoop), totalF(totalF) {}
