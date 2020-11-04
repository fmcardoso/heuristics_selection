//
// Created by fmcardoso on 9/15/17.
//
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>

#ifndef COOPERATION_CONTROL_GRAPH_H
#define COOPERATION_CONTROL_GRAPH_H

struct Edge {
    int u;
    int v;
};

struct Neighbors {
    const long n;
    const int *neighbors;
};

class Graph {
    int size;
    std::vector<std::vector<int>> adj_list;
    std::vector<Edge> edges;
    bool regenerate_edges;
protected:
    bool directed;
public:
    Graph (int, bool);
    int add_node();
    std::vector<Edge> get_edges();
    void remove_edge(int, int);
    void add_edge(int, int);

    virtual std::vector<int> get_neighbors(int);
    virtual void wmUpdate() {

    };

    bool isNeighbor(int u, int v);

    std::vector<int> get_degreeDist();

    int getSize();

    void clearEdges();

    static const double DBL_MAX;
    static const int INT_MAX_M;

    bool checkConnectedness();
};


class RRN: public Graph{
public:
    RRN(int, int);
};

class Lattice: public Graph{
public:
    Lattice(int L);
};

Graph readEdgeList(const char *file_path, bool directed, bool one_indexed);


#endif //COOPERATION_CONTROL_GRAPH_H