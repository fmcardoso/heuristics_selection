//
// Created by fmcardoso on 9/15/17.
//
#include <memory>
#include "Graph.h"
#include <iostream>
#include <cassert>
#include <deque>
#include "MyRandom.h"
#include <cmath>
#include <algorithm>
#include <numeric>

const double Graph::DBL_MAX = std::numeric_limits<double>::max() / 10; // TODO - find a better way to get a MAX
const int Graph::INT_MAX_M = std::numeric_limits<int>::max() / 10; // The division is used to avoid overflow when summing

Graph::Graph(const int size, const bool directed) {
    this->regenerate_edges = false;
    this->size = size;
    this->adj_list.resize(size);
    this->directed = directed;
}

/**
 * Return the edges of the graph. If the regenerate_edges flag is set, it will
 * regenerate the edges of the graph in O(E).
 * @return
 */
std::vector<Edge> Graph::get_edges() {
    if (regenerate_edges) {
        this->edges.clear();

        for (int u = 0; u < size; u++) {
            for (auto &v : this->adj_list[u]) {

                if (u < 0 or 0 > v or u >= size or v >= size) {
                    printf("Problem -----%d %d", u, v);
                }

                if (u < v) {
                    Edge e = {u, v};
                    this->edges.push_back(e); // To assure that only add once
                }
            }
        }
        regenerate_edges = false;
    }

    return this->edges;
}

std::vector<int> Graph::get_degreeDist(){
    std::vector<int> degrees(this->size, 0);

    for (int u = 0; u < size; u++) {
        degrees[u] = static_cast<int>(adj_list[u].size());
    }

    return degrees;
}

void Graph::clearEdges() {
    this->edges.clear();

    for (int u = 0; u < size; u++) {
        adj_list[u].clear();
    }
}

void Graph::remove_edge(int u, int v) {
    std::vector<int> ulist = this->adj_list[u];
    std::vector<int> vlist = this->adj_list[v];

//    std::cout <<  this->adj_list[u].size() << this->adj_list[v].size() << " - ";

    ulist.erase(std::remove(ulist.begin(), ulist.end(), v), ulist.end());
    this->adj_list[u] = ulist;

    if (!this->directed) {
        vlist.erase(std::remove(vlist.begin(), vlist.end(), u), vlist.end());
        this->adj_list[v] = vlist;
    }

//    std::cout <<  this->adj_list[u].size() << this->adj_list[v].size() << std::endl;

//    this->edges.clear();
    this->regenerate_edges = true; // Don'T change, some functions assume that the edges are not updated upon removal
}

void Graph::add_edge(int u, int v) {
    if (u >= this->size or v >= this->size) {
        throw "Node id >= size!";
    }

    this->adj_list[u].push_back(v);

    if (!this->directed) {
        this->adj_list[v].push_back(u);
        if (u < v) {
            Edge e = {u, v};
            this->edges.push_back(e);
        } else {
            Edge e = {v, u};
            this->edges.push_back(e);
        }
    } else {
        Edge e = {u, v};
        this->edges.push_back(e);
    }

}

std::vector<int> Graph::get_neighbors(int u) {
//    Neighbors ngb = {adj_list[u].size(), adj_list[u].data()};
    return adj_list[u];
}

bool Graph::isNeighbor(int u, int v) {
    return std::find(this->adj_list[u].begin(), this->adj_list[u].end(), v) != this->adj_list[u].end();
//    for (auto v_ : this->adj_list[u]){
//        if (v == v_){
//            return true;
//        }
//    }
//    return false;
}

int Graph::add_node() {
    int u = this->size;
    size++;
    std::vector<int> adj_list;
    this->adj_list.push_back(adj_list);

    return u;
}

int Graph::getSize() {
    return size;
}

/**
 * Construct a Random Regular Graph, it executes until it obtains a simple graph
 * @param n Size of the graph, i.e., number of nodes.
 * @param k average degree which is fixed for all nodes.
 */
RRN::RRN(int n, int k) : Graph(n, false) {
    std::vector<int> nodes(((unsigned long) n * k));

    MyRandom rng(rand());

    int j = 0;
    int u, v;

    for (u = 0; u < n; u++) {
        for (int i = 0; i < k; i++) {
            nodes[j] = u;
            j++;
        }
    }

    // Executes until a simple graph is obtained
    while (true) {
        this->clearEdges();
        RandomHelper::shuffle(&nodes, &rng);

        for (int i = 0; i < (n * k) / 2; i++) {
            u = nodes[i * 2];
            v = nodes[(i * 2) + 1];

            if (u == v or this->isNeighbor(u, v)) {
                break;
            }

            this->add_edge(u, v);
        }

        if (get_edges().size() == n * k / 2) {
            return;
        }
    }
}

/**
 * Consider an unweighted edge list
 * @param file_path
 * @param directed
 * @param one_indexed - If the nodes in the file starts in 1.
 * @return
 */
Graph readEdgeList(const char *file_path, bool directed, bool one_indexed) {
    std::ifstream infile(file_path);
    Graph g(0, directed);

    assert(infile.good());

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) {
            std::cout << "Error reading file!" << std::endl;
            break;
        } // error

        if (one_indexed) {
            u--;
            v--;
        }
        // Add nodes if they are not included in graph
        while (g.getSize() < u + 1 || g.getSize() < v + 1) {
            g.add_node();
        }

        g.add_edge(u, v);
    }

    return g;
}

bool Graph::checkConnectedness() {
    int u = 0;
    std::vector<bool> visited(this->size, false);
    std::deque<int> toVisit;
    toVisit.push_back(u);
    visited[u] = true;
    int nVisited = 1;

    while (not toVisit.empty()){
        u = toVisit.front();
        toVisit.pop_front();

        auto ngb = this->get_neighbors(u);

        for (auto v:ngb){
            if (not visited[v]){
                visited[v] = true;
                nVisited++;
                toVisit.push_back(v);
            }
        }
    }

//    std::cout << nVisited <<  " - " << std::accumulate(visited.begin(), visited.end(), 0) << std::endl;

    assert(nVisited == std::accumulate(visited.begin(), visited.end(), 0));

    return nVisited == this->getSize();
}

/**
 * Build a square lattice with periodic boundary conditions.
 * @param L dimension
 */
Lattice::Lattice(int L) : Graph(L * L, false) {
    int N = L * L;
    int u;

    for (u = 0; u < N; u++) {
        if ((u + 1) % L != 0){
            this->add_edge(u, u + 1);
        }else{
            this->add_edge(u, u + 1 - L);
        }

        if (u + L < N){
            this->add_edge(u, u + L);
        }
        else{
            this->add_edge(u, u % L);
        }
    }

}
