#ifndef TREE_ViLa_H
#define TREE_ViLa_H

#endif //TREE_ViLa_H

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "../graph/graph.h"
#include <utility/utils.h>
using boost::dynamic_bitset;

//private List<Map<Integer, Set<Node>>> ViLa;

struct ViLa_Edge {
    int nei_id;
    dynamic_bitset<> lifespan;
    ViLa_Edge(){};
    ViLa_Edge(int nei_) {
        nei_id = nei_;
        lifespan = {};
    }
    bool operator<(const uint r) const {
        return nei_id < r;
    }
};

struct ViLa_Node {
    int id;
    std::map<int, dynamic_bitset<>> labels;
    ViLa_Node(){ id = -1; };
    ViLa_Node(int id_) {
        id = id_;
        labels = {};
    }
};

class ViLa_Graph {
public:
    int T;
    std::vector<ViLa_Node> nodes;
    uint node_nums;
    std::vector<std::vector<ViLa_Edge>> ViLaEdge;
    // std::vector<std::map<int, std::set<ViLa_Node>>> ViLa;

    ViLa_Graph() {
        node_nums = 0;
        nodes = {};
        ViLaEdge = {};
        // ViLa = {};
    }

    void build(int K, std::string query_path);
    void addEdge(int a, int b, int time);
    void addEdge(int a, int b, int la, int lb, int time);

    void query(std::string input_file, bool fun);
};

