#ifndef GRAPH_H
#define GRAPH_H

#include "type.h"

#include <string>
#include <unordered_map>
#include <vector>
class GraphDelta {
private:
    ui vertices_count_; // 点个数
    ui first_vertices_count_;
    ui last_vertices_count_;
    ui edges_count_; // 边个数
    ui sumT_ = 0;
    std::string dataid_;

    vector<vector<ui>> neighbor; // 每个点邻居偏移量 大小vertices_count_+1

    vector<vector<ui>> last_neighbor;

    vector<vector<pair<ui, ui>>> addEdge_first;
    vector<vector<pair<ui, ui>>> deleteEdge_first;

    vector<vector<pair<ui, ui>>> addEdge_last;
    vector<vector<pair<ui, ui>>> deleteEdge_last;

public:
    GraphDelta() {
        vertices_count_ = 0;
        edges_count_ = 0;
        sumT_ = 0;
    }
    ~GraphDelta() {}

    void lordDateset(const std::string &file_path);
    void getGraphL(ui L, vector<pair<ui, ui>> &ans);
    vector<pair<ui, ui>> queryAnd(ui l, ui r);
    vector<pair<ui, ui>> queryOr(ui l, ui r);
    double getMemoryMB();

    const ui getVerticesCount() const {
        return vertices_count_;
    }
    const ui getEdgesCount() const {
        return edges_count_;
    }
};

#endif