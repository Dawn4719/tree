#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "ViLa.h"

void ViLa_Graph::build(int K, std::string query_path) {
    T = K;
    for (int snap_i = 0; snap_i < K; ++snap_i) {
        std::string path =  query_path + std::to_string(snap_i);
        if (!io::file_exists(path.c_str()))
        {
            std::cout << "Failed to open: " << query_path << std::endl;
            exit(-1);
        }
        std::ifstream ifs(path);

        std::string type;
        uint from_id, to_id;
        while (ifs >> from_id >> to_id)
        {
            addEdge(from_id, to_id, snap_i);
            // }
        }
        ifs.close();
    }
}

void ViLa_Graph::addEdge(int a, int b, int time) {
    if (a + 1 > ViLaEdge.size()) ViLaEdge.resize(a + 1);
    if (ViLaEdge[a].empty()) {
        ViLaEdge[a] = {b};
        ViLaEdge[a][0].lifespan.resize(T + 1);
        ViLaEdge[a][0].lifespan.set(time);
    }
    else {
        auto lower = std::lower_bound(ViLaEdge[a].begin(), ViLaEdge[a].end(), b);
        int idx = lower - ViLaEdge[a].begin();
        if (lower != ViLaEdge[a].end() && lower->nei_id == b) {
            ViLaEdge[a][idx].lifespan.set(time);
        }
        else {
            ViLaEdge[a].insert(lower, ViLa_Edge(b));
            ViLaEdge[a][idx].lifespan.resize(T + 1);
            ViLaEdge[a][idx].lifespan.set(time);
        }
    }
}

void ViLa_Graph::addEdge(int a, int b, int la, int lb, int time) {
    if (a + 1 > nodes.size()) nodes.resize(a + 1);
    if (nodes[a].id == -1) {
        node_nums++;
        nodes[a] = ViLa_Node(node_nums);
    }
    if (b + 1 > nodes.size()) nodes.resize(b + 1);
    if (nodes[b].id == -1) {
        node_nums++;
        nodes[b] = ViLa_Node(node_nums);
    }
    if (time >= nodes[a].labels[la].size()) {
        nodes[a].labels[la].resize(time + 1);
    }
    nodes[a].labels[la].set(time);
    if (time >= nodes[b].labels[lb].size()) {
        nodes[b].labels[lb].resize(time + 1);
    }
    nodes[b].labels[lb].set(time);
}

void ViLa_Graph::query(std::string input_file, bool fun){
    std::vector<std::pair<int, int> > res;
    res.reserve(10000);
    std::fstream input(input_file, std::ios::in);
    int L, R;
    while (input >> L >> R) {
        std::cout << "Query: " << L << " " << R;
        res.clear();
        res.reserve(100000);
        bool pas = true;
        dynamic_bitset<uint32_t> bit(R - L + 1);
        for (int mp = 0; mp < ViLaEdge.size(); mp++) {
            for (const auto &ed: ViLaEdge[mp]) {
                std::string s;
                to_string(ed.lifespan, s);
                s = s.substr(s.size() - 1 - R, R - L + 1);
                // s = s.substr(L, R - L + 1);
                if (fun && s.find('1') != std::string::npos) {
                    res.emplace_back(mp, ed.nei_id);
                }
                if (!fun && s.find('0') == std::string::npos) {
                    res.emplace_back(mp, ed.nei_id);
                    pas = false;
                }
            }
        }
        int all_edges = 0;
        for (auto& i : res )
            all_edges++;
        std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;

    }
    input.close();
}
