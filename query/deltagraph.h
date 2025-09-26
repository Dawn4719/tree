#ifndef MATCHING_DELTAGRAPH
#define MATCHING_DELTAGRAPH

#include <map>
#include <vector>

struct EdgeDG {
    uint val;
    std::vector<uint> csr;
    std::vector<std::string> bm;
    EdgeDG()=default;
    EdgeDG(uint a, std::vector<uint> b, std::vector<std::string> c={}): val(a), csr(b), bm(c) {};
    bool operator<(const int r) const {
        return val < r;
    }
};

struct NodeDG {
    uint id;
    uint fa;
    int isMat=0;
    int ednums=0;
    std::map<int, std::vector<EdgeDG>> Neighbors;
    std::map<int, int> Neighbor_delta;
    std::vector<uint> son;
};

class DeltaGraph {
public:
    // Graph cur_graph {};
    std::vector<NodeDG> skeleton;
    int k;
    int leafs;
    int depth;
    int node_nums;
    int pool_bias;
    int pooll, poolr;

    std::string fdif;
    std::vector<std::vector<unsigned long long>> dist;
    std::vector<std::vector<uint>> path;

    std::vector<std::vector<uint>> gpvertex;
    std::vector<uint> p;

    DeltaGraph(int k_, int n);
//    ~DeltaGraph(){};
    void build(std::string path, int n);
    void dif(int l, int r);
    void cha(int big, int small, int store, int big_idx=-1, int store_idx=-1);
    void pool_cha(int big, int small, int store, int big_idx=-1, int store_idx=-1);
    void pool_cup(int big, int small, int store, int big_idx=-1, int store_idx=-1);
    void cha(NodeDG& skeleton, NodeDG& store, int idx);
    void cap(NodeDG&, NodeDG&);
    void graphpool(int level);
    void cup(NodeDG&, NodeDG&, int idx=-1);
    void GetDist();
    std::vector<std::vector<uint>> Prim(int l, int r);
    void skewed();
    void balence();
    void Rskewed();
    void Lskewed();
    void mixed();
    void empty();
    int find(int u);
    void get_snap();

    void query(std::string input_file, int K, bool fun, int level);
};

#endif //MATCHING_DELTAGRAPH
