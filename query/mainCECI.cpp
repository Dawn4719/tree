#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <bit>
#include "../utils/CLI11.hpp"
#include "../utils/globals.h"
#include "../utils/types.h"
#include <boost/dynamic_bitset.hpp>
// #include "../graph/graph.h"
#include "deltagraph.h"
#include "ViLa.h"
#include <fstream>
#include <random>
#include <ctime>
#include <immintrin.h>
#include <omp.h>
#include "../utils/pod.h"
// #include <rocksdb/db.h>
// #include <rocksdb/options.h>
// #include <rocksdb/slice.h>
#include "staticcore.hpp"
#include "reach.hpp"
#include "LLAMA.hpp"
#include <utility/utils.h>
#include "CECI.hpp"
#include "clique.hpp"
// using namespace rocksdb;
using namespace std;
int K = 18;
int level;
bool fun;
size_t EDGE_MAX = 3e4 * 256;

std::string query_path, dataset, input_file;
bool CHECK = true;
size_t all_edges;

bool USEDB;
int depth;

// vector<DB*> bv_db;
// DB* tree_db;

struct TNode {
    TNode() : lson(-1), rson(-1), fa(-1) {
        //        bv.resize(N, false);
        dep = 0;
        id = 0;
        prex = 0;
    };
    dynamic_bitset<uint32_t> bv;
    int id;
    int prex;
    int lson;
    int rson;
    vector<int> sons;
    int fa;
    int dep;

    ~TNode()=default;
};

size_t edge_idx;
std::vector<std::vector<std::pair<uint, uint> > > CSR;
// std::vector<std::vector<uint>> CSR;
struct egg{
    unsigned l, r;
    int id;
    egg() = default;
    egg(unsigned l, unsigned r, int id) : l(l), r(r), id(id) {}
};
// std::vector<egg> id2edge;
std::vector<pair<int, int>> id2edge;

struct Tree {
    std::vector<TNode> tr;
    int leaf;
};

std::vector<Tree> tree;
std::vector<uint> snap_idx;

size_t sum;

auto lb(std::vector<std::pair<uint, uint> > &V, uint val) {
    size_t L = 0, R = V.size();
    while (L < R) {
        size_t MID = (L + R) >> 1;
        if (V[MID].first >= val) R = MID;
        else L = MID + 1;
    }
    return V.begin() + L;
}

// auto lb(std::vector<uint> &V, uint val) {
//     size_t L = 0, R = V.size();
//     while (L < R) {
//         size_t MID = (L + R) >> 1;
//         if (V[MID] >= val) R = MID;
//         else L = MID + 1;
//     }
//     return V.begin() + L;
// }

map<string, size_t> get_index_mem() {
    FILE *fp = fopen("/proc/self/status", "r");
    char line[128];
    map<string, size_t> res;
    while (fgets(line, 128, fp) != NULL) {
        //        if (strncmp(line, "VmPeak", 2) == 0)
        //        {
        //            cout << line << endl;
        ////            printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
        //        }
        if (strncmp(line, "VmRSS:", 6) == 0) {
            string p = line;
            res["now"] = size_t(stoull(p.substr(6)));
            cout << line;
        }
        if (strncmp(line, "VmPeak:", 7) == 0) {
            string p = line;
            res["pk"] = size_t(stoull(p.substr(7)));
            cout << line;
        }
    }
    fclose(fp);
    return res;
}
bool useRoaring;
void build(Tree &tre, bool f) {
    auto &tr = tre.tr;

    size_t nums = tr.size();

    //cout << "depth:" << depth << "nums:" << nums << endl;
    int cur_depth = 1;
    int n = EDGE_MAX;

    for (size_t i = 0; i < nums; ++i) tr[i].id = i, tr[i].dep = 0;
    tr[nums - 1].fa = nums - 1;
    int cnt = n;

    while (cur_depth < depth) {
        // std::cout << "cur_depth_node_num:" << cnt << std::endl;
        int i = n - cnt + 2;
        cnt = 0;
        for (; i <= n; i += 2) {

            tr[n + cnt].lson = i - 2;
            tr[n + cnt].rson = i - 1;
            tr[n + cnt].dep = cur_depth;
            tr[i - 2].fa = n + cnt;
            tr[i - 1].fa = n + cnt;

            if (tr[i - 2].bv.empty() || tr[i - 1].bv.empty()) {
                cnt++;
                continue;
            }

            if (USEDB) {
                if (fun) {
                    if (tr[i - 2].bv.size() <= tr[i - 1].bv.size()) {
                        tr[n + cnt].bv = tr[i - 2].bv;
                        tr[n + cnt].bv.resize(tr[i - 1].bv.size());
                        tr[n + cnt].bv |= tr[i - 1].bv;
                    }
                    else {
                        tr[n + cnt].bv = tr[i - 1].bv;
                        tr[n + cnt].bv.resize(tr[i - 2].bv.size());
                        tr[n + cnt].bv |= tr[i - 2].bv;
                    }
                }
                else {
                    if (tr[i - 2].bv.size() <= tr[i - 1].bv.size()) {
                        tr[n + cnt].bv = tr[i - 1].bv;
                        tr[n + cnt].bv.resize(tr[i - 2].bv.size());
                        tr[n + cnt].bv &= tr[i - 2].bv;
                    }
                    else {
                        tr[n + cnt].bv = tr[i - 2].bv;
                        tr[n + cnt].bv.resize(tr[i - 1].bv.size());
                        tr[n + cnt].bv &= tr[i - 1].bv;
                    }
                }
                /*if (fun) {
                    if (tr[i - 1].prex * 32 + tr[i - 1].bv.size() >= tr[i - 2].prex * 32 + tr[i - 2].bv.size()) {
                        tr[n + cnt].bv = tr[i - 1].bv;
                        tr[n + cnt].prex = min(tr[i - 1].prex, tr[i - 2].prex);
                        int commonSize = tr[i - 1].prex - tr[i - 2].prex;
                        if (commonSize > 0)
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() + 32 * commonSize);

                        for (int l = tr[i - 2].bv.num_blocks() - 1, j = tr[n + cnt].bv.num_blocks() - 1 + min(0, commonSize); ~l && ~j; l--, j--)
                            tr[n + cnt].bv.m_bits[j] |= tr[i - 2].bv.m_bits[l];

                        assert(tr[n + cnt].bv.m_bits.back() != 0);
                    }
                    else{
                        tr[n + cnt].bv = tr[i - 2].bv;
                        tr[n + cnt].prex = min(tr[i - 1].prex, tr[i - 2].prex);
                        int commonSize = tr[i - 1].prex - tr[i - 2].prex;
                        if (commonSize > 0)
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() + 32 * commonSize);

                        for (int l = tr[i - 1].bv.num_blocks() - 1, j = tr[n + cnt].bv.num_blocks() - 1 + min(0, commonSize); ~l && ~j; l--, j--)
                            tr[n + cnt].bv.m_bits[j] |= tr[i - 1].bv.m_bits[l];

                        assert(tr[n + cnt].bv.m_bits.back() != 0);
                    }
                }
                else {
                    if (tr[i - 1].prex * 32 + tr[i - 1].bv.size() <= tr[i - 2].prex * 32 + tr[i - 2].bv.size()) {
                        tr[n + cnt].bv = tr[i - 1].bv;
                        tr[n + cnt].prex = min(tr[i - 1].prex, tr[i - 2].prex);
                        int commonSize = tr[i - 2].prex - tr[i - 1].prex;
                        int l = tr[i - 2].bv.num_blocks() - 1 + min(commonSize, 0), j = tr[n + cnt].bv.num_blocks() - 1;
                        tr[n + cnt].bv.resize(tr[n + cnt].bv.size() - min(0, commonSize) * 32);
                        int zeroCnt = 0;
                        commonSize = max(0, commonSize);
                        for (; commonSize > 0; j--, commonSize--) {
                            tr[n + cnt].bv.m_bits[j] = 0;
                            zeroCnt++;
                        }

                        for (; ~l && ~j; l--, j--)
                            tr[n + cnt].bv.m_bits[j] &= tr[i - 2].bv.m_bits[l];
                        for (j = tr[n + cnt].bv.num_blocks() - 1; ~j; j--) {
                            if (tr[n + cnt].bv.m_bits[j] != 0) {
                                tr[n + cnt].prex += tr[n + cnt].bv.num_blocks() - j - 1;
                                tr[n + cnt].bv.resize((j + 1) * 32);
                                tr[n + cnt].bv.shrink_to_fit();
                                break;
                            }
                        }
                        // assert(tr[n + cnt].bv.m_bits.back() != 0);
                    }
                    else {
                        tr[n + cnt].bv = tr[i - 2].bv;
                        tr[n + cnt].prex = min(tr[i - 1].prex, tr[i - 2].prex);
                        int commonSize = tr[i - 1].prex - tr[i - 2].prex;
                        int l = tr[i - 1].bv.num_blocks() - 1 + min(commonSize, 0), j = tr[n + cnt].bv.num_blocks() - 1;
                        tr[n + cnt].bv.resize(tr[n + cnt].bv.size() - min(0, commonSize) * 32);
                        // tr[n + cnt].prex -= min(0, commonSize);
                        int zeroCnt = 0;
                        commonSize = max(0, commonSize);
                        for (; commonSize > 0; j--, commonSize--) {
                            tr[n + cnt].bv.m_bits[j] = 0;
                            zeroCnt++;
                        }
                        for (; ~l && ~j; l--, j--)
                            tr[n + cnt].bv.m_bits[j] &= tr[i - 1].bv.m_bits[l];

                        for (j = tr[n + cnt].bv.num_blocks() - 1; ~j; j--) {
                            if (tr[n + cnt].bv.m_bits[j] != 0) {
                                tr[n + cnt].prex += tr[n + cnt].bv.num_blocks() - j - 1;
                                tr[n + cnt].bv.resize((j + 1) * 32);
                                tr[n + cnt].bv.shrink_to_fit();
                                break;
                            }
                        }
                        // assert(tr[n + cnt].bv.m_bits.back() != 0);
                    }
                }*/
            }
            cnt++;
        }
        /*if (i - 2 != n) {
            if (tr[n + cnt].lson == -1) {
                tr[n + cnt].lson = n - 1;
                tr[n + cnt].dep = cur_depth;
                if (USEDB) {
                    tr[n + cnt].bv = tr[n - 1].bv;
                    tr[n + cnt].prex = tr[n - 1].prex;
                }
                tr[n - 1].fa = n + cnt;
            } else {
                tr[n + cnt].rson = n - 1;
                tr[n + cnt].dep = cur_depth;
                tr[n - 1].fa = n + cnt;
                if (USEDB) {
                    if (fun) {
                        if (tr[n - 1].bv.size() >= tr[tr[n + cnt].lson].bv.size()) {
                            tr[n + cnt].bv = tr[tr[n + cnt].lson].bv;
                            tr[n + cnt].bv.resize(tr[n - 1].bv.size());
                            tr[n + cnt].bv |= tr[n - 1].bv;
                        }
                        else {
                            tr[n + cnt].bv = tr[n - 1].bv;
                            tr[n + cnt].bv.resize(tr[tr[n + cnt].lson].bv.size());
                            tr[n + cnt].bv |= tr[tr[n + cnt].lson].bv;
                        }
                    }
                    else {
                        if (fun) {
                            tr[n + cnt].bv = tr[n - 1].bv;
                            tr[n + cnt].bv.resize(tr[tr[n + cnt].lson].bv.size());
                            tr[n + cnt].bv &= tr[tr[n + cnt].lson].bv;
                        }
                        else {
                            tr[n + cnt].bv = tr[tr[n + cnt].lson].bv;
                            tr[n + cnt].bv.resize(tr[n - 1].bv.size());
                            tr[n + cnt].bv &= tr[n - 1].bv;
                        }
                    }
                    /*if (fun) {
                        if (tr[n - 1].prex * 32 + tr[n - 1].bv.size() >= tr[tr[n + cnt].lson].prex * 32 + tr[tr[n + cnt].lson].bv.size()) {
                            tr[n + cnt].bv = tr[n - 1].bv;
                            tr[n + cnt].prex = min(tr[n - 1].prex, tr[tr[n + cnt].lson].prex);
                            int commonSize = tr[n - 1].prex - tr[tr[n + cnt].lson].prex;
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() + 32 * commonSize);
                            for (int l = tr[tr[n + cnt].lson].bv.num_blocks() - 1, j = tr[n + cnt].bv.num_blocks() - 1; ~l && ~j; l--, j--)
                                tr[n + cnt].bv.m_bits[j] |= tr[tr[n + cnt].lson].bv.m_bits[l];
                            assert(tr[n + cnt].bv.m_bits.back() != 0);
                        }
                        else if (tr[n - 1].prex * 32 + tr[n - 1].bv.size() < tr[tr[n + cnt].lson].prex * 32 + tr[tr[n + cnt].lson].bv.size()) {
                            tr[n + cnt].bv = tr[tr[n + cnt].lson].bv;
                            tr[n + cnt].prex = min(tr[n - 1].prex, tr[tr[n + cnt].lson].prex);
                            int commonSize = tr[tr[n + cnt].lson].prex - tr[n - 1].prex;
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() + 32 * commonSize);
                            for (int l = tr[n - 1].bv.num_blocks() - 1, j = tr[n + cnt].bv.num_blocks() - 1; ~l && ~j; l--, j--)
                                tr[n + cnt].bv.m_bits[j] |= tr[n - 1].bv.m_bits[l];
                            assert(tr[n + cnt].bv.m_bits.back() != 0);
                        }
                    }
                    else {
                        if (tr[n - 1].prex * 32 + tr[n - 1].bv.size() <= tr[tr[n + cnt].lson].prex * 32 + tr[tr[n + cnt].lson].bv.size()) {
                            tr[n + cnt].bv = tr[n - 1].bv;
                            tr[n + cnt].prex = max(tr[n - 1].prex, tr[tr[n + cnt].lson].prex);
                            int commonSize = tr[n - 1].prex - tr[tr[n + cnt].lson].prex;

                            int l = tr[n - 1].bv.num_blocks() - 1 + max(0, commonSize), j = tr[n + cnt].bv.num_blocks() - 1;
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() - min(0, commonSize) * 32);
                            int zeroCnt = 0;
                            commonSize = max(0, commonSize);
                            for (; commonSize > 0; j--, commonSize--) {
                                tr[n + cnt].bv.m_bits[j] = 0;
                                zeroCnt++;
                            }

                            for (; ~l && ~j; l--, j--)
                                tr[n + cnt].bv.m_bits[j] &= tr[tr[n + cnt].lson].bv.m_bits[l];
                            for (j = tr[n + cnt].bv.num_blocks() - 1; ~j; j--) {
                                if (tr[n + cnt].bv.m_bits[j] != 0) {
                                    tr[n + cnt].prex += tr[n + cnt].bv.num_blocks() - j - 1;
                                    tr[n + cnt].bv.resize((j + 1) * 32);
                                    tr[n + cnt].bv.shrink_to_fit();
                                    break;
                                }
                            }
                            // assert(tr[n + cnt].bv.m_bits.back() != 0);
                        }
                        else {
                            tr[n + cnt].bv = tr[tr[n + cnt].lson].bv;
                            tr[n + cnt].prex = max(tr[n - 1].prex, tr[tr[n + cnt].lson].prex);
                            int commonSize = tr[tr[n + cnt].lson].prex - tr[n - 1].prex;

                            int l = tr[n - 1].bv.num_blocks() - 1 + min(0, commonSize), j = tr[n + cnt].bv.num_blocks() - 1;
                            tr[n + cnt].bv.resize(tr[n + cnt].bv.size() - min(0, commonSize) * 32);
                            int zeroCnt = 0;
                            commonSize = max(0, commonSize);
                            for (; commonSize > 0; j--, commonSize--) {
                                tr[n + cnt].bv.m_bits[j] = 0;
                                zeroCnt++;
                            }
                            for (; ~l && ~j; l--, j--)
                                tr[n + cnt].bv.m_bits[j] &= tr[n - 1].bv.m_bits[l];
                            for (j = tr[n + cnt].bv.num_blocks() - 1; ~j; j--) {
                                if (tr[n + cnt].bv.m_bits[j] != 0) {
                                    tr[n + cnt].prex += tr[n + cnt].bv.num_blocks() - j - 1;
                                    tr[n + cnt].bv.resize((j + 1) * 32);
                                    tr[n + cnt].bv.shrink_to_fit();
                                    break;
                                }
                            }
                            // assert(tr[n + cnt].bv.m_bits.back() != 0);
                        }
                    }#1#
                }
            }
            // cout << n + cnt << " " << tr[n + cnt].bv.count() << endl;
            cnt++;
        }*/
        n += cnt;
        cur_depth++;
    }
}

void build(Tree &tre, bool f, int k) {
    auto &tr = tre.tr;
    tre.leaf = tr.size();
    int depth = 0;

    size_t nums = 0;
    size_t nn = tr.size();
    while (nn > k) {
        depth++;
        nums += nn;

        if (nn % k != 0)
            nn = nn / k + 1;
        else
            nn /= k;
    }
    depth += 2;
    nums += nn + 1;
    // cout << "depth:" << depth << "nums:" << nums << endl;
    int cur_depth = 1;
    int n = tr.size();
    int half = depth - 2;
    if (half == -1)
        half = 0;
    // leaf_size.emplace_back(1 << half);
    sum += n;
    tr.resize(nums);
    for (size_t i = 0; i < nums; ++i) tr[i].id = i, tr[i].dep = 0;
    tr[nums - 1].fa = nums - 1;
    int cnt = n;
    int sz = 0;
    while (cur_depth < depth) {
        //// std::cout << "cur_depth_node_num:" << cnt << std::endl;
        int i = n - cnt + k;
        cnt = 0;
        for (; i <= n; i += k) {
            tr[n + cnt].sons.resize(k);
            // cout << "F: " << n + cnt << "- ";
            for (int j = 0; j < k; ++j) {
                tr[n + cnt].sons[j] = i - k + j;
                // cout << i - k + j << " ";
                tr[i - k + j].fa = n + cnt;
            }
            // tr[n + cnt].lson = i - 2;
            // tr[n + cnt].rson = i - 1;
            tr[n + cnt].dep = cur_depth;

            if (!f) {
                tr[n + cnt].bv = tr[i - 1].bv;
                for (int j = i - 2; j >= i - k; j--) {
                    // if (tr[n + cnt].bv.size() < tr[j].bv.size()) {
                        tr[n + cnt].bv.resize(tr[j].bv.size());
                        tr[n + cnt].bv &= tr[j].bv;
                    // }
                    // else {
                    //     sz = tr[j].bv.size();
                    //     tr[j].bv.resize(tr[n + cnt].bv.size());
                    //     tr[n + cnt].bv &= tr[j].bv;
                    //     tr[j].bv.resize(sz);
                    //     tr[j].bv.shrink_to_fit();
                    // }
                }
            }
            else {
                tr[n + cnt].bv = tr[i - k].bv;
                for (int j = i - k + 1; j <= i - 1; j++) {
                    // if (tr[n + cnt].bv.size() < tr[j].bv.size()) {
                        tr[n + cnt].bv.resize(tr[j].bv.size());
                        tr[n + cnt].bv |= tr[j].bv;
                    // }
                    // else {
                    //     sz = tr[j].bv.size();
                    //     tr[j].bv.resize(tr[n + cnt].bv.size());
                    //     tr[n + cnt].bv |= tr[j].bv;
                    //     tr[j].bv.resize(sz);
                    //     tr[j].bv.shrink_to_fit();
                    // }
                }
                // tr[n + cnt].bv = tr[i - 2].bv | tr[i - 1].bv;
            }
            cnt++;
        }
        if (i - k != n) {
            tr[n + cnt].sons.reserve(n - i + k);
            tr[n + cnt].dep = cur_depth;
            // cout << "F: " << n + cnt << "- ";
            for (int j = i - k; j < n; j++) {
                tr[n + cnt].sons.emplace_back(j);
                // cout << j << " ";
                tr[j].fa = n + cnt;
            }

            if (!f) {
                tr[n + cnt].bv = tr[tr[n + cnt].sons.back()].bv;
                for (int j = n - 2; j >= i - k; j--) {
                    // if (tr[n + cnt].bv.size() < tr[j].bv.size()) {
                        tr[n + cnt].bv.resize(tr[j].bv.size());
                        tr[n + cnt].bv &= tr[j].bv;
                    // }
                    // else {
                    //     sz = tr[j].bv.size();
                    //     tr[j].bv.resize(tr[n + cnt].bv.size());
                    //     tr[n + cnt].bv &= tr[j].bv;
                    //     tr[j].bv.resize(sz);
                    //     tr[j].bv.shrink_to_fit();
                    // }
                }
            } else {
                tr[n + cnt].bv = tr[tr[n + cnt].sons[0]].bv;
                for (int j = i - k + 1; j < n; j++) {
                    // if (tr[n + cnt].bv.size() < tr[j].bv.size()) {
                        tr[n + cnt].bv.resize(tr[j].bv.size());
                        tr[n + cnt].bv |= tr[j].bv;
                    // }
                    // else {
                    //     sz = tr[j].bv.size();
                    //     tr[j].bv.resize(tr[n + cnt].bv.size());
                    //     tr[n + cnt].bv |= tr[j].bv;
                    //     tr[j].bv.resize(sz);
                    //     tr[j].bv.shrink_to_fit();
                    // }
                }
                // tr[n + cnt].bv = tr[i - 2].bv | tr[i - 1].bv;
            }
            cnt++;
        }
        n += cnt;
        cur_depth++;
    }
}

void DG_method(bool fun) {
    auto start = Get_Time();
    // Graph query_graph{};
    // query_graph.LoadFromFile(query_path);
    // query_graph.PrintMetaData();
    int k = 4;
    // K--;
    DeltaGraph DG(k, K);
    if (fun) DG.fdif = "union";
    else DG.fdif = "intersection";
    std::cout << "----------- Building DeltaGraph ------------" << std::endl;
    DG.build(query_path, K);

    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    NodeDG res;

    level = 2;
    if (DG.depth <= 2) level = 1;
    if (level != 0)
        DG.graphpool(level);

    get_index_mem();
    cout << input_file << endl;
    fstream input(input_file, ios::in);
    fstream CSV("../result.csv", ios::app);
    uint L, R;
    int CNT = 0;

    auto getdist_time = Get_Time();
    std::cout << "----------- GetDist ------------" << std::endl;
    DG.GetDist();
    Print_Time("getDistTime: ", getdist_time);

    CSV << "DG" << ',' << fun << ',' << dataset << ',' << K << ',' << Duration(start) / 1000 << ',' << mem::getValue() / 1024 << ',';
    std::cout << "Build Over: " << mem::getValue() / 1024 << "Mb" << std::endl;

    size_t intSz = 0;
    for (const auto& i : DG.skeleton) {
        intSz += i.son.capacity();
        intSz += 4;
        intSz += 17 * 1.0 / 4 * i.Neighbors.size();
        intSz += 17 * 1.0 / 4 * i.Neighbor_delta.size();
        for (const auto& j : i.Neighbors) {
            intSz += j.second.capacity();
            for (auto k : j.second) {
                intSz += k.csr.capacity();
            }
        }
    }

    auto query_start_time = Get_Time();
    size_t asd = 0;
    int super_root = DG.skeleton.size() - 1;
    while (input >> L >> R) {
        auto begin_q = Get_Time();
        cout << "Query: " << L << " " << R << endl;
        CNT += 1;

        auto prim_time = Get_Time();

        std::cout << "----------- Prim ------------" << std::endl;
        auto p = DG.Prim(L, R);
        // for (const auto &i: p) {
        //     for (auto j: i) {
        //         cout << j << " ";
        //     }
        //     cout << endl;
        // }

        Print_Time("primTime: ", prim_time);

        double time1 = 0;
        double time2 = 0;
        NodeDG res;
        if (level == 0)
            res.Neighbors[-1] = DG.skeleton[DG.skeleton.back().son[0]].Neighbors[-1];
        map<int, bool> st;

        for (int i_ = 0; i_ < p.size(); i_++) {
            auto i = p[i_];
            if (i[0] == super_root) {
                if (i.size() > 2) {
                    NodeDG restmp = move(res);
                    // if (restmp.Neighbors.size() > 0) {
                    //     int ct = 0;
                    //     for (const auto& asd : restmp.Neighbors[-1]) {
                    //         ct += asd.csr.size();
                    //     }
                    //     cout << "776 " << ct << endl;
                    // }
                    if (res.Neighbors.empty() ) {
                        res.Neighbors[-1] = DG.skeleton[i[1]].Neighbors[-1];
                    }

                    // res.Neighbors[-1] = DG.skeleton[i[1]].Neighbors[-1];
                    // cout << " 479->" << i[1] << endl;
                    // int ct = 0;
                    // for (const auto& asd : res.Neighbors[-1])  ct += asd.csr.size();
                    //
                    // cout << "785 " << ct << endl;
                    bool quickL = false, quickR = false;
                    for (int j = 1; j + 1 < i.size(); j++) {
                        // if (fun) {
                            if (i[j] >= K) {
                                if (i_ == p.size() - 1 && DG.skeleton[i[j]].son[0] == i[j + 1] && R >= DG.skeleton[i[j]].son.back()) {
                                    quickL = true;
                                    // cout << "quickL" << endl;
                                    break;
                                }
                                if (i_ == 0 && DG.skeleton[i[j]].son.back() == i[j + 1] && L <= DG.skeleton[i[j]].son[0]) {
                                    quickR = true;
                                    // cout << "quickR" << endl;
                                    break;
                                }

                                if (fun) DG.cha(DG.skeleton[i[j]], res, i[j + 1]);
                                else DG.cup(res, DG.skeleton[i[j]], i[j + 1]);

                                // cout << i[j] << " 484->" << i[j + 1] << endl;

                            }
                            else {
                                if (fun) DG.cup(res, DG.skeleton[i[j]], i[j + 1]);
                                else DG.cha(DG.skeleton[i[j]], res, i[j + 1]);

                                if (DG.skeleton[i[j + 1]].fa != DG.skeleton[i[j]].fa) {
                                    if (fun) DG.cha(DG.skeleton[i[j + 1]], res, i[j]);
                                    else DG.cup(res, DG.skeleton[i[j + 1]], i[j]);
                                }

                                // cout << i[j] << " 487->" << i[j + 1] << endl;
                            }
                        // int ct = 0;
                        // for (auto asd : res.Neighbors[-1]) {
                        //     ct += asd.csr.size();
                        // }
                        // cout << ct << endl;
                        // }
                    }
                    if (i.back() == R && !quickR)   {
                        auto& ve =  DG.skeleton[DG.skeleton[i.back()].fa].son;
                        for (auto j = ve.size() - 1; ~j; j--) {
                            if (ve[j] < R) {
                                if (fun) DG.cup(res, DG.skeleton[ve[j + 1]], ve[j]);
                                else DG.cha(DG.skeleton[ve[j + 1]], res, ve[j]);
                                // cout << ve[j + 1] << " 496->" << ve[j] << endl;
                                // int ct = 0;
                                // for (auto asd : res.Neighbors[-1]) {
                                //     ct += asd.csr.size();
                                // }
                                // cout << ct << endl;
                            }
                        }
                    }
                    if (i.back() == L && !quickL) {
                        auto& ve =  DG.skeleton[DG.skeleton[i.back()].fa].son;
                        for (auto j = 0; j + 1 < ve.size(); j++) {
                            if (L <= ve[j] && ve[j + 1] < R) {
                                if (fun) DG.cup(res, DG.skeleton[ve[j]], ve[j + 1]);
                                else DG.cha(DG.skeleton[ve[j]], res, ve[j + 1]);

                                // cout << ve[j] << " 505->" << ve[j + 1] << endl;
                                // int ct = 0;
                                // for (auto asd : res.Neighbors[-1]) {
                                //     ct += asd.csr.size();
                                // }
                                // cout << ct << endl;
                            }
                        }
                    }
                    if (fun) DG.cup(res, restmp, -1);
                    else {
                        if (!restmp.Neighbors.empty())
                            DG.cap(res, restmp);
                    }
                    // int ct = 0;
                    // for (auto asd : res.Neighbors[-1]) {
                    //     ct += asd.csr.size();
                    // }
                    // cout << ct << endl;
                }
                if (i.size() == 2) {
                    // assert(res.Neighbors.empty());
                    // st[i[1]] = true;
                    if (fun) DG.cup(res, DG.skeleton[i[1]], -1);
                    else DG.cap(res, DG.skeleton[i[1]]);
                    // cout << " 514->" << i[1] << endl;
                    // int ct = 0;
                    // for (auto asd : res.Neighbors[-1]) {
                    //     ct += asd.csr.size();
                    // }
                    // cout << ct << endl;
                }
            }
            else {
                for (int j = 0; j + 1 < i.size(); ++j) {
                    if (DG.skeleton[j].id > K) continue;
                    // st[i[j + 1]] = true;
                    if (fun)
                        DG.cup(res, DG.skeleton[i[j]], i[j + 1]);
                    else {
                        DG.cha(DG.skeleton[i[j]], res, i[j + 1]);
                    }
                    // int ct = 0;
                    // for (auto asd : res.Neighbors[-1]) {
                    //     ct += asd.csr.size();
                    // }
                    // cout << ct << endl;
                    // cout << i[j] << " 523->" << i[j + 1] << endl;
                }
            }
        }

        for (const auto &i: res.Neighbors[-1]) all_edges += i.csr.size();

        Print_Time("queryTime: ", begin_q);
        std::cout << "----------------------------------------------Res: " << all_edges << std::endl;
        cout << "peek memoty: " << mem::getValue() / 1024 << "Mb" << std::endl;

        all_edges = 0;
        if (CNT == 15) {
            CSV << Duration(query_start_time) / 1000 << '\t';
        }
    }
    input.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_start_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << std::endl;
    CSV << Duration(query_start_time) / 1000 << ',' << mem::getValue() / 1024;
    // if (level != 0)
    //     CSV << ',' << level;
    CSV << "," << intSz << std::endl;
    CSV.close();
}

void Base_method(bool fun) {
    auto start = Get_Time();
    std::vector<std::vector<std::pair<int, int>>> snaps(K);
    size_t csr_edge_count = 0;
    //    Graph* snap = new Graph[K];
    for (int snap_i = 0; snap_i < K; ++snap_i) {
        // snaps[snap_i].LoadFromFile(query_path + std::to_string(snap_i), true);
        string path = query_path + std::to_string(snap_i);;
        if (!io::file_exists(path.c_str()))
        {
            std::cout << "Failed to open: " << path << std::endl;
            exit(-1);
        }
        std::ifstream ifs(path);

        std::string type;
        uint from_id, to_id;
        while (ifs >> from_id >> to_id)
        {
            snaps[snap_i].resize(snaps[snap_i].size() + 1);
            if (from_id > to_id)
                snaps[snap_i].back() = {to_id, from_id};
            else
                snaps[snap_i].back() = {from_id, to_id};
        }
        ifs.close();

        sort(snaps[snap_i].begin(), snaps[snap_i].end());

    }
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    get_index_mem();

    uint L, R;
    cout << input_file << endl;
    fstream input(input_file, ios::in);
    fstream CSV("/home/qsl/exp/tree/SM.csv", ios::app);

    int CNT = 0;

    size_t intSz = 0;
    for (int snap_i = 0; snap_i < K; ++snap_i) {
        intSz += snaps[snap_i].size() * 2;
    }

    auto query_time = Get_Time();
    auto query_start_time = Get_Time();
    while (input >> L >> R) {
        cout << "Query: " << L << " " << R;
        CNT++;
        std::vector<std::pair<int, int>> res;
        for (int snap_i = L; snap_i <= R; ++snap_i) {
            if (fun) {
                // for (int i = 0; i < snaps[snap_i].neighbors_.size(); ++i) {
                    std::vector<pair<int, int>> tmp;
                    std::set_union(res.begin(), res.end(),
                                   snaps[snap_i].begin(), snaps[snap_i].end(),
                                   std::back_inserter(tmp));
                    res = move(tmp);
                // }
            }
            else {
                // for (int i = 0; i < snaps[snap_i].neighbors_.size(); ++i) {
                    if (snap_i == L) {
                        res = snaps[snap_i];
                    } else {
                        vector<pair<int, int>> tmp;
                        std::set_intersection(res.begin(), res.end(),
                                              snaps[snap_i].begin(), snaps[snap_i].end(),
                                              std::back_inserter(tmp));
                        res = move(tmp);
                    }
                // }
            }
        }

        for (const auto &i: res) {
            all_edges++;
        }
        if (res.size() >= 3) {
            // string s = "CECI," + dataset + "," + to_string(L) + "," + to_string(R) + ",";
            // auto cecit = Get_Time();
            auto ceciRes = CECI(res);
            auto Vcount = ceciRes[0];
            auto Ecount = ceciRes[1];
            auto embedding_count = ceciRes[2];
            assert(embedding_count % 6 == 0);
            embedding_count /= 6;
            // CECI_information[L - R + 1].emplace_back();
            // std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
            //        CSV << "TR" << ',' << L << ',' << R << ',' << all_edges << endl;
            // cout << "peek memoty: " << mem::getValue() / 1024 << "Mb" << std::endl;
            CSV << L << "," << R << "," << R - L + 1 << "," << Vcount << "," << Ecount << "," << embedding_count << "," << Duration(query_time) << endl;
            // cout << L << "," << R << "," << R - L + 1 << "," << Vcount << "," << Ecount << "," << embedding_count << "," << Duration(query_time) << endl;
            // cout << embedding_count << endl;
        }
        std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
        //        CSV << "TR" << ',' << L << ',' << R << ',' << all_edges << endl;
        all_edges = 0;
    }
    input.close();

    CSV.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb";
}

void ViLa_method(bool fun) {
    auto start = Get_Time();
    ViLa_Graph ViLa;
    ViLa.build(K, query_path);
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    get_index_mem();

    uint L, R;
    fstream input(input_file, ios::in);
    fstream CSV("/home/qsl/exp/tree/SM.csv", ios::app);

    // CSV << "ViLa" << ',' << fun << ',' << dataset << ',' << K << ',' << Duration(start) / 1000 << ',' << mem::getValue() / 1024 << ',';

    int CNT = 0;
    size_t intSz = 0;
    for (auto& i : ViLa.ViLaEdge) {
        for (auto& j: i) {
            intSz ++;
            intSz += j.lifespan.num_blocks();
        }
    }

    auto query_time = Get_Time();
    auto query_start_time = Get_Time();
    cout << "begin Query" << endl;
    vector<pair<int, int> > res;
    res.reserve(10000);
    map<int, vector<int>> vmp;
    map<int, map<int, vector<int>>> clique_ifm;
    for (int L = 0; L < K; L++)
        // while (input >> L >> R)
            for (int R = L + 1; R < K; R++)
                if (R - L + 1 < 7)
        {
        // cout << "Query: " << L << " " << R;
        CNT++;
        res.clear();
        res.reserve(100000);
        bool pas = true;
        dynamic_bitset<uint32_t> bit(R - L + 1);
        for (int mp = 0; mp < ViLa.ViLaEdge.size(); mp++) {
            for (const auto &ed: ViLa.ViLaEdge[mp]) {
                string s;
                to_string(ed.lifespan, s);
                s = s.substr(s.size() - 1 - R, R - L + 1);
                // s = s.substr(L, R - L + 1);
                if (fun && s.find('1') != string::npos) {
                    res.emplace_back(mp, ed.nei_id);
                }
                if (!fun && s.find('0') == string::npos) {
                    res.emplace_back(mp, ed.nei_id);
                    pas = false;
                }
            }
        }
        // sort(res.begin(), res.end());
        // res.erase(unique(res.begin(), res.end()), res.end());
        int all_edges = 0;
        for (auto& i : res )
            all_edges++;
        std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;

        if (res.size() >= 3) {
            // string s = "CECI," + dataset + "," + to_string(L) + "," + to_string(R) + ",";
            // auto cecit = Get_Time();
            for (int cliquek = 3; cliquek <= 7; cliquek++) {
                auto ceciRes = clique(res, cliquek);
                // auto embedding_count = ceciRes[2];
                // assert(embedding_count % 6 == 0);
                // embedding_count /= 6;
                clique_ifm[R - L + 1][cliquek].emplace_back(ceciRes);
                // CECI_information[L - R + 1].emplace_back();
                // std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
                //        CSV << "TR" << ',' << L << ',' << R << ',' << all_edges << endl;
                // cout << "peek memoty: " << mem::getValue() / 1024 << "Mb" << std::endl;
                CSV << L << "," << R << "," << R - L + 1 << "," << cliquek << "," << ceciRes << endl;
                // cout << L << "," << R << "," << R - L + 1 << "," << cliquek << "," << ceciRes << endl;
                // if (L == 40 && R == 44 && cliquek == 7) {
                //     cout << endl;
                // }
            }
            // cout << L << "," << R << "," << R - L + 1 << "," << Vcount << "," << Ecount << "," << embedding_count << "," << Duration(query_time) << endl;
            // cout << embedding_count << endl;
        }

        // if (CNT == 15) {
        //     CSV << Duration(query_start_time) / 1000 << ',';
        //     // cout << "ViLa" << ',' << mem::getValue() / 1024 << "mb" << ',' << Duration(query_start_time) / 1000 << "s"
        //     //         << std::endl;
        // }
    }
    input.close();
    for (auto i : clique_ifm) {
        for (auto j : i.second) {
            int sum = 0;
            for (auto l : j.second) {
                sum += l;
            }
            CSV << i.first << "," << j.first << "," << sum << "," << 1.0*sum / j.second.size() << endl;
        }
    }
    // std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000
    //         << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << std::endl;
    // CSV << Duration(query_start_time) / 1000 << ',' << mem::getValue() / 1024 << "," << intSz * 4 << std::endl;

    CSV.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb";
}

mt19937 rd(time(nullptr));

void getK() {
    if (dataset == "mo1H") K = 34920;
        if (dataset == "mo12H") K = 4693;
        if (dataset == "mo1D") K = 2350;
        if (dataset == "mo3D") K = 784;
        if (dataset == "mo7D") K = 336;
        if (dataset == "mo15D") K = 157;
        if (dataset == "mo1M") K = 79;
        if (dataset == "mo2M") K = 40;
        if (dataset == "mo4M") K = 20;
        if (dataset == "au1H") K = 41804;
        if (dataset == "au12H") K = 4080;
        if (dataset == "au1D") K = 2046;
        if (dataset == "au3D") K = 684;
        if (dataset == "au7D") K = 294;
        if (dataset == "au15D") K = 138;
        if (dataset == "au1M") K = 70;
        if (dataset == "au2M") K = 36;
        if (dataset == "au4M") K = 19;
        if (dataset == "su1H") K = 56232;
        if (dataset == "su12H") K = 4892;
        if (dataset == "su1D") K = 2493;
        if (dataset == "su3D") K = 866;
        if (dataset == "su7D") K = 387;
        if (dataset == "su15D") K = 185;
        if (dataset == "su1M") K = 93;
        if (dataset == "su2M") K = 47;
        if (dataset == "su4M") K = 24;
        if (dataset == "so1H") K = 66445;
        if (dataset == "so12H") K = 5548;
        if (dataset == "so1D") K = 2775;
        if (dataset == "so3D") K = 925;
        if (dataset == "so7D") K = 397;
        if (dataset == "so15D") K = 185;
        if (dataset == "so1M") K = 93;
        if (dataset == "so2M") K = 47;
        if (dataset == "so4M") K = 24;
}

int GET_CNT_Read = 0;
double GET_CNT_TIME_Read = 0;
int GET_CNT_Write = 0;
double GET_CNT_TIME_Write = 0;

// inline void GetValFromDB(int DBidx, int key, TNode& tr) {
//     Status s;
//     Options op;
//     // cout << "../../DB/" + dataset + "/bv" + to_string(key) + "_db" << " " << key << endl;
//     // TNode as;
//
//     s = DB::Open(op, "../../DB/" + dataset + "/bv" + to_string(DBidx) + "_db", &bv_db[DBidx]);
//     if (!s.ok()) { cout << "bv_db open failed" << endl; cout << "../../DB/" + dataset + "/bv" + to_string(DBidx) + "_db" << endl; exit(1); }
//     uint l = 0;
//     string value;
//     auto getTi = Get_Time();
//     bv_db[DBidx]->Get(ReadOptions(), to_string(key), &value);
//     GET_CNT_TIME_Read += Duration(getTi);
//     GET_CNT_Read++;
//
// //    cout << value << endl;
//     int idx = 0;
//     // tr.prex = -1;
//     for (int sidx = 0; sidx < value.size(); ++sidx) {
//         if (value[sidx] == ',') {
//             // if (tr.prex == -1) tr.prex = l;
//             // else {
//                 // cout << l << " " << tr.
//                 tr.bv.append(l);
//                 idx++;
//             // }
//             // cout << l << endl;
//             l = 0;
//             continue;
//         }
//         l = l * 10 + value[sidx] - '0';
//     }
//     tr.bv.append(l);
//     tr.bv.shrink_to_fit();
//     bv_db[DBidx]->Close();
// }

int main(int argc, char *argv[]) {
    auto memst = mem::getValue();
    cout << memst << endl;

    time_t timep;
    time(&timep);
    printf("%s", ctime(&timep));
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
    CLI::App app{"App description"};
    string method = "1";
    dataset = "so1M";
    level = 0;
    fun = 0;
    int k = 2;
    EDGE_MAX = 64;
    app.add_option("-d,--dataset", dataset, "query graph path")->required();
    app.add_option("-m,--method", method, "method")->required();
    app.add_option("-f,--fun", fun, "method")->required();
    app.add_option("-k,--K", k, "K");
    // app.add_option("-l,--mat", level, "level");
    app.add_option("-e,--edgemax", EDGE_MAX, "edgemax");
    // CLI11_PARSE(app, argc, argv);

    getK();

    USEDB = true;
    query_path = "/home/qsl/exp/tree/dataset/" + dataset + "/q";
    input_file = "/home/qsl/exp/tree/dataset/" + dataset + "/inputCECI.txt";

    // input_file = "../dataset/" + dataset + "/inputcore.txt";
    // input_file = "../dataset/" + dataset + "/inputreach.txt";

    std::chrono::high_resolution_clock::time_point start, lstart;
    cout << query_path << " " << input_file << " " << method << " K=" << K << " " << EDGE_MAX << endl;
    start = Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;
    if (method == "2") {
        std::cout << "delta graph" << std::endl;
        DG_method(fun);
        return 0;
    }
    if (method == "0") {
        std::cout << "base" << std::endl;
        Base_method(fun);
        return 0;
    }
    if (method == "3") {
        std::cout << "ViLa" << std::endl;
        ViLa_method(fun);
        return 0;
    }
    if (method == "4") {
        return 0;
    }
    if (method == "5") {

        return 0;
    }
    std::cout << "tree" << std::endl;
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;

    // Status s;
    // Options op;
    // op.create_if_missing = true;
    // op.IncreaseParallelism();
    // op.OptimizeLevelStyleCompaction();

    int cnt = 0;
    size_t cur_edges = 0;
    size_t csr_edge_count = 0;
    size_t snap_sum = 0;

    snap_idx.emplace_back(0);
    tree.resize((K + EDGE_MAX - 1) / EDGE_MAX);
    for (auto& i : tree) i.tr.resize((k * EDGE_MAX - 1) / (k - 1));
    int tree_idx = 0;

    auto tmp = EDGE_MAX;
    while(tmp) {
        tmp /= k;
        depth++;
    }

    for (int snap_i = 0; snap_i < K; ++snap_i) {
        string path = query_path + std::to_string(snap_i);
        if (!io::file_exists(path.c_str()))
        {
            std::cout << "Failed to open: " << path << std::endl;
            exit(-1);
        }

        auto &b = tree[tree_idx].tr[cnt].bv;
        tree.back().tr.back().lson = -1;
        tree.back().tr.back().rson = -1;
        tree.back().tr.back().fa = -1;
        size_t edge_count_ = 0;
        std::ifstream ifs(path);

        uint v1, v2;
        while (ifs >> v1 >> v2)
        {
            edge_count_++;
            all_edges++;
            if (v1 >= CSR.size()) CSR.resize(v1 + 1);

            auto lower = lb(CSR[v1], v2);
            if (lower != CSR[v1].end() && (*lower).first == v2) {
                if ((*lower).second > b.size()) b.resize((*lower).second);
                b.set((*lower).second - 1);
            } else {
                CSR[v1].insert(lower, {v2, ++csr_edge_count});
                if (csr_edge_count > b.size()) b.resize(csr_edge_count);
                b.set(csr_edge_count - 1);
                id2edge.resize(id2edge.size() + 1);
                id2edge.back() = {v1, v2};
            }

            // auto lower = lower_bound(CSR[v1].begin(), CSR[v1].end(), v2);
            // if (lower != CSR[v1].end() && (*lower) == v2) {
            //     if ( > b.size()) b.resize((*lower));
            //     b.set((*lower) - 1);
            // } else {
            //     csr_edge_count++;
            //     if (csr_edge_count > b.size()) b.resize(csr_edge_count);
            //     b.set(csr_edge_count - 1);
            //     id2edge.resize(id2edge.size() + 1);
            //     // id2edge.reserve(id2edge.size() + 1);
            //     // auto p = lower - CSR[v1].begin();
            //     // cout << lower - CSR[v1].begin() << endl;
            //     id2edge.back() = {v1, static_cast<unsigned>((size_t)(lower - CSR[v1].begin())), (int)csr_edge_count};
            //     CSR[v1].insert(lower, v2);
            // }
        }
        ifs.close();
        cur_edges++;
        cnt++;
        // b.resize((b.size() + 31) / 32 * 32);
        // cout << edge_count_ << " " << b.count() << " " << b.num_blocks() << "->";

        // reverse(b.m_bits.begin(), b.m_bits.end());
        // for (int j = b.num_blocks() - 1; ~j; j--) {
        //     if (b.m_bits[j] != 0) {
        //         tree.back().tr.back().prex = b.num_blocks() - j - 1;
        //         b.resize((j + 1) * sizeof(uint32_t) * 8);
        //         b.shrink_to_fit();
        //         break;
        //     }
        // }
        // cout << b.num_blocks() << " " << b.count() << " " << tree.back().tr.back().prex << endl;

        /*for (int i = 0; i < snap.neighbors_.size(); ++i) {
            uint v1 = snap.GetVertexLabel(i);
            for (auto j: snap.neighbors_[i]) {
                uint v2 = snap.GetVertexLabel(j);
                if (v1 >= CSR.size()) CSR.resize(v1 + 1);

                auto lower = lb(CSR[v1], v2);
                if (lower != CSR[v1].end() && (*lower).first == v2) {
                    if ((*lower).second > b.size()) b.resize((*lower).second);
                    b.set((*lower).second - 1);
                } else {
                    CSR[v1].insert(lower, {v2, ++csr_edge_count});
                    if (csr_edge_count > b.size()) b.resize(csr_edge_count);
                    b.set(csr_edge_count - 1);
                    id2edge.resize(id2edge.size() + 1);
                    id2edge.back() = {v2, v1};
                }
            }
        }*/

        if (cur_edges >= EDGE_MAX) {
            // std::cout << "----------- Building tree ------------" << std::endl;
            // cout << "tree size" << " " << tree.back().tr.size() << endl;
            snap_sum += cnt;
            snap_idx.emplace_back(snap_sum);
            // build(tree[tree_idx], fun, k);
            build(tree[tree_idx], fun);
            cnt = 0;
            cur_edges = 0;
            tree_idx++;
            // std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
        }
    }

    if (cur_edges) {
        // std::cout << "----------- Building tree ------------" << std::endl;
        // cout << "tree size" << " " << tree.back().tr.size() << endl;
        snap_sum += cnt;
        snap_idx.emplace_back(snap_sum);

        // build(tree[tree_idx], fun, k);
        build(tree[tree_idx], fun);
        cur_edges = 0;
    }

    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;

    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;

    std::cout << "----------- Snap Query ------------" << std::endl;
    uint L, R;

    cout << input_file << endl;
    fstream input(input_file, ios::in);
    fstream CSV("/home/qsl/exp/tree/result.csv", ios::app);

    int CNT = 0;

    size_t one;
    // map<int, vector<pair<int, int>>> CECI_information;
    // map<int, map<int, vector<int>>> clique_ifm;
    //
    // for (int L = 0; L < K; L++)
    while (input >> L >> R)
        // for (int R = L + 1; R < K; R++)
            // if (R - L + 1 <= 7)
    // L = 40, R = 45;
    {
        auto query_time = Get_Time();
        // cout << "Query: " << L << " " << R << " ";
        assert(L < K);
        assert(R < K);
        CNT++;
        uint Ltree, Rtree; // L,R: snap idx Ltree,Rtree: tree idx about snap
        auto lower1 = std::lower_bound(snap_idx.begin(), snap_idx.end(), L);
        if (*lower1 != L) Ltree = lower1 - snap_idx.begin() - 1;
        else Ltree = lower1 - snap_idx.begin();
        auto lower2 = std::lower_bound(snap_idx.begin(), snap_idx.end(), R);
        if (*lower2 != R) Rtree = lower2 - snap_idx.begin() - 1;
        else Rtree = lower2 - snap_idx.begin();
        // cout << L << " " << Ltree << " " << R << " " << Rtree << " ";
        dynamic_bitset<uint32_t> lres;

        if (Ltree == Rtree) { // 一起往上跳
            auto &tr = tree[Ltree].tr;
            int lidx = L - snap_idx[Ltree], ridx = R - snap_idx[Rtree];
            // cout << snap_idx[Ltree] << " " << snap_idx[Rtree] << endl;
            if (ridx - lidx + 1 == snap_idx[Rtree + 1] - snap_idx[Ltree]) {
                lres = tr.back().bv;
                // cout << "quick " << tr.back().id << endl;
                goto RES;
            }

            while (true) {
                if (lidx + 1 == ridx) {
                    if (tr[lidx].fa != tr[ridx].fa) {
                        if (lres.empty()) {
                            lres = tr[lidx].bv;
                            int sz = -1;
                            if (fun) {
                                if (lres.size() > tr[ridx].bv.size()) sz = tr[ridx].bv.size(), tr[ridx].bv.resize(lres.size());
                                else if (lres.size() < tr[ridx].bv.size()) lres.resize(tr[ridx].bv.size());
                                lres |= tr[ridx].bv;
                            }
                            else {
                                lres.resize(tr[ridx].bv.size());
                                lres &= tr[ridx].bv;
                            }
                            if (sz != -1) tr[ridx].bv.resize(sz);
                        }
                        else {
                            int sz = -1;
                            if (fun) {
                                if (lres.size() > tr[lidx].bv.size()) sz = tr[lidx].bv.size(), tr[lidx].bv.resize(lres.size());
                                else if (lres.size() < tr[lidx].bv.size()) lres.resize(tr[lidx].bv.size());
                                lres |= tr[lidx].bv;
                            }
                            else {
                                lres.resize(tr[lidx].bv.size());
                                lres &= tr[lidx].bv;
                            }
                            if (sz != -1) tr[lidx].bv.resize(sz);

                            sz = -1;
                            if (fun) {
                                if (lres.size() > tr[ridx].bv.size()) sz = tr[ridx].bv.size(), tr[ridx].bv.resize(lres.size());
                                else if (lres.size() < tr[ridx].bv.size()) lres.resize(tr[ridx].bv.size());
                                lres |= tr[ridx].bv;
                            }
                            else {
                                lres.resize(tr[ridx].bv.size());
                                lres &= tr[ridx].bv;
                            }
                            if (sz != -1) tr[ridx].bv.resize(sz);
                        }
                        break;
                    }
                }
                if (lidx == ridx) {
                    if (lres.empty()) lres = tr[lidx].bv;
                    else {
                        if (fun) {
                            if (lres.size() < tr[lidx].bv.size()) lres.resize(tr[lidx].bv.size());
                            else if (lres.size() > tr[lidx].bv.size()) tr[lidx].bv.resize(lres.size());
                            lres |= tr[lidx].bv;
                        }
                        else {
                            lres.resize(tr[lidx].bv.size());
                            lres &= tr[lidx].bv;
                        }
                    }
                    break;
                }
                if (tr[tr[lidx].fa].rson == lidx) {
                    if (lres.empty()) lres = tr[lidx].bv;
                    else {
                        if (fun) {
                            if (lres.size() < tr[lidx].bv.size()) lres.resize(tr[lidx].bv.size());
                            else if (lres.size() > tr[lidx].bv.size()) tr[lidx].bv.resize(lres.size());
                            lres |= tr[lidx].bv;
                        }
                        else {
                            lres.resize(tr[lidx].bv.size());
                            lres &= tr[lidx].bv;
                        }
                    }
                    lidx = tr[lidx].fa + 1;
                }
                else lidx = tr[lidx].fa;
                if (tr[tr[ridx].fa].lson == ridx) {
                    if (lres.empty()) lres = tr[ridx].bv;
                    else {
                        if (fun) {
                            if (lres.size() < tr[ridx].bv.size()) lres.resize(tr[ridx].bv.size());
                            else if (lres.size() > tr[ridx].bv.size()) tr[ridx].bv.resize(lres.size());
                            lres |= tr[ridx].bv;
                        }
                        else {
                            lres.resize(tr[ridx].bv.size());
                            lres &= tr[ridx].bv;
                        }
                    }
                    ridx = tr[ridx].fa - 1;
                }
                else {
                    ridx = tr[ridx].fa;
                }
            }

            goto RES;
        }

        // 左边
        // if (L != snap_idx[Ltree])
        {
            auto &tr = tree[Ltree].tr;
            uint idx = L - snap_idx[Ltree];
            uint bias = 0;
            uint root = tr.back().id;
            if (idx == 0) {
                lres = tr[root].bv;
                // cout << "Lcal: " << idx << endl;
            } else if (idx == snap_idx[Ltree + 1] - snap_idx[Ltree] - 1) {
                lres = tr[idx].bv;
                // cout << "Lcal: " << idx << endl;
            } else {
                while (1) {
                    bool brk = true;
                    bool brk2 = false;
                    if (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].rson == tr[idx].id && tr[idx].id != tr.back().id && tr[tr[idx].fa + 1].dep == tr[tr[idx].fa].dep) {
                        while (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].rson == tr[idx].id && tr[idx].id != tr.back().id) {
                            if (lres.empty()) {lres = tr[idx].bv;}
                            else {
                                if (fun) {
                                    if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                                    else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                                    lres |= tr[idx].bv;
                                }
                                else {
                                    lres.resize(tr[idx].bv.size());
                                    lres &= tr[idx].bv;
                                }
                            }
                            // all_edges = 0;
                            // one = lres.find_first();
                            // while (one != lres.npos) {
                            //     all_edges++;
                            //     one = lres.find_next(one);
                            // }
                            // cout << all_edges << endl;
                            // cout << "Lcal: " << idx << endl;
                            if (tr[tr[idx].fa + 1].dep != tr[tr[idx].fa].dep) {
                                brk2 = true;
                                break;
                            }
                            idx = tr[idx].fa + 1;
                        }
                        brk = false;
                    }
                    if (tr[idx].fa < tr.back().id && tr[tr[idx].fa].lson == tr[idx].id && tr[idx].id != tr.back().id) {
                        while (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].lson == tr[idx].id && tr[idx].id != tr.back().id) {
                            idx = tr[idx].fa;
                        }
                        if (lres.empty()) {lres = tr[idx].bv;}
                        else {
                            if (fun) {
                                if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                                else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                                lres |= tr[idx].bv;
                            }
                            else {
                                lres.resize(tr[idx].bv.size());
                                lres &= tr[idx].bv;
                            }
                        }
                        brk = false;
                        // cout << "Lcal: " << idx << endl;
                    }
                    if (brk) break;
                    if (brk2) break;
                }
            }
        }

        // 中间
        for (uint i = Ltree + 1; i < Rtree; ++i) {
            auto tmp = tree[i].tr.back().bv;

            size_t sz = -1;
            if (fun) {
                if (lres.size() > tmp.size()) {
                    sz = tmp.size();
                    tmp.resize(lres.size());
                }
                else if (lres.size() < tmp.size()) lres.resize(tmp.size());
                lres |= tmp;
            }
            else {
                lres.resize(tmp.size());
                lres &= tmp;
            }
            if (sz != -1) tmp.resize(sz);
            // cout << "cal mid" << endl;
        }
        // all_edges = 0;
        // one = lres.find_first();
        // while (one != lres.npos) {
        //     all_edges++;
        //     one = lres.find_next(one);
        // }
        // std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
        // 右边
        // if (R != snap_idx[Rtree + 1] - 1)
        {
            //            cout << "Right" << endl;
            auto &tr = tree[Rtree].tr;
            uint idx = R - snap_idx[Rtree];
            if (idx == (tr.size() + 1) / 2) {
                if (fun) {
                    if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                    else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                    lres |= tr[tr.back().id].bv;
                }
                else {
                    lres.resize(tr[idx].bv.size());
                    lres &= tr[tr.back().id].bv;
                }
                // cout << "Rcal: " << idx << endl;
            } else if (idx == 0) {
                if (fun) {
                    if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                    else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                    lres |= tr[idx].bv;
                }
                else {
                    lres.resize(tr[idx].bv.size());
                    lres &= tr[idx].bv;
                }
                // cout << "Rcal: " << idx << endl;
            } else {
                while (1) {
                    bool brk = true;
                    bool brk2 = true;
                    if (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].rson == tr[idx].id && tr[idx].id != tr.back().id) {
                        while (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].rson == tr[idx].id && tr[idx].id != tr.back().id) {
                            idx = tr[idx].fa;
                        }
                        brk = false;
                        if (fun) {
                            if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                            else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                            lres |= tr[idx].bv;
                        }
                        else {
                            lres.resize(tr[idx].bv.size());
                            lres &= tr[idx].bv;
                        }
                        // cout << "Rcal: " << idx << endl;
                    }
                    if (tr[idx].fa < tr.back().id && tr[tr[idx].fa].lson == tr[idx].id && tr[idx].id != tr.back().id &&
                        tr[tr[idx].fa - 1].dep == tr[tr[idx].fa].dep) {
                        while (tr[idx].fa <= tr.back().id && tr[tr[idx].fa].lson == tr[idx].id && tr[idx].id != tr.back().id) {
                            if (brk) {
                                if (fun) {
                                    if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                                    else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                                    lres |= tr[idx].bv;
                                }
                                else {
                                    lres.resize(tr[idx].bv.size());
                                    lres &= tr[idx].bv;
                                }
                                // cout << "Rcal: " << idx << endl;
                            } else
                                brk = true;
                            if (tr[tr[idx].fa - 1].dep != tr[tr[idx].fa].dep) {
                                brk2 = true;
                                break;
                            }
                            idx = tr[idx].fa - 1;
                        }
                        if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                        else if (lres.size() > tr[idx].bv.size()) tr[idx].bv.resize(lres.size());
                        brk = false;
                        }
                    if (brk) break;
                }
            }
        }
        RES:
            all_edges = 0;
        one = lres.find_first();
        // std::vector<std::pair<int, int>> E;
        while (one != lres.npos) {
            // E.emplace_back(id2edge[one]);
            all_edges++;
            one = lres.find_next(one);
        }
        // cout << " " << all_edges << endl;
        // if (E.size() >= 3) {
        //     // string s = "CECI," + dataset + "," + to_string(L) + "," + to_string(R) + ",";
        //     // auto cecit = Get_Time();
        //     for (int cliquek = 3; cliquek <= 7; cliquek++) {
        //         auto ceciRes = clique(E, cliquek);
        //         // auto embedding_count = ceciRes[2];
        //         // assert(embedding_count % 6 == 0);
        //         // embedding_count /= 6;
        //         clique_ifm[R - L + 1][cliquek].emplace_back(ceciRes);
        //         // CECI_information[L - R + 1].emplace_back();
        //         // std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
        //         //        CSV << "TR" << ',' << L << ',' << R << ',' << all_edges << endl;
        //         // cout << "peek memoty: " << mem::getValue() / 1024 << "Mb" << std::endl;
        //         CSV << L << "," << R << "," << R - L + 1 << "," << cliquek << "," << ceciRes << endl;
        //         cout << L << "," << R << "," << R - L + 1 << "," << cliquek << "," << ceciRes << endl;
        //     }
        //     // cout << L << "," << R << "," << R - L + 1 << "," << Vcount << "," << Ecount << "," << embedding_count << "," << Duration(query_time) << endl;
        //     // cout << embedding_count << endl;
        // }
    }
    // for (auto i : clique_ifm) {
    //     for (auto j : i.second) {
    //         int sum = 0;
    //         for (auto l : j.second) {
    //             sum += l;
    //         }
    //         CSV << i.first << "," << j.first << "," << sum << "," << j.second.size() << "," << 1.0*sum / j.second.size() << endl;
    //     }
    // }
    input.close();
    return 0;
}
