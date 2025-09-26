#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <bitset>
#include <cmath>
#include <algorithm>
#include "../utils/CLI11.hpp"
#include "../utils/globals.h"
#include "../utils/types.h"
#include <boost/dynamic_bitset.hpp>
#include "../graph/graph.h"
#include "matching.h"
#include "deltagraph.h"
#include "TiLa.h"
#include <fstream>
#include <rocksdb/db.h>
#include <rocksdb/options.h>
#include <rocksdb/slice.h>
//#include "Mem.h"
//#include "deltagraph.h"
using namespace std;
using namespace rocksdb;
int K = 18;
int level;
bool fun;
size_t EDGE_MAX = 3e4 * 256;
const int N = 100000000;
std::string query_path, dataset, input_file;
bool CHECK = true;
size_t all_edges;
using boost::dynamic_bitset;
struct TNode
{
    TNode() : lson(-1), rson(-1), fa(-1){
//        bv.resize(N, false);
    };
    dynamic_bitset<uint32_t> bv;
    int id;
    int glo_id;
    int lson;
    int rson;
    int fa;
    vector<int> ver_set;
    ~TNode() {
        bv.clear();
        bv.shrink_to_fit();
    }
};
size_t edge_idx;
std::vector<std::vector<std::pair<uint, uint>>> CSR;
std::vector<std::pair<uint, uint>> id2edge;
struct Tree{
    std::vector<TNode> tr;
};
std::vector<Tree> tree;
std::vector<int> one_index;
std::vector<uint> snap_idx;
std::vector<std::pair<uint, uint>> leaf_size;
size_t sum;
auto lb(std::vector<std::pair<uint, uint>>& V, uint val) {
    size_t L = 0, R = V.size();
    while (L < R) {
        size_t MID = (L + R) >> 1;
        if (V[MID].first >= val) R = MID;
        else L = MID + 1;
    }
    return V.begin() + L;
}

map<string, size_t> get_index_mem() {
    FILE* fp = fopen("/proc/self/status", "r");
    char line[128];
    map<string, size_t> res;
    while (fgets(line, 128, fp) != NULL)
    {
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

void build(Tree& tre, bool f) {
    auto &tr = tre.tr;
    int depth = log2(tr.size());
    if ((tr.size() & (tr.size()-1)) == 0)
        depth  ++;
    else
        depth += 2;
    size_t nums = 0;
    size_t nn = tr.size();
    while (nn) {
        nums += nn;
        if (nn == 1)
            break;
        nn = (nn + 1) / 2;
    }
    //cout << "depth:" << depth << "nums:" << nums << endl;
    int cur_depth = 1;
    int n = tr.size();
    int half = depth - 2;
    if (half == -1)
        half = 0;
    leaf_size.emplace_back(sum, 1 << half);
    sum += n;
    tr.resize(nums);
    for (size_t i = 0; i < nums; ++i) tr[i].id = i;
    tr[nums - 1].fa = nums - 1;
    int cnt = n;

    while (cur_depth < depth) {
        //// std::cout << "cur_depth_node_num:" << cnt << std::endl;
        int i = n - cnt + 2;
        cnt = 0;
        for (; i <= n; i += 2) {
            tr[n + cnt].lson = i - 2;
            tr[n + cnt].rson = i - 1;
            tr[i - 2].fa = n + cnt;
            tr[i - 1].fa = n + cnt;
            size_t sz = std::max(tr[i - 2].bv.size(), tr[i - 1].bv.size());
            //std::cout << "sz: " << sz << std::endl;
            int fl = 0;
            size_t sz_bk;
            if (tr[i - 2].bv.size() < sz) {
                sz_bk = tr[i - 2].bv.size();
                tr[i - 2].bv.resize(sz);
                fl = 1;
            }
            if (tr[i - 1].bv.size() < sz) {
                sz_bk = tr[i - 1].bv.size();
                tr[i - 1].bv.resize(sz);
                fl = 2;
            }
            if (tr[n + cnt].bv.size() < sz) {
                tr[n + cnt].bv.resize(sz);
            }
            if (!f) {
                tr[n + cnt].bv = tr[i - 2].bv & tr[i - 1].bv;
            } else {
                tr[n + cnt].bv = tr[i - 2].bv | tr[i - 1].bv;
            }
            if (fl == 1) {
                tr[i - 2].bv.resize(sz_bk);
            }
            if (fl == 2) {
                tr[i - 1].bv.resize(sz_bk);
            }
            int sum = 0;
////             std::cout << "fa:" << n + cnt << " " << tr[n+cnt].bv  << " lson:" << tr[n+cnt].lson << " " << tr[tr[n+cnt].lson].bv  << " rson:" << tr[n+cnt].rson << " " << tr[tr[n+cnt].rson].bv << std::endl;
            //std::cout << "fa:" << n + cnt << " lson:" << tr[n+cnt].lson << " rson:" << tr[n+cnt].rson << std::endl;
            //std::cout << "---" << std::endl;
            cnt++;
        }
        if (i - 2 != n) {
            if (tr[n+cnt].lson == -1) {
                tr[n+cnt].lson = n - 1;
                if (tr[n + cnt].bv.size() < tr[n - 1].bv.size()) {
                    tr[n + cnt].bv.resize(tr[n - 1].bv.size());
                }
                tr[n + cnt].bv = tr[n - 1].bv;
                tr[n - 1].fa = n + cnt;
            } else {
                tr[n+cnt].rson = n - 1;
                tr[n - 1].fa = n + cnt;
                size_t sz = std::max(tr[n - 1].bv.size(), tr[tr[n+cnt].lson].bv.size());
                int fl = 0;
                size_t sz_bk = 0;
                if (tr[n - 1].bv.size() < sz) tr[n - 1].bv.resize(sz), fl = 1;
                if (tr[tr[n+cnt].lson].bv.size() < sz) tr[tr[n+cnt].lson].bv.resize(sz), fl = 2;
                if (tr[n + cnt].bv.size() < sz) tr[n + cnt].bv.resize(sz);
                if (!f) {
                    tr[n + cnt].bv = tr[n - 1].bv & tr[tr[n+cnt].lson].bv;
                } else {
                    tr[n + cnt].bv = tr[n - 1].bv | tr[tr[n+cnt].lson].bv;
                }
                if (fl == 1) {
                    tr[n - 1].bv.resize(sz_bk);
                }
                if (fl == 2) {
                    tr[tr[n+cnt].lson].bv.resize(sz_bk);
                }
            }
////             std::cout << "fa:" << n + cnt << " " << tr[n+cnt].bv << " lson:" << tr[n+cnt].lson << " " << tr[tr[n+cnt].lson].bv  << " rson:" << tr[n+cnt].rson << " " << (tr[n+cnt].rson==-1 ? 0 : tr[tr[n+cnt].rson].bv) << std::endl;
            //std::cout << "fa:" << n + cnt << " lson:" << tr[n+cnt].lson << " rson:" << tr[n+cnt].rson << std::endl;
            //std::cout << "---" << std::endl;
            cnt ++;
        }
        n += cnt;
        cur_depth++;
    }
}

void DG_method(bool fun) {
    auto start = Get_Time();
    Graph query_graph {};
    // query_graph.LoadFromFile(query_path);
    // query_graph.PrintMetaData();
    int k = 2;
    DeltaGraph DG(k, K);
    if (fun) DG.fdif = "union";
    else DG.fdif = "intersection";
    std::cout << "----------- Building DeltaGraph ------------" << std::endl;
    DG.build(query_path, K);

    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;

    if (level != 0)
        DG.graphpool(level);

    get_index_mem();

    fstream input(input_file, ios::in);
    fstream CSV("../result.csv", ios::app);
    uint L, R;
    int CNT = 0;
    auto query_start_time = Get_Time();
    CSV << "DG" << '\t' << dataset << '\t' << K << '\t' << Duration(start) / 1000 << '\t' <<  mem::getValue() / 1024 << '\t';
    auto query_time = Get_Time();
    size_t asd = 0;

    std::cout << "Build Over: " << mem::getValue() / 1024  << "Mb" << std::endl;

    auto getdist_time = Get_Time();
    std::cout << "----------- GetDist ------------" << std::endl;
    DG.GetDist();
    Print_Time("getDistTime: ", getdist_time);

    while (input >> L >> R) {
        auto begin_q = Get_Time();
        cout << "Query: " << L << " " << R << endl;
        CNT += 1;

        auto prim_time = Get_Time();
        std::cout << "----------- Prim ------------" << std::endl;
        auto p = DG.Prim(L, R);
        for (const auto& i : p) {
            for (auto j : i) {
                cout << j << " ";
            }
            cout << endl;
        }
        Print_Time("primTime: ", prim_time);

        double time1 = 0;
        double time2 = 0;
        Node res;
        if (level == 0)
            res.Neighbors[-1] = DG.skeleton[DG.skeleton.back().son[0]].Neighbors[-1];
        bool has_mat = false;
        int idx_mat = 0;
        for (const auto& i: p) {
            for (auto j: i) {
                if (!has_mat && DG.skeleton[j].isMat) {
                    idx_mat = j;
                    has_mat = true;
                    break;
                }
            }
        }
        bool meet_mat = !has_mat;
//        cout << "has_mat " << has_mat << endl;
//        cout << "meet_mat " << meet_mat << endl;
//        cout << "idx_mat " << idx_mat << endl;
        for (size_t i = 1; i < p.size(); ++i) {
            if (p[i].back() == DG.skeleton.size() - 2) {
                auto t1 = Get_Time();
                for (size_t j = p[i].size() - 1; j; --j) {
                    if (level != 0 && DG.skeleton[p[i][j]].isMat == -1) continue;
                    if (level != 0 && DG.skeleton[p[i][j]].isMat == 1) res.Neighbors[-1] = DG.skeleton[p[i][j]].Neighbors[-1];;
                    if (fun) DG.cha(DG.skeleton[p[i][j]], res, p[i][j - 1]);
                    else DG.cup(res, DG.skeleton[p[i][j]], p[i][j - 1]);
                }
                time1 += Duration(t1);
            } else {
                auto t2 = Get_Time();
                if (p[i][0] == p[i - 1][0]) {}
                else reverse(p[i].begin(), p[i].end());
                for (size_t j = 0; j < p[i].size() - 1; ++j) {
                    if (level != 0 && DG.skeleton[p[i][j]].isMat == -1) continue;
                    if (level != 0 && DG.skeleton[p[i][j]].isMat == 1) res.Neighbors[-1] = DG.skeleton[p[i][j]].Neighbors[-1];;
                    if (fun) DG.cup(res, DG.skeleton[p[i][j]], p[i][j + 1]);
                    else DG.cha(DG.skeleton[p[i][j]], res, p[i][j + 1]);
                }
                time2 += Duration(t2);
            }
            // for (const auto& eac: res.Neighbors[-1]) all_edges += eac.csr.size();
            // cout << all_edges << endl;
            // all_edges = 0;
        }

        cout << "time1: " << time1 << "ms" << endl;
        cout << "time2: " << time2 << "ms" << endl;

//        fstream f("../res.txt", ios::app);
        for (const auto& i: res.Neighbors[-1]) all_edges += i.csr.size();

        Print_Time("queryTime: ", begin_q);
        std::cout << "----------------------------------------------Res: " << all_edges << std::endl;
        cout << "peek memoty: " << mem::getValue() / 1024  << "Mb" << std::endl;

        all_edges = 0;
        if (CNT == 5) {
            CSV << Duration(query_start_time) / 1000 << '\t';
            cout << "DG" << '\t' << mem::getValue() / 1024 << "mb" << '\t' << Duration(query_start_time) / 1000 << "s" << std::endl;
            query_start_time = Get_Time();
            CNT = 0;
        }
    }
    input.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000 << "s" << '\t'<< "ALL Time: " << Duration(start) / 1000 << std::endl;
    CSV << Duration(query_time) / 1000 << '\t' << Duration(start) / 1000 << '\t' << mem::getValue() / 1024;
    if (level != 0)
        CSV << '\t' << level;
    CSV << std::endl;
    CSV.close();
}

void Base_method(bool fun) {
    auto start = Get_Time();
    std::vector<Graph> snaps(K);
    size_t csr_edge_count = 0;
//    Graph* snap = new Graph[K];
    for (int snap_i = 0; snap_i < K; ++snap_i) {
        snaps[snap_i].LoadFromFile(query_path + std::to_string(snap_i));
    }
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    get_index_mem();

    uint L, R;
    fstream input(input_file, ios::in);
    fstream CSV("../result.csv", ios::app);

    int CNT = 0;
    CSV << "BS" << '\t' << dataset << '\t' << K << '\t' << Duration(start) / 1000 << '\t' <<  mem::getValue() / 1024 << '\t';
    auto query_time = Get_Time();
    auto query_start_time = Get_Time();
    while (input >> L >> R) {
        cout << "Query: " << L << " " << R << endl;
        CNT++;
        std::vector<Edge> res;
        map<int, int> st;
        for (int snap_i = L; snap_i <= R; ++snap_i) {
            if (fun) {
                for (int i = 0; i < snaps[snap_i].neighbors_.size(); ++i) {
                    auto v1 = snaps[snap_i].GetVertexLabel(i);
                    auto lower = std::lower_bound(res.begin(), res.end(), v1);
                    uint idx = lower - res.begin();
                    if (lower != res.end() && lower->val == v1) {}
                    else {
                        res.insert(lower, Edge (v1, vector<uint>{}));
                    }
                    auto true_csr = snaps[snap_i].neighbors_[i];
                    for (auto& j : true_csr) j = snaps[snap_i].GetVertexLabel(j);
                    std::sort(true_csr.begin(), true_csr.end());
                    std::vector<uint> tmp;
                    std::set_union(res[idx].csr.begin(), res[idx].csr.end(),
                                   true_csr.begin(), true_csr.end(),
                                   std::back_inserter(tmp));
                    for (auto j : res[idx].csr);
                    res[idx].csr.assign(tmp.begin(), tmp.end());
                }
            }
            else {
                for (int i = 0; i < snaps[snap_i].neighbors_.size(); ++i) {
                    if (snap_i == L) {
                        auto v1 = snaps[snap_i].GetVertexLabel(i);
                        auto true_csr = snaps[snap_i].neighbors_[i];
                        for (auto& j : true_csr) j = snaps[snap_i].GetVertexLabel(j);
                        res.insert(res.end(), Edge (v1, true_csr));
                        st[v1] = 0;
                    }
                    else {
                        auto v1 = snaps[snap_i].GetVertexLabel(i);
                        cout << "v1: " << v1 << endl;
                        if (st.find(v1) == st.end()) continue;
                        auto lower = std::lower_bound(res.begin(), res.end(), v1);
                        uint idx = lower - res.begin();
                        if (lower == res.end() || lower->val != v1) continue;
                        auto true_csr = snaps[snap_i].neighbors_[i];
                        for (auto& j : true_csr) j = snaps[snap_i].GetVertexLabel(j);
                        vector<uint> tmp;
                        std::set_intersection(res[idx].csr.begin(), res[idx].csr.end(),
                            true_csr.begin(), true_csr.end(),
                            std::back_inserter(tmp));
                        res[idx].csr = tmp;
                        st[v1]++;
                    }
                }
            }
        }
        for (auto i : st)
            cout << i.first << " " << i.second << endl;
        for (auto i : st) {
            if (i.second == 0) {
                auto lower = std::lower_bound(res.begin(), res.end(), i.first);
                uint idx = lower - res.begin();
                if (lower == res.end() || lower->val != i.first) continue;
                res[idx].csr.clear();
            }
        }
        for (const auto& i: res) {
            all_edges += i.csr.size();
        }

        std::cout << "----------------------------------------------Res: " << all_edges << std::endl;
//        CSV << "TR" << '\t' << L << '\t' << R << '\t' << all_edges << endl;
        all_edges = 0;
        cout << "peek memoty: " << mem::getValue() / 1024  << "Mb" << std::endl;

        if (CNT == 5) {
            CSV << Duration(query_start_time) / 1000 << '\t';
            cout << "BS" << '\t' << mem::getValue() / 1024 << "mb" << '\t' << Duration(query_start_time) / 1000 << "s" << std::endl;
            query_start_time = Get_Time();
            CNT = 0;
        }
    }
    input.close();

    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000 << "s" << '\t'<< "ALL Time: " << Duration(start) / 1000 << std::endl;
    CSV << Duration(query_time) / 1000 <<'\t' << Duration(start) / 1000 <<'\t'<< mem::getValue() / 1024 << std::endl;

    CSV.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb";
}

void TiLa_method(bool fun) {
    auto start = Get_Time();
    TiLa_Graph tila;
    tila.build(K, query_path);
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    get_index_mem();

    uint L, R;
    fstream input(input_file, ios::in);
    fstream CSV("../result.csv", ios::app);

    int CNT = 0;
    CSV << "TiLa" << '\t' << dataset << '\t' << K << '\t' << Duration(start) / 1000 << '\t' <<  mem::getValue() / 1024 << '\t';
    auto query_time = Get_Time();
    auto query_start_time = Get_Time();
    while (cin >> L >> R) {
        cout << "Query: " << L << " " << R << endl;
        CNT++;
        vector<pair<int, int>> res;
        dynamic_bitset<uint32_t> bit(R - L + 1);
        for (auto mp : tila.TiLaEdge) {
            for (const auto& ed : mp.second) {
                string s;
                to_string(ed.lifespan, s);
                s = s.substr(s.size() - 1 - R, R - L + 1);
                if (fun && s.find('1') != string::npos)res.emplace_back(mp.first, ed.nei_id);
                if (!fun && count(s.begin(), s.end(), '1') == s.size()) res.emplace_back(mp.first, ed.nei_id);
            }
        }
        // sort(res.begin(), res.end());
        // res.erase(unique(res.begin(), res.end()), res.end());
        std::cout << "----------------------------------------------Res: " << res.size() << std::endl;
        cout << "peek memoty: " << mem::getValue() / 1024  << "Mb" << std::endl;

//        CSV << "TR" << '\t' << L << '\t' << R << '\t' << all_edges << endl;
        all_edges = 0;
        if (CNT == 5) {
            CSV << Duration(query_start_time) / 1000 << '\t';
            cout << "TiLa" << '\t' << mem::getValue() / 1024 << "mb" << '\t' << Duration(query_start_time) / 1000 << "s" << std::endl;
            query_start_time = Get_Time();
            CNT = 0;
        }
    }
    input.close();

    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000 << "s" << '\t'<< "ALL Time: " << Duration(start) / 1000 << std::endl;
    CSV << Duration(query_time) / 1000 <<'\t' << Duration(start) / 1000 <<'\t'<< mem::getValue() / 1024 << std::endl;

    CSV.close();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb";
}

int this_tree_leaf_size;
inline int f(int x) { return x / 2 + this_tree_leaf_size; }
inline int ls(int x) { return 2 * f(x) - 2 * this_tree_leaf_size; }
inline int rs(int x) { return 2 * f(x) - 2 * this_tree_leaf_size + 1; }
inline int lf(int x) { return f(x) - 1; }
inline int rf(int x) { return f(x) + 1; }

int main(int argc, char *argv[])
{
    // DB* tree_db0, * tree_db1;
    // DB* bv_db0, * bv_db1;
    // DB* node_db0, * node_db1;
    // Status s;
    // Options op;
    // op.create_if_missing = true;
    // s = DB::Open(op, "../tree_db", &tree_db0);
    // if (!s.ok()) { cout << " tree_db0 open failed" << endl; return 0; }
    // s = DB::Open(op, "../tree_db", &tree_db1);
    // if (!s.ok()) { cout << " tree_db1 open failed" << endl; return 0; }
    // s = DB::Open(op, "../bv_db", &bv_db0);
    // if (!s.ok()) { cout << "bv_db open failed" << endl; return 0; }
    // s = DB::Open(op, "../node_db", &bv_db1);
    // if (!s.ok()) { cout << "node_db open failed" << endl; return 0; }
    //
    // int idx = 0;
    // s = bv_db0->Put(WriteOptions(), to_string(idx++), "111111111111");
    // s = bv_db0->Put(WriteOptions(), to_string(idx++), "1111001111110011");
    // s = bv_db0->Put(WriteOptions(), to_string(idx++), "1111111111111111");
    // idx = 0;
    // s = bv_db1->Put(WriteOptions(), to_string(idx++), "11110000001111110011");
    // s = bv_db1->Put(WriteOptions(), to_string(idx++), "111100000000001111110011");
    // s = bv_db1->Put(WriteOptions(), to_string(idx++), "111111110000001111110011");
    // idx = 0;
    // // in_tree_idx tr_idx glo_id
    // tree_db0->Put(WriteOptions(), to_string(idx++), "0");
    // tree_db0->Put(WriteOptions(), to_string(idx++), "1");
    // tree_db0->Put(WriteOptions(), to_string(idx++), "2");
    // idx = 0;
    // tree_db1->Put(WriteOptions(), to_string(idx++), "3");
    // tree_db1->Put(WriteOptions(), to_string(idx++), "4");
    // tree_db1->Put(WriteOptions(), to_string(idx++), "5");
    // idx = 0;
    // node_db0->Put(WriteOptions(), to_string(idx++), "-1 -1 2");
    // node_db0->Put(WriteOptions(), to_string(idx++), "-1 -1 2");
    // node_db0->Put(WriteOptions(), to_string(idx++), "0 1 2");
    // node_db1->Put(WriteOptions(), to_string(idx++), "-1 -1 2");
    // node_db1->Put(WriteOptions(), to_string(idx++), "-1 -1 2");
    // node_db1->Put(WriteOptions(), to_string(idx++), "0 1 2");
    //
    //     // s = tree_db->Get(ReadOptions(), key, &get_val);
    //
    // // Iterator* it = tree_db->NewIterator(rocksdb::ReadOptions());
    // // for (it->SeekToFirst(); it->Valid(); it->Next()) {
    // //     cout << it->key().ToString() << ": " << it->value().ToString() << endl;
    // //     tree_db->Delete(WriteOptions(), it->key());
    // // }
    //
    // tree_db0->Close();
    // tree_db1->Close();
    // bv_db->Close();
    // delete tree_db0;
    // delete tree_db1;
    // delete bv_db;
    //
    // return 0;
    time_t timep;time(&timep);printf("%s", ctime(&timep));
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
    CLI::App app{"App description"};
    string method = "1";
    dataset = "mo";
    level = 1;
    fun = true;
    // app.add_option("-d,--dataset", dataset, "query graph path")->required();
    // app.add_option("-m,--method", method, "method")->required();
    // app.add_option("-k,--K", K, "K");
    // app.add_option("-l,--mat", level, "level");
    // app.add_option("-e,--edgemax", EDGE_MAX, "edgemax");
    // CLI11_PARSE(app, argc, argv);

    if (dataset == "mo") K = 18, EDGE_MAX = 3e4 * 4;
    if (dataset == "mo2") K = 5, EDGE_MAX = 50 * 4;
    if (dataset == "su") K = 30, EDGE_MAX = 3e4 * 8;
    if (dataset == "wi") K = 105, EDGE_MAX = 3e4 * 8;
    if (dataset == "so") K = 92, EDGE_MAX = 3e4 * 256;
//    if (dataset == "so5m") K = 18, EDGE_MAX = 2e6;
    if (dataset == "test") K = 4, EDGE_MAX = 24;
    if (dataset == "test2") K = 8, EDGE_MAX = 200;
    if (dataset == "so2") K = 2773, EDGE_MAX = 3e4 * 256;

    query_path = "../dataset/" + dataset + "/q";
    input_file = "../dataset/" + dataset + "/input.txt";
    std::chrono::high_resolution_clock::time_point start, lstart;
    cout << query_path << " " << input_file << " " << method << " K=" << K << " " << level << endl;
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
         std::cout << "TiLa" << std::endl;
         TiLa_method(fun);
         return 0;
     }
    std::cout << "tree" << std::endl;
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;

    int cnt = 1;
    size_t cur_edges = 0;
    size_t csr_edge_count = 0;
    size_t snap_sum = 0;

    snap_idx.emplace_back(0);
    tree.resize(1);

    for (int snap_i = 0; snap_i < K; ++snap_i, ++cnt) {
        Graph snap;
        snap.LoadFromFile(query_path + std::to_string(snap_i));
        tree.back().tr.resize(cnt);
        one_index.resize(snap_i + 1);
        auto& b = tree.back().tr.back().bv;
        b.resize(snap.edge_count_);
        tree.back().tr.back().lson = -1;
        tree.back().tr.back().rson = -1;
        tree.back().tr.back().fa = -1;
        cur_edges += snap.edge_count_;

        for (int i = 0; i < snap.neighbors_.size(); ++i) {
            uint v1 = snap.GetVertexLabel(i);
            for (auto j : snap.neighbors_[i]) {
                uint v2 = snap.GetVertexLabel(j);
                if (v1 >= CSR.size()) CSR.resize(v1 + 1);

                auto lower = lb(CSR[v1], v2);
                if (lower != CSR[v1].end() && (*lower).first == v2) {
                    if ((*lower).second > b.size()) b.resize((*lower).second);
                    b.set((*lower).second - 1);
                } else {
                    CSR[v1].insert(lower, {v2, ++csr_edge_count});
                    if (csr_edge_count > b.size()) b.resize(csr_edge_count);
                    one_index.back()++;
                    b.set(csr_edge_count - 1);
                    id2edge.emplace_back(v2,v1);
                }
            }
        }

        tree.back().tr.back().glo_id = snap_i;
        // b.resize(b.size() - one_index.back());
        // cout << "one count: " << one_index.back() << " ";
        // one_index.back() += one_index[snap_i - 1];
        // cout << one_index.back() << endl;
//        cout << b << endl;
        if (cur_edges >= EDGE_MAX) {
            std::cout << "----------- Building tree ------------" << std::endl;
            cout << "tree size" << " " << tree.back().tr.size() << endl;
            snap_sum += tree.back().tr.size();
            snap_idx.emplace_back(snap_sum);
            // if (tree.size() > 1) {
            //     int len = tree.back().tr[0].bv.size();
            //     for (auto leaf : tree.back().tr) {
            //         auto b_tmp = leaf.bv;
            //         cout << leaf.bv << endl;
            //         leaf.bv >>= len;
            //         cout << leaf.bv << " " << one_index[leaf.glo_id] - one_index[tree.back().tr[0].glo_id - 1] - one_index[leaf.glo_id] + one_index[leaf.glo_id - 1] << endl;
            //         leaf.bv.resize(one_index[leaf.glo_id] - one_index[tree.back().tr[0].glo_id - 1] - one_index[leaf.glo_id] + one_index[leaf.glo_id - 1]);
            //         cout << leaf.bv << endl;
            //         leaf.bv.shrink_to_fit();
            //         b_tmp.resize(one_index[leaf.glo_id] - one_index[tree.back().tr[0].glo_id - 1]);
            //     }
            // }
            build(tree.back(), fun);

            cnt = 0;
            cur_edges = 0;
            if (snap_i + 1 != K) tree.resize(tree.size() + 1);
            std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
        }
    }
// 000000000000000000111111111111
// 000000000000001111001111110011
// 000000000000001111111111111111
// 000000000011110000001111110011
// 000000111100000000001111110011
// 000000111111110000001111110011
// 111111000000000000000000000000
    if (cur_edges) {
        std::cout << "----------- Building tree ------------" << std::endl;
        cout << "tree size" << " " << tree.back().tr.size() << endl;
        snap_sum += tree.back().tr.size();
        snap_idx.emplace_back(snap_sum);
        build(tree.back(), fun);
        cur_edges = 0;
        std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
    }

    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;
    // cout << tree[0].tr[0].bv << endl;
    // cout << tree.size() << endl;
    // cout << tree[0].tr.size() << endl;
    // cout << tree[1].tr.size() << endl;
    // fstream fs("../res.txt", ios::out);
    // fs << tree[0].tr[0].bv << endl;
    // fs << tree[0].tr[1].bv << endl;
    // fs << tree[0].tr[2].bv << endl;
    // fs << tree[0].tr[3].bv << endl;
    // fs << tree[0].tr[4].bv << endl;
    // fs << tree[0].tr[5].bv << endl;
    // fs << tree[0].tr[6].bv << endl;
    // fs << tree[0].tr[7].bv << endl;
    // fs << endl;
    // fs << tree[1].tr[0].bv << endl;
    // fs << tree[1].tr[1].bv << endl;
    // fs << tree[1].tr[2].bv << endl;
    // fs << tree[1].tr[3].bv << endl;
    // fs << tree[1].tr[4].bv << endl;
    // fs << tree[1].tr[5].bv << endl;
    // fs << tree[1].tr[6].bv << endl;
    // // fs << tree[1].tr[7].bv << endl;
    // fs << endl;
    // fs.close();

    // for (const auto& i : tree) {
    //     for (const auto& j : i.tr)
    //         cout << j.bv << endl;
    // }
    // return 0;

    std::cout << "----------- Snap Query ------------" << std::endl;
    uint L, R;
    fstream input(input_file, ios::in);
    fstream CSV("../result.csv", ios::app);

    int CNT = 0;
    CSV << "TR" << ',' << dataset << ',' << K << ',' << Duration(start) / 1000 << ',' << mem::getValue() / 1024 << ",";

    auto query_time = Get_Time();
    auto query_start_time = Get_Time();
    while (input >> L >> R) {
        cout << "Query: " << L << " " << R << endl;
        CNT++;
        uint Ltree, Rtree; // L,R: snap idx Ltree,Rtree: tree idx about snap
        auto lower1 = std::lower_bound(snap_idx.begin(), snap_idx.end(), L);
        if (*lower1 != L) Ltree = lower1 - snap_idx.begin() - 1;
        else Ltree = lower1 - snap_idx.begin();
        auto lower2 = std::lower_bound(snap_idx.begin(), snap_idx.end(), R);
        if (*lower2 != R) Rtree = lower2 - snap_idx.begin() - 1;
        else Rtree = lower2 - snap_idx.begin();
        cout << L << " " << Ltree << " " << R << " " << Rtree << endl;
        dynamic_bitset<uint32_t> lres;

        if (Ltree == Rtree) {
            auto &tr = tree[Ltree].tr;
            int lidx = L - leaf_size[Ltree].first, ridx = R - leaf_size[Rtree].first;
            if (R - lidx + 1 == snap_idx[Rtree + 1] - snap_idx[Ltree]) {
                lres = tr.back().bv;
//                cout << "quick " << tr.back().id << endl;
                goto RES;
            }
            if ((ridx & 1) && lidx + 1 == ridx) {
                if (tr[lidx].bv.size() > tr[ridx].bv.size()) {
                    auto s1 = tr[ridx].bv.size();
                    tr[ridx].bv.resize(tr[lidx].bv.size());
                    if (fun) lres = tr[lidx].bv | tr[ridx].bv;
                    else lres = tr[lidx].bv & tr[ridx].bv;
                    tr[ridx].bv.resize(s1);
                }
                else if (tr[lidx].bv.size() < tr[ridx].bv.size()) {
                    auto s1 = tr[lidx].bv.size();
                    tr[lidx].bv.resize(tr[ridx].bv.size());
                    if (fun) lres = tr[lidx].bv | tr[ridx].bv;
                    else lres = tr[lidx].bv & tr[ridx].bv;
                    tr[lidx].bv.resize(s1);
                }
                else {
                    if (fun) lres = tr[lidx].bv | tr[ridx].bv;
                    else lres = tr[lidx].bv & tr[ridx].bv;
                }
//                std::cout << "nei quick" << std::endl;
                goto RES;
            }
            uint ROOT = tr.back().id;
            uint HALF = leaf_size[Ltree].second;
            uint bias = 0;

            uint bias_bk = 0;
            if (lidx >= HALF && ridx >= HALF) bias = HALF, ROOT = tr[ROOT].rson, HALF /= 2;
            bias_bk = bias;
            while (!(lidx - bias_bk < HALF && ridx - bias_bk >= HALF)) {
                if (lidx - bias_bk < HALF) {
                    ROOT = tr[ROOT].lson;
                } else {
                    ROOT = tr[ROOT].rson;
                    bias_bk += HALF;
                }
                HALF /= 2;
            }
            //cout << ROOT << " " << HALF << endl;
            if ((R & 1) && HALF * 2 == (R - L + 1) && ((L & (L-1)) == 0)) {
                lres = tr[ROOT].bv;
//                std::cout << "quick" << std::endl;
                goto RES;
            }

            {
                uint idx = L - leaf_size[Ltree].first;
////            cout << "half leaf: " << leaf_size[Ltree].second << " " << idx  << endl;
                auto root = tr[ROOT].lson;
                auto half = HALF / 2;
                uint bias2 = 0;
                if (idx >= leaf_size[Ltree].second)
                    bias2 = leaf_size[Ltree].second;

                if (idx - bias2 >= half) {
                    bias_bk = half;
                    root = tr[root].rson;
                    half /= 2;
//                    cout << "right" << endl;
                } else {
//                    cout << "left" << endl;
                }
                if ((!(L & 1)) && R >= L + half * 2) {
                    lres = tr[root].bv;
//                    cout << "quick " << root << endl;
                    if (R - L + 1 == half * 2)
                        goto RES;
                }
                else {
                    while (idx - bias_bk < half) {
                        if (tr[tr[root].lson].id == idx || tr[tr[root].rson].id == idx) break;
                        int fl = 0;
                        size_t sz=0;
                        if (lres.size() > tr[tr[root].rson].bv.size()) {
                            fl = 1;
                            sz = tr[tr[root].rson].bv.size();
                            tr[tr[root].rson].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[root].rson].bv.size()) {
                            if (lres.empty()) lres.resize(tr[tr[root].rson].bv.size(), !fun);
                            else lres.resize(tr[tr[root].rson].bv.size());
                        }
                        if (fun) lres |= tr[tr[root].rson].bv;
                        else lres &= tr[tr[root].rson].bv;
                        if (fl == 1) tr[tr[root].rson].bv.resize(sz);
//                        cout << "cal:" << tr[tr[root].rson].id << endl;
                        root = tr[root].lson;
                        half /= 2;
//                cout << root << " " << half << endl;
                    }
                    if (idx & 1) {
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[idx].bv.size()) {
                            fl = 1;
                            sz = tr[idx].bv.size();
                            tr[idx].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[idx].bv.size()) {
                            if (lres.empty()) lres.resize(tr[idx].bv.size(), !fun);
                            else lres.resize(tr[idx].bv.size());
                        }
                        if (fun) lres |= tr[idx].bv;
                        else lres &= tr[idx].bv;

                        if (fl == 1) tr[idx].bv.resize(sz);
//                        cout << "cal:" << tr[idx].id << endl;
                    } else {
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[tr[idx].fa].bv.size()) {
                            fl = 1;
                            sz = tr[tr[idx].fa].bv.size();
                            tr[tr[idx].fa].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[idx].fa].bv.size()) {
                            if (lres.empty()) lres.resize(tr[tr[idx].fa].bv.size(), !fun);
                            else lres.resize(tr[tr[idx].fa].bv.size());
                        }
                        if (fun) lres |= tr[tr[idx].fa].bv;
                        else lres &= tr[tr[idx].fa].bv;
                        if (fl == 1) tr[tr[idx].fa].bv.resize(sz);

//                        cout << "cal:" << tr[tr[idx].fa].id << endl;
                    }
                }
            }
            {
                uint idx = R - leaf_size[Rtree].first;
////            cout << "half leaf: " << leaf_size[Rtree].second << " " << idx << endl;

                auto root = tr[ROOT].rson;
                auto half = HALF / 2;
                bias = HALF;
                bias_bk = leaf_size[Rtree].second;
                if (idx - bias >= half) {
//                    cout << "right" << endl;

                    if ((!(L & 1)) && (R & 1) && R >= L + half * 2 && (((R + 1) & R) == 0)) {
                        if (tr[root].bv.size() > lres.size()) {
                            if (lres.empty()) lres.resize(tr[root].bv.size(), !fun);
                            else lres.resize(tr[root].bv.size());
                        }
                        if (fun) lres |= tr[root].bv;
                        else lres &= tr[root].bv;
//                        cout << "quick " << root << endl;
                        goto RES;
                    }

                    bias = half;
                    if (tr[root].rson != -1) {
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[tr[root].lson].bv.size()) {
                            fl = 1;
                            sz = tr[tr[root].lson].bv.size();
                            tr[tr[root].lson].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[root].lson].bv.size()) {
                            if (lres.empty()) lres.resize(tr[tr[root].lson].bv.size(), !fun);
                            else lres.resize(tr[tr[root].lson].bv.size());
                        }

                        if (fun) lres |= tr[tr[root].lson].bv;
                        else lres &= tr[tr[root].lson].bv;
                        if (fl == 1) tr[tr[root].lson].bv.resize(sz);
//                        cout << "cal:" << tr[tr[root].lson].id << endl;
                        root = tr[root].rson;
                    }
                    else {
                        root = tr[root].lson;
                    }
                    half /= 2;
                } else {
//                    cout << "left" << endl;
                }

                bool f = true;
                while (tr[tr[root].lson].lson != -1 && half > 1) {
                    if (tr[tr[root].lson].id == idx || tr[tr[root].rson].id == idx) break;

                    while (idx - bias < half) {
                        //                lres |= tr[tr[root].lson].bv;
                        root = tr[root].lson;
                        if (idx - bias == half - 1) {
                            int fl = 0;
                            size_t sz = 0;
                            if (lres.size() > tr[root].bv.size()) {
                                fl = 1;
                                sz = tr[root].bv.size();
                                tr[root].bv.resize(lres.size());
                            }
                            else if (lres.size() < tr[root].bv.size()) {
                                if (lres.empty()) lres.resize(tr[root].bv.size(), !fun);
                                else lres.resize(tr[root].bv.size());
                            }

                            if (fun) lres |= tr[root].bv;
                            else lres &= tr[root].bv;
                            if (fl == 1) tr[root].bv.resize(sz);
//                            cout << "cal:" << tr[root].id << endl;
                            f = false;
                            break;
                        }
                        half /= 2;
////                    cout << root << " " << half << endl;
                    }
                    if (!f) break;
                    if (tr[tr[root].lson].lson == -1)
                        break;
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[tr[root].lson].bv.size()) {
                        fl = 1;
                        sz = tr[tr[root].lson].bv.size();
                        tr[tr[root].lson].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[tr[root].lson].bv.size()) {
                        if (lres.empty()) lres.resize(tr[tr[root].lson].bv.size(), !fun);
                        else lres.resize(tr[tr[root].lson].bv.size());
                    }

                    if (fun) lres |= tr[tr[root].lson].bv;
                    else lres &= tr[tr[root].lson].bv;
                    if (fl == 1)
                        tr[tr[root].lson].bv.resize(sz);
//                    cout << "cal:" << tr[tr[root].lson].id << endl;
                    if (tr[root].rson == -1) {
                        f = false;
                        break;
                    }
                    root = tr[root].rson;
                    bias += half;
                    half /= 2;
//                cout << root << " " << half << endl;
                }
                if (f) {
                    if (idx & 1) {
//                        cout << "cal:" << tr[tr[idx].fa].id << endl;
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[tr[idx].fa].bv.size()) {
                            fl = 1;
                            sz = tr[tr[idx].fa].bv.size();
                            tr[tr[idx].fa].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[idx].fa].bv.size()) {
                            if (lres.empty()) lres.resize(tr[tr[idx].fa].bv.size(), !fun);
                            else lres.resize(tr[tr[idx].fa].bv.size());
                        }

                        if (fun) lres |= tr[tr[idx].fa].bv;
                        else lres &= tr[tr[idx].fa].bv;
                        if (fl == 1) tr[tr[idx].fa].bv.resize(sz);
                    } else {
//                        cout << "cal:" << tr[idx].id << endl;
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[idx].bv.size()) {
                            fl = 1;
                            sz = tr[idx].bv.size();
                            tr[idx].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[idx].bv.size()) {
                            if (lres.empty()) lres.resize(tr[idx].bv.size(), !fun);
                            else lres.resize(tr[idx].bv.size());
                        }

                        if (fun) lres |= tr[idx].bv;
                        else lres &= tr[idx].bv;
                        if (fl == 1) tr[idx].bv.resize(sz);
                    }
                }
            }
            goto RES;
        }

        /*// 左边
        if (L != snap_idx[Ltree]) {
            auto &tr = tree[Ltree].tr;
            uint idx = L - leaf_size[Ltree].first;
////            cout << "half leaf: " << leaf_size[Ltree].second << " " << idx  << endl;
            uint bias = 0;
            uint half = leaf_size[Ltree].second;
            uint root = tr.back().id;
            if (idx >= half) {
                bias = half;
                root = tr[root].rson;
                half /= 2;
                //cout << "right" << endl;
            } else {
                //cout << "left" << endl;
            }
            bool pass = false;
            if (idx == bias) {
                int fl = 0;
                size_t sz = 0;
                if (lres.size() > tr[root].bv.size()) {
                    fl = 1;
                    sz  = tr[root].bv.size();
                    tr[root].bv.resize(lres.size());
                }
                else if (lres.size() < tr[root].bv.size()) {
                    if (lres.empty()) lres.resize(tr[root].bv.size(), !fun);
                    else lres.resize(tr[root].bv.size());
                }
                if (fun) lres |= tr[root].bv;
                else lres &= tr[root].bv;
                if (fl == 1) tr[root].bv.resize(sz);
//                cout << "cal:" << tr[root].id << endl;
                pass = true;
            }

            while (!pass && half > 1) {
                while (idx - bias < half ) {
                    if (tr[root].rson != -1) {
                        int fl = 0;
                        size_t sz = 0;
                        if (lres.size() > tr[tr[root].rson].bv.size()) {
                            fl = 1;
                            sz = tr[tr[root].rson].bv.size();
                            tr[tr[root].rson].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[root].rson].bv.size()) {
                            if (lres.empty()) lres.resize(tr[tr[root].rson].bv.size(), !fun);
                            else lres.resize(tr[tr[root].rson].bv.size());
                        }

                        if (fun) lres |= tr[tr[root].rson].bv;
                        else lres &= tr[tr[root].rson].bv;
                        if (fl == 1) tr[tr[root].rson].bv.resize(sz);
//                        cout << "cal:" << tr[tr[root].rson].id << endl;
                        root = tr[root].lson;
                        half /= 2;
//                        bias += half;
                    } else {
                        break;
                    }
                }
                if (idx > half)
                    bias += half;
                else if (idx == half) {
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[tr[root].rson].bv.size()) {
                        fl = 1;
                        sz = tr[tr[root].rson].bv.size();
                        tr[tr[root].rson].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[tr[root].rson].bv.size()) {
                        if (lres.empty()) lres.resize(tr[tr[root].rson].bv.size(), !fun);
                        else lres.resize(tr[tr[root].rson].bv.size());
                    }

                    if (fun) lres |= tr[tr[root].rson].bv;
                    else lres &= tr[tr[root].rson].bv;
                    if (fl == 1) tr[tr[root].rson].bv.resize(sz);
//                    cout << "cal:" << tr[tr[root].rson].id << endl;
                    pass = true;
                    break;
                }
            }
            if (!pass) {
                if (idx & 1) {
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[idx].bv.size()) {
                        fl = 1;
                        sz = tr[idx].bv.size();
                        tr[idx].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[idx].bv.size()) {
                        if (lres.empty()) lres.resize(tr[idx].bv.size(), !fun);
                        else lres.resize(tr[idx].bv.size());
                    }

                    if (fun) lres |= tr[idx].bv;
                    else lres &= tr[idx].bv;
                    if (fl == 1) tr[idx].bv.resize(sz);
//                    cout << "cal:" << tr[idx].id << endl;
                } else {
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[tr[idx].fa].bv.size()) {
                        fl = 1;
                        sz = tr[tr[idx].fa].bv.size();
                        tr[tr[idx].fa].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[tr[idx].fa].bv.size()) {
                        if (lres.empty()) lres.resize(tr[tr[idx].fa].bv.size(), !fun);
                        else lres.resize(tr[tr[idx].fa].bv.size());
                    }

                    if (fun) lres |= tr[tr[idx].fa].bv;
                    else lres &= tr[tr[idx].fa].bv;
                    if (fl == 1) tr[tr[idx].fa].bv.resize(sz);
//                    cout << "cal:" << tr[tr[idx].fa].id << endl;
                }
            }
        }
        else {
            int fl = 0;
            size_t sz = 0;
            if (lres.size() > tree[Ltree].tr.back().bv.size()) {
                fl = 1;
                sz = tree[Ltree].tr.back().bv.size();
                tree[Ltree].tr.back().bv.resize(lres.size());
            }
            else if (lres.size() < tree[Ltree].tr.back().bv.size()) lres.resize(tree[Ltree].tr.back().bv.size(), !fun);

            if (fun) lres |= tree[Ltree].tr.back().bv;
            else lres &= tree[Ltree].tr.back().bv;
            if (fl == 1) tree[Ltree].tr.back().bv.resize(sz);
            cout << "cal:" << tree[Ltree].tr.back().id << endl;
//            cout << "all tree" << endl;
        }*/
        if (L != snap_idx[Ltree]) {
            auto &tr = tree[Ltree].tr;
            this_tree_leaf_size = snap_idx[Ltree + 1] - snap_idx[Ltree];
            int level_size = this_tree_leaf_size / 2;
            int stop_idx = this_tree_leaf_size;

            while (1) {
                if (L & 1) {
                    if (lres.size() > tr[L].bv.size()) {
                        tr[L].bv.resize(lres.size());
                        lres |= tr[L].bv;
                    }
                    else if (lres.size() < tr[L].bv.size()) {
                        lres.resize(tr[L].bv.size());
                        lres |= tr[L].bv;
                    }
                    if (stop_idx - 1 == L)
                        break;
                    stop_idx += level_size;
                    level_size /= 2;
                    L = rf(L);
                } else {
                    while (!(L & 1)) {
                        L = f(L);
                        stop_idx += level_size;
                        level_size /= 2;
                    }
                }
            }
        }
        else {
            int fl = 0;
            size_t sz = 0;
            if (lres.size() > tree[Ltree].tr.back().bv.size()) {
                fl = 1;
                sz = tree[Ltree].tr.back().bv.size();
                tree[Ltree].tr.back().bv.resize(lres.size());
            }
            else if (lres.size() < tree[Ltree].tr.back().bv.size()) lres.resize(tree[Ltree].tr.back().bv.size(), !fun);
            if (fun) lres |= tree[Ltree].tr.back().bv;
            else lres &= tree[Ltree].tr.back().bv;
            if (fl == 1) tree[Ltree].tr.back().bv.resize(sz);
        }
        // 中间
        for (uint i = Ltree + 1; i < Rtree; ++i) {
            int fl = 0;
            size_t sz = 0;
            auto tmp = tree[i].tr.back().bv;
            if (lres.size() > tmp.size()) {
                fl = 1;
                sz = tmp.size();
                tmp.resize(lres.size());
            }
            else if (lres.size() < tmp.size()) lres.resize(tmp.size());
            if (fun) lres |= tmp;
            else lres &= tmp;
            if (fl == 1) tmp.resize(sz);
        }
        // 右边
        /*if (R != snap_idx[Rtree + 1] - 1) {
//            cout << "Right" << endl;
            auto &tr = tree[Rtree].tr;
            uint idx = R - leaf_size[Rtree].first;
            if (idx == 0) {
                int fl = 0;
                size_t sz = 0;
                if (lres.size() > tr[0].bv.size()) {
                    fl = 1;
                    sz = tr[0].bv.size();
                    tr[0].bv.resize(lres.size());
                }
                else if (lres.size() < tr[0].bv.size()) lres.resize(tr[0].bv.size());
                if (fun) lres |= tr[0].bv;
                else lres &= tr[0].bv;
                if (fl == 1) tr[0].bv.resize(sz);
//                cout << "cal:" << tr[0].id << endl;
                goto RES;
            }

////            cout << "half leaf: " << leaf_size[Rtree].second << " " << idx << endl;
            uint bias = 0;
            uint half = leaf_size[Rtree].second;
            uint root = tr.back().id;
            if (idx >= half) {
                bias = half;
                int fl = 0;
                size_t sz = 0;
                if (lres.size() > tr[tr[root].lson].bv.size()) {
                    fl = 1;
                    sz = tr[tr[root].lson].bv.size();
                    tr[tr[root].lson].bv.resize(lres.size());
                }
                else if (lres.size() < tr[tr[root].lson].bv.size()) lres.resize(tr[tr[root].lson].bv.size());

                if (fun) lres |= tr[tr[root].lson].bv;
                else lres &= tr[tr[root].lson].bv;
                if (fl == 1) tr[tr[root].lson].bv.resize(sz);
//                cout << "cal:" << tr[tr[root].lson].id << endl;
                if (tr[root].rson == -1) root = tr[root].lson;
                else root = tr[root].rson;
                half /= 2;
                //cout << "right" << endl;
            } else {
                //cout << "left" << endl;
            }
            while (tr[root].rson == -1)
                root = tr[root].lson;
            bool f = true;
            while (tr[root].lson != -1 && tr[tr[root].lson].lson != -1) {
                while (idx - bias < half) {
                    //                lres |= tr[tr[root].lson].bv;

                    root = tr[root].lson;
                    if (idx + 1 == half) {
                        int fl = 0;
                        size_t  sz = 0;
                        if (lres.size() > tr[root].bv.size()) {
                            fl = 1;
                            sz = tr[root].bv.size();
                            tr[root].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[root].bv.size()) lres.resize(tr[root].bv.size());
                        if (fun) lres |= tr[root].bv;
                        else lres &= tr[root].bv;
                        if (fl == 1) tr[root].bv.resize(sz);
//                        cout << "cal: " << tr[root].id << endl;
                        goto RES;
                    }
                    if (tr[root].lson == -1 || (tr[root].lson != -1 && tr[tr[root].lson].lson == -1)) break;
                    if (idx - bias == half - 1) {
                        int fl = 0;
                        size_t  sz = 0;
                        if (lres.size() > tr[tr[root].lson].bv.size()) {
                            fl = 1;
                            sz = tr[tr[root].lson].bv.size();
                            tr[tr[root].lson].bv.resize(lres.size());
                        }
                        else if (lres.size() < tr[tr[root].lson].bv.size()) lres.resize(tr[tr[root].lson].bv.size());

                        if (fun) lres |= tr[tr[root].lson].bv;
                        else lres &= tr[tr[root].lson].bv;
                        if (fl == 1) tr[tr[root].lson].bv.resize(sz);
//                        cout << "cal:" << tr[tr[root].lson].id << endl;
                        f = false;
                        break;
                    }
                    half /= 2;
////                    cout << root << " " << half << endl;
                }
                if (!f) break;
                if (tr[root].lson == -1 || (tr[root].lson != -1 && tr[tr[root].lson].lson == -1))
                    break;
                int fl = 0;
                size_t  sz = 0;
                if (lres.size() > tr[tr[root].lson].bv.size()) {
                    fl = 1;
                    sz = tr[tr[root].lson].bv.size();
                    tr[tr[root].lson].bv.resize(lres.size());
                }
                else if (lres.size() < tr[tr[root].lson].bv.size()) lres.resize(tr[tr[root].lson].bv.size());

                if (fun) lres |= tr[tr[root].lson].bv;
                else lres &= tr[tr[root].lson].bv;
                if (fl == 1) tr[tr[root].lson].bv.resize(sz);
//                cout << "cal:" << tr[tr[root].lson].id << endl;
                if (tr[root].rson == -1) break;
                root = tr[root].rson;
                bias += half;
                half /= 2;
////                cout << root << " " << half << endl;
            }
            if (f) {
                if ((idx - 1) & 1) {
//                    cout << "cal:" << tr[idx].id << endl;
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[idx].bv.size()) {
                        fl = 1;
                        sz = tr[idx].bv.size();
                        tr[idx].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[idx].bv.size()) lres.resize(tr[idx].bv.size());
                    if (fun) lres |= tr[idx].bv;
                    else lres &= tr[idx].bv;
                    if (fl == 1) tr[idx].bv.resize(sz);
                } else {
//                    cout << "cal:" << tr[tr[idx].fa].id << endl;
                    int fl = 0;
                    size_t sz = 0;
                    if (lres.size() > tr[tr[idx].fa].bv.size()) {
                        fl = 1;
                        sz = tr[tr[idx].fa].bv.size();
                        tr[tr[idx].fa].bv.resize(lres.size());
                    }
                    else if (lres.size() < tr[tr[idx].fa].bv.size()) lres.resize(tr[tr[idx].fa].bv.size());
                    if (fun) lres |= tr[tr[idx].fa].bv;
                    else lres &= tr[tr[idx].fa].bv;
                    if (fl == 1) tr[tr[idx].fa].bv.resize(sz);
                }
            }
        }
        else {
//            cout << "Right" << endl;
            int fl = 0;
            size_t sz = 0;
            if (lres.size() > tree[Rtree].tr.back().bv.size()) {
                fl = 1;
                sz = tree[Rtree].tr.back().bv.size();
                tree[Rtree].tr.back().bv.resize(lres.size());
            }
            else if (lres.size() < tree[Rtree].tr.back().bv.size()) lres.resize(tree[Rtree].tr.back().bv.size());

            if (fun) lres |= tree[Rtree].tr.back().bv;
            else lres &= tree[Rtree].tr.back().bv;
            if (fl == 1) tree[Rtree].tr.back().bv.resize(sz);
//            cout << "cal:" << tree[Rtree].tr.back().id << endl;
            //cout << "one leaf" << endl;
        }*/
        if (R != snap_idx[Rtree + 1] - 1) {
            auto &tr = tree[Rtree].tr;
            this_tree_leaf_size = snap_idx[Rtree + 1] - snap_idx[Rtree];
            int level_size = this_tree_leaf_size;
            int stop_idx = 0;

            while (1) {
                if (R & 1) {
                    while (R & 1) {
                        R = f(R);
                        stop_idx += level_size;
                        level_size /= 2;
                    }
                } else {
                    if (lres.size() > tr[R].bv.size()) {
                        tr[R].bv.resize(lres.size());
                        lres |= tr[R].bv;
                    }
                    else if (lres.size() < tr[R].bv.size()) {
                        lres.resize(tr[R].bv.size());
                        lres |= tr[R].bv;
                    }
                    if (stop_idx == R)
                        break;
                    stop_idx += level_size;
                    level_size /= 2;
                    R = lf(R);
                }
            }
        }
        else {
            int fl = 0;
            size_t sz = 0;
            if (lres.size() > tree[Rtree].tr.back().bv.size()) {
                fl = 1;
                sz = tree[Rtree].tr.back().bv.size();
                tree[Rtree].tr.back().bv.resize(lres.size());
            }
            else if (lres.size() < tree[Rtree].tr.back().bv.size()) lres.resize(tree[Rtree].tr.back().bv.size());

            if (fun) lres |= tree[Rtree].tr.back().bv;
            else lres &= tree[Rtree].tr.back().bv;
            if (fl == 1) tree[Rtree].tr.back().bv.resize(sz);
        }
        RES:

        auto one = lres.find_first();
//        fstream f("../restree.txt", ios::app);
        while (one != lres.npos) {
////            cout << one << endl;
//            cout << id2edge[one].first << " " << id2edge[one].second << endl;
//            f << id2edge[one].first << " " << id2edge[one].second << endl;
            one = lres.find_next(one);
            all_edges++;
        }
//        f.close();
        std::cout << "----------------------------------------------Res: " << all_edges << std::endl;
//        CSV << "TR" << ',' << L << ',' << R << ',' << all_edges << endl;
        cout << "peek memoty: " << mem::getValue() / 1024  << "Mb" << std::endl;

        all_edges = 0;
        if (CNT == 5) {
            CSV << Duration(query_start_time) / 1000 << ',';
            cout << "TR" << '\t' << mem::getValue() / 1024 << "mb" << '\t' << Duration(query_start_time) / 1000 << "s" << std::endl;
            query_start_time = Get_Time();
            CNT = 0;
        }
    }
    input.close();

    // std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000 << "s" << '\t'<< "ALL Time: " << Duration(start) / 1000 << std::endl;
    CSV << Duration(query_time) / 1000 <<',' << Duration(start) / 1000 <<','<< mem::getValue() / 1024 << ',';

    if (CHECK) {
        size_t sz = 0;
        size_t sz2 = 0;
        for (auto i : tree) {
            for (auto j : i.tr) {
                sz += j.bv.size();
                sz2 += 4;
            }
        }
        cout << "Tree memory: " << '\t'  << (long double)(sz / 8 + sz2 * 4) / 1024 / 1024 << endl;
        CSV << (long double)(sz / 8 + sz2 * 4) / 1024 / 1024 << ',';
        sz = 0, sz2 = 0;
        for (const auto& i : CSR) sz += i.size() * 2 * sizeof(unsigned int);
        cout << "CSR memory: " << '\t'  << (long double)sz / 1024 / 1024 << endl;
        CSV << (long double)sz / 1024 / 1024 << ',';
        sz = id2edge.size() * 2 * sizeof(unsigned int);
        cout << "2edge memory: " << '\t' << (long double)sz / 1024 / 1024 << endl;
        CSV << (long double)sz / 1024 / 1024 << ',';
//        cout << "now memory: " << ',' << get_index_mem()["now"] / 1024 << endl;
    }
    CSV << "base" << std::endl; CSV.close();
    return 0;
}
