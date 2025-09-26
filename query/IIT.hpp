#include <boost/dynamic_bitset.hpp>
// #include "staticcore.hpp"
// #include "CECI.hpp"
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

std::vector<pair<int, int>> id2edge;

struct Tree {
    std::vector<TNode> tr;
    int leaf;
};


std::vector<uint> snap_idx;

auto lb(std::vector<std::pair<uint, uint> > &V, uint val) {
    size_t L = 0, R = V.size();
    while (L < R) {
        size_t MID = (L + R) >> 1;
        if (V[MID].first >= val) R = MID;
        else L = MID + 1;
    }
    return V.begin() + L;
}

class IIT {
public:
    int alpha;
    int beta;
    int K;
    int depth;
    string query_path;
    string dataset;
    string input_file;
    bool fun;
    std::vector<std::vector<std::pair<uint, uint> > > CSR;
    std::vector<Tree> tree;
    bool full;
    int tree_idx;
    int cnt;
    size_t cur_edges = 0;
    size_t csr_edge_count = 0;
    size_t snap_sum = 0;
    IIT(int K_, int alpha_, int beta_, bool fun_, string query_path_, string dataset_, string input_file_) {
        K = K_;
        alpha = alpha_;
        beta = beta_;
        query_path = query_path_;
        dataset = dataset_;
        input_file = input_file_;
        fun = fun_;
        full = false;
    }

    void load() {
        cnt = 0;

        snap_idx.emplace_back(0);

        auto tmp = beta;
        while(tmp) {
            tmp /= alpha;
            depth++;
        }

        tree.resize((K + beta - 1) / beta);
        for (auto& i : tree) i.tr.resize((alpha * beta - 1) / (alpha - 1));

        tree_idx = 0;

        double t2 = 0;
        for (int snap_i = 0; snap_i < K; ++snap_i) {
            string path = query_path + std::to_string(snap_i);
            if (!io::file_exists(path.c_str()))
            {
                std::cout << "Failed to open: " << path << std::endl;
                exit(-1);
            }

            auto &b = tree[tree_idx].tr[cnt].bv;
            size_t edge_count_ = 0;
            std::ifstream ifs(path);

            uint v1, v2;
            while (ifs >> v1 >> v2)
            {
                edge_count_++;
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
            }
            ifs.close();
            cur_edges++;
            cnt++;

            if (cur_edges >= beta) {
                snap_sum += cnt;
                snap_idx.emplace_back(snap_sum);
                build(tree[tree_idx], fun);
                cnt = 0;
                cur_edges = 0;
                tree_idx++;
            }
        }
        if (cur_edges) {
            snap_sum += cnt;
            snap_idx.emplace_back(snap_sum);

            build(tree[tree_idx], fun);

            cur_edges = 0;
            if (cnt == beta) {
                tree_idx++;
                full = true;
            }
        }
    }

    double getIndexMemory() {
        size_t intSz = 0;
        for (auto& i : tree) {
            for (auto& j : i.tr) {
                intSz += 5;
                intSz += j.bv.num_blocks();
            }
        }

        return 1.0 * intSz * 4 / 1024 / 1024;
    }

    void build(Tree &tre, bool f) {
        auto &tr = tre.tr;

        size_t nums = tr.size();

        int cur_depth = 1;
        int n = beta;

        for (size_t i = 0; i < nums; ++i) tr[i].id = i, tr[i].dep = 0;
        tr[nums - 1].fa = nums - 1;
        int idx = n;

        while (cur_depth < depth) {
            int i = n - idx + 2;
            idx = 0;
            for (; i <= n; i += 2) {

                tr[n + idx].lson = i - 2;
                tr[n + idx].rson = i - 1;
                tr[n + idx].dep = cur_depth;
                tr[i - 2].fa = n + idx;
                tr[i - 1].fa = n + idx;

                if (tr[i - 2].bv.empty() || tr[i - 1].bv.empty()) {
                    idx++;
                    continue;
                }
                if (fun) {
                    if (tr[i - 2].bv.size() <= tr[i - 1].bv.size()) {
                        tr[n + idx].bv = tr[i - 2].bv;
                        tr[n + idx].bv.resize(tr[i - 1].bv.size());
                        tr[n + idx].bv |= tr[i - 1].bv;
                    }
                    else {
                        tr[n + idx].bv = tr[i - 1].bv;
                        tr[n + idx].bv.resize(tr[i - 2].bv.size());
                        tr[n + idx].bv |= tr[i - 2].bv;
                    }
                }
                else {
                        if (tr[i - 2].bv.size() <= tr[i - 1].bv.size()) {
                            tr[n + idx].bv = tr[i - 1].bv;
                            tr[n + idx].bv.resize(tr[i - 2].bv.size());
                            tr[n + idx].bv &= tr[i - 2].bv;
                        }
                        else {
                            tr[n + idx].bv = tr[i - 2].bv;
                            tr[n + idx].bv.resize(tr[i - 1].bv.size());
                            tr[n + idx].bv &= tr[i - 1].bv;
                        }
                    }

                idx++;
            }
            n += idx;
            cur_depth++;
        }
    }

    void update(int UPDATECOUNT) {
        K += UPDATECOUNT;
        if (K > beta * tree.size()) {
            tree.resize((K + beta - 1) / beta);
            for (auto& i : tree)
                if (i.tr.empty()) i.tr.resize((alpha * beta - 1) / (alpha - 1));
        }
        if (!full) {
            snap_sum -= cnt;
            snap_idx.pop_back();
        }
        full = false;
       for (int snap_i = K - UPDATECOUNT; snap_i < UPDATECOUNT; ++snap_i) {
            string path = query_path + std::to_string(snap_i);
            if (!io::file_exists(path.c_str()))
            {
                std::cout << "Failed to open: " << path << std::endl;
                exit(-1);
            }

            auto &b = tree[tree_idx].tr[cnt].bv;

            size_t edge_count_ = 0;
            std::ifstream ifs(path);

            uint v1, v2;
            while (ifs >> v1 >> v2)
            {
                edge_count_++;
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

            }
            ifs.close();

            cur_edges++;
            cnt++;
            if (cnt >= beta) {

                snap_sum += cnt;
                snap_idx.emplace_back(snap_sum);
                build(tree[tree_idx], fun);
                cnt = 0;
                cur_edges = 0;
                tree_idx++;
            }
        }
        if (cur_edges) {
            snap_sum += cnt;
            snap_idx.emplace_back(snap_sum);
            build(tree[tree_idx], fun);
            cur_edges = 0;
            if (cnt == beta) {
                tree_idx++;
                full = true;
            }
        }
    }

    void query(string moreQuery) {
        fstream input(input_file, ios::in);
        int L, R;
        cout << input_file << endl;
        while (input >> L >> R) {
        {
            cout << "Query: " << L << " " << R << endl;
            assert(L < K);
            assert(R < K);
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
            int all_edges = 0;
            auto one = lres.find_first();
            if (moreQuery.empty()) {
                while (one != lres.npos) {
                    all_edges++;
                    one = lres.find_next(one);
                }
            }
            else if (moreQuery == "DSM") {
                std::vector<std::pair<int, int>> E;
                while (one != lres.npos) {
                    E.emplace_back(id2edge[one]);
                    all_edges++;
                    one = lres.find_next(one);
                }

                if (E.size() >= 3) {
                    // string s = "CECI," + dataset + "," + to_string(L) + "," + to_string(R) + ",";
                    auto cecit = Get_Time();
                    auto ceciRes = CECI(E);
                    auto embedding_count = ceciRes[2];
                    assert(embedding_count % 6 == 0);
                    embedding_count /= 6;
                    cout << L << "," << R << ", |V| " << ceciRes[0] << ", |E| " << ceciRes[1] << ", Count " << embedding_count << " , " << Duration(cecit) << endl;;
                }
            }
            else if (moreQuery == "KCORE") {
                std::vector<std::pair<int, int>> E;
                while (one != lres.npos) {
                    E.emplace_back(id2edge[one]);
                    all_edges++;
                    one = lres.find_next(one);
                }

                int k;
                input >> k;

                cout << L << "," << R << ", " << k << ", # of KCORE: " << KCORE(E, k) << endl;
            }

            std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;

        }
    }
        input.close();
    }
};