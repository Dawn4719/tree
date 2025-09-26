#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <cassert>
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
using SnapshotId = uint32_t;

const int N = 1 << 8;

struct Page {
    std::vector<int> elems;
    // std::vector<int> loc;
    // std::vector<int> degree;
    std::vector<std::pair<int, int>> idx;
    // int Id;
    // size_t size{};
    size_t vsize;
    // size_t refcount;
    // const int n = 1 << 9;

    Page(){
        vsize = 0;
    }
    void init(int l) {
        elems.resize(N);
        for (int i = 0; i < N; i++) {
            elems[i] = l*N+i;
        }
        // loc.resize(N);
        // degree.resize(n);
        idx.resize(N);
        for (auto& i : idx) i = {-1, -1};
    }
    Page& operator=(const Page& other) {
        if (this != &other) {
            elems = other.elems;
            // loc = other.loc;
            // degree = other.degree;
            idx = other.idx;
            // Id = other.Id;
            // size = other.size;
            // vsize = other.vsize;
        }
        return *this;
    }

    bool operator<(const int& m) const {
        return elems[0] < m;
    }
};

class LAMA {
public:
    std::vector<Page> pages;
    int vertexCount{};
    int maxV;
    size_t capacity_{};
    unsigned page_exp_{};
    size_t page_size_{};

    std::vector<std::pair<int, int>> indirectionTables;
    std::vector<int> degree;
    // std::vector<int> trueD;
    std::vector<int> pre;
    std::vector<int> ls;

    LAMA(size_t capacity=0, unsigned page_exp=1)
        : capacity_(capacity), page_exp_(page_exp), page_size_(N) {
        maxV = 0;
    }

    ~LAMA() {

    }
};

struct EDGE {
    int v;
    bool is_con;
    int egtableIdx;
    int len;
    bool del;
    std::vector<bool> st;
    int stfalse;
    EDGE() {
        v = 0;
        stfalse = 0;
        is_con = false;
        len = 0;
        del = false;
        egtableIdx = -1;
        st = {};
    }

    EDGE(int l) {
        v = l;
        is_con = false;
        len = 0;
        stfalse = 0;
        del = false;
        egtableIdx = -1;
        st = {};
    }
};

struct continuation {
    int a, b;
    int len;
    std::vector<bool> st;
    continuation()=default;
    continuation(int a, int b, int len): a(a), b(b), len(len) {
        st.assign(len, false);
    }

    bool operator<(const int& m) const {
        return b < m;
    }
};

class EdgeTable {
public:
    int edgeCount;

    std::vector<EDGE> edge_table;
};

class LLAMA {
public:
    std::vector<LAMA> lms;
    std::vector<EdgeTable> egs;
    int EdgeCount{};

    int s_;
    int K;

    LLAMA(int K_) {
        K = K_;
        lms.resize(1);
        egs.resize(1);
    }
    void load(std::string query_path) {
        std::vector<std::vector<std::pair<int,int>>> full_snaps(K+0);
        for(int s=0;s<K;++s){
                std::string path = query_path + std::to_string(s);
                std::ifstream ifs(path);
                if(!ifs){ std::cerr<<"Cannot open "<<path<<std::endl; return; }
                int u,v;
                while(ifs>>u>>v) full_snaps[s].emplace_back(u,v);
                // sort(full_snaps[s].begin(), full_snaps[s].end());
                // full_snaps[s].erase(unique(full_snaps[s].begin(), full_snaps[s].end()), full_snaps[s].end());
                // cout<<"Input snapshot "<<s<<" edges="<<full_snaps[s].size()<<endl;
                if (s == 0) {
                    ADDBatch(full_snaps[s]);
                }
                else {
                    const auto &prev = full_snaps[s-1];
                    const auto &cur  = full_snaps[s];
                    std::vector<std::pair<int,int>> added, removed;
                    set_difference(cur.begin(), cur.end(), prev.begin(), prev.end(), back_inserter(added));
                    set_difference(prev.begin(), prev.end(), cur.begin(), cur.end(), back_inserter(removed));
                    lms.resize(lms.size() + 1);
                    egs.resize(egs.size() + 1);

                    lms.back().indirectionTables = lms[lms.size() - 2].indirectionTables;
                    lms.back().pages.resize(lms.back().indirectionTables.size());

                    int maxV = 0;
                    for (auto i : cur) if (i.first > maxV) maxV = i.first;

                    build(added, removed, maxV, s);
                    full_snaps[s-1].clear();
                    full_snaps[s-1].shrink_to_fit();
                }
            }
        full_snaps.back().clear();
        full_snaps.back().shrink_to_fit();
    }

    void query(std::string input_file, bool fun) {
        std::fstream input(input_file, std::ios::in);
        int L, R;
        while (input >> L >> R) {
            std::cout << L << " " << R << std::endl;
            if (fun == true)
                query(L, R, true);
            else
                query(L, R, false);
        }
    }

    ~LLAMA() = default;
    void ADDBatch(std::vector<std::pair<int, int>>& E) {
        auto& lm = lms.back();
        int maxV = 0;
        for (auto [u, v] : E) {
            if (u > maxV) maxV = u;
            if (u >= lm.degree.size()) {
                lm.degree.resize(u + 1);
            }
        }
        lm.maxV = maxV;
        lm.vertexCount = lm.degree.size();
        lm.pre.resize(lm.degree.size() + 1);

        lm.ls.resize(lm.vertexCount, -1);

        // lm.trueD.resize(lm.vertexCount);

        int vCount = 0;
        for (int i = 0; i < lm.vertexCount; i++) {lm.ls[i] = i;}

        for (auto [u, v] : E) {
            lm.degree[lm.ls[u]]++;
        }


        for (int i = 1; i <= lm.vertexCount; i++) {
            lm.pre[i] = lm.pre[i - 1] + lm.degree[i - 1];
        }

        lm.pages.resize((lm.vertexCount + lm.page_size_ - 1) / lm.page_size_);
        for (int i = 0; i < lm.pages.size(); i++) lm.pages[i].init(i);
        lm.indirectionTables.resize(lm.pages.size());
        lm.indirectionTables[0].first = lms.size() - 1;

        for (int i = 0, j = 0, idx = 0; i < lm.vertexCount; i++) {
            if (j == lm.page_size_) {
                j = 0;
                idx++;
                // lm.pages[idx].Id = idx;
                // if (lm.pages[idx].elems[j]) {
                    lm.indirectionTables[idx].first = lms.size() - 1;
                    lm.indirectionTables[idx].second = idx;
                // }
            }
            lm.pages[idx].elems[j] = i;
            if (lm.degree[lm.pages[idx].elems[j]]) {
                lm.pages[idx].idx[j].first = lms.size() - 1;
                lm.pages[idx].idx[j].second = lm.pre[i];
            }
            lm.pages[idx].vsize++;

            j++;
        }
        if (lm.pages.back().vsize < N) {
            for (int i = lm.pages.back().vsize, j = 1; i < N; i++, j++) {
                lm.pages.back().elems[i] = lm.pages.back().elems[lm.pages.back().vsize - 1] + j;
            }
        }
        lm.degree.assign(lm.vertexCount, 0);

        auto& eg = egs.back();
        eg.edgeCount = lm.pre.back();
        eg.edge_table.resize(eg.edgeCount);

        for (auto [u, v] : E) {
            eg.edge_table[lm.pre[u] + lm.degree[u]++] = {v};
        }
        // lm.trueD = lm.degree;
    }

    std::map<int, int> tmp;

    void getReadd(int u, int u2snapidx, int ls_u, std::vector<int>& adj, std::map<int, std::vector<int>>& readd, bool& haveReadd) {
        int conSize = 0;
        for (int l = lms[u2snapidx].pre[ls_u + 1] - 1; l >= lms[u2snapidx].pre[ls_u]; l--) {
            if (egs[u2snapidx].edge_table[l].is_con) conSize++;
            else break;
        }

        for (int i = lms[u2snapidx].pre[ls_u], idx = 0; i < lms[u2snapidx].pre[ls_u + 1] - conSize; i++, idx++) {
            auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
            if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {

            }
            else
            {
                readd[u].emplace_back(egs[u2snapidx].edge_table[i].v);
                haveReadd = true;
            }
        }

        for (int i = lms[u2snapidx].pre[ls_u + 1] - conSize; i < lms[u2snapidx].pre[ls_u + 1]; i++) {
            auto eg = egs[u2snapidx].edge_table[i];
            auto& st = egs[u2snapidx].edge_table[i].st;
            assert(eg.is_con);
            for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                if (egs[eg.egtableIdx].edge_table[l].is_con) {
                    auto& eg2 = egs[eg.egtableIdx].edge_table[l];
                    auto ls_u2 = lms[eg2.egtableIdx].ls[u];
                    getReadd(u, eg2.egtableIdx, ls_u2, adj, readd, haveReadd);
                    continue;
                }
                auto lb = lower_bound(adj.begin(), adj.end(), egs[eg.egtableIdx].edge_table[l].v);
                if (st.empty() || !st[idx])
                    if (lb == adj.end() || *lb != egs[eg.egtableIdx].edge_table[l].v) {
                        readd[u].emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                        haveReadd = true;
                    }
            }
        }
    }

    void getAdj(int u, int u2snapidx, int ls_u, std::vector<int>& adj) {
        int conSize = 0;
        for (int l = lms[u2snapidx].pre[ls_u + 1] - 1; l >= lms[u2snapidx].pre[ls_u]; l--) {
            if (egs[u2snapidx].edge_table[l].is_con) conSize++;
            else break;
        }

        for (int i = lms[u2snapidx].pre[ls_u], idx = 0; i < lms[u2snapidx].pre[ls_u + 1] - conSize; i++, idx++) {
            auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
            if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {

            }
            else
            {
                adj.emplace_back(egs[u2snapidx].edge_table[i].v);
            }
        }

        for (int i = lms[u2snapidx].pre[ls_u + 1] - conSize; i < lms[u2snapidx].pre[ls_u + 1]; i++) {
            auto eg = egs[u2snapidx].edge_table[i];
            auto& st = egs[u2snapidx].edge_table[i].st;
            assert(eg.is_con);
            for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                if (egs[eg.egtableIdx].edge_table[l].is_con) {
                    auto& eg2 = egs[eg.egtableIdx].edge_table[l];
                    auto ls_u2 = lms[eg2.egtableIdx].ls[u];
                    getAdj(u, eg2.egtableIdx, ls_u2, adj);
                    continue;
                }
                auto lb = lower_bound(adj.begin(), adj.end(), egs[eg.egtableIdx].edge_table[l].v);
                if (st.empty() || !st[idx])
                    if (lb == adj.end() || *lb != egs[eg.egtableIdx].edge_table[l].v) {
                        adj.emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                    }
            }
        }
    }

    void build(std::vector<std::pair<int, int>>& ADD, std::vector<std::pair<int, int>>& DEL, int maxV, int s) {
        int la = lms.size() - 2;
        s_ = s;
        lms.back().indirectionTables = lms[la].indirectionTables;

        std::vector<int> stADD;
        std::vector<int> stDEL;

        for (auto [u, v] : ADD) stADD.emplace_back(u);
        for (auto [u, v] : DEL) stDEL.emplace_back(u);
        sort(stADD.begin(), stADD.end());
        sort(stDEL.begin(), stDEL.end());
        stADD.erase(std::unique(stADD.begin(), stADD.end()), stADD.end());
        stDEL.erase(std::unique(stDEL.begin(), stDEL.end()), stDEL.end());

        for (int indirectionTablesIdx = 0; indirectionTablesIdx < lms.back().indirectionTables.size(); indirectionTablesIdx++) {
            auto [snapid, pgid] = lms.back().indirectionTables[indirectionTablesIdx];
            if (snapid == -1) continue;
            lms.back().pages[pgid] = lms[snapid].pages[pgid];
        }

        if (maxV >= lms.back().indirectionTables.size()) {
            lms.back().indirectionTables.resize((maxV + lms.back().page_size_) / lms.back().page_size_, {-1, -1});
            lms.back().pages.resize(lms.back().indirectionTables.size());
        }


        int itADD = 0;
        int itDEL = 0;
        int DELidx = 0;
        int ADDidx = 0;
        int indirectionTablesIdx = 0;
        int vCount = 0;
        std::map<int, EDGE> deledge;
        std::map<int, std::vector<int>> readd;

        for (; indirectionTablesIdx < lms.back().indirectionTables.size() && itADD < ADD.size() && itDEL < DEL.size(); indirectionTablesIdx++) {
            auto& [snapid, pgid] = lms.back().indirectionTables[indirectionTablesIdx];
            // std::cout << indirectionTablesIdx << std::endl;

            if (N * (indirectionTablesIdx + 1) - 1 < stADD[itADD] && N * (indirectionTablesIdx + 1) - 1 < stDEL[itDEL]) {
                if (snapid == -1) continue;
                lms.back().pages[pgid] = lms[snapid].pages[pgid];
                continue;
            }

            std::vector<int> adjADD;
            std::vector<int> adjDEL;
            int bk = snapid;
            if (snapid == -1) {
                snapid = lms.size() - 1;
                pgid = indirectionTablesIdx;
                lms.back().pages[pgid].init(indirectionTablesIdx);
            }

            if (N * (pgid + 1) >= stDEL[itDEL]) {
                for (int j = 0; j < N; j++) {
                    int u = N * pgid + j;
                    if (s == 278 && u == 586581)
                        s = 278;
                    if (s == 279 && u == 586581)
                        s = 279;
                    if (s == 280 && u == 586581)
                        s = 280;
                    auto lb = lower_bound(stDEL.begin(), stDEL.end(), u);

                    // 这个点在删除边中
                    if (lb < stDEL.end() && (*lb) == u) {
                        auto [u2snapidx, u2idx] = lms[snapid].pages[pgid].idx[j];
                        auto ls_u = lms[u2snapidx].ls[u];

                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }

                        std::vector<int> adj;
                        for (; DELidx < DEL.size(); DELidx++) { if (DEL[DELidx].first != u) break; adj.emplace_back(DEL[DELidx].second); }
                        assert(lms[u2snapidx].pre.size() > ls_u + 1);
                        assert(ls_u + 1 < lms[u2snapidx].pre.size());
                        // 这个点有con链，处理con链
                        if (egs[u2snapidx].edge_table[lms[u2snapidx].pre[ls_u + 1] - 1].is_con) {


                            // 点的边删了一部分
                            // EDGE tmp;
                            // tmp.is_con = true;
                            // tmp.v = lms[u2snapidx].pre[ls_u];
                            // tmp.egtableIdx = u2snapidx;
                            // tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                            // tmp.st.resize(lms[u2snapidx].pre[ls_u + 1] - conSize - lms[u2snapidx].pre[ls_u]);
                            bool haveReadd = false;
                            getReadd(u, u2snapidx, ls_u, adj, readd, haveReadd);

                            // if (conSize != 1 && lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u] > 1)
                            //     deledge[u] = std::move(tmp);
                            if (haveReadd)
                                lms.back().pages[pgid].idx[j] = {lms.size() - 1, -1};
                            else
                                lms.back().pages[pgid].idx[j] = { - 1, -1};
                        }
                        //这个点无con链
                        else {
                            // 点的边全删了
                            assert(ls_u + 1 < lms[u2snapidx].pre.size());
                            if (adj.size() == lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u]) {
                                lms.back().pages[pgid].idx[j] = {-1, -1};
                            }
                            else {
                                // 点的边删了一部分
                                EDGE tmp;
                                tmp.is_con = true;

                                tmp.v = lms[u2snapidx].pre[ls_u];
                                tmp.egtableIdx = u2snapidx;

                                assert(ls_u != -1);
                                assert(ls_u + 1 < lms[u2snapidx].pre.size());
                                tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                                tmp.st.resize(lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u]);
                                for (int i = lms[u2snapidx].pre[ls_u], idx = 0; i < lms[u2snapidx].pre[ls_u + 1]; i++, idx++) {
                                    auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
                                    if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {
                                        // egs[u2snapidx].edge_table[i].del = true;
                                        tmp.st[idx] = true;
                                        tmp.stfalse++;
                                    }
                                }
                                tmp.stfalse = tmp.st.size() - tmp.stfalse;
                                deledge[u] = std::move(tmp);
                            }
                        }
                    }
                    // 这个点不在删除边中
                    else {
                        auto [u2snapidx, u2idx] = lms[snapid].pages[pgid].idx[j];
                        if (u2snapidx == -1 || u >= lms[u2snapidx].ls.size()) continue; // 这个点还未存在
                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        auto ls_u = lms[u2snapidx].ls[u];
                        if (ls_u == -1 || lms[u2snapidx].degree[ls_u] == 0) continue;
                        EDGE tmp;
                        tmp.is_con = true;
                        tmp.v = lms[u2snapidx].pre[ls_u];
                        tmp.egtableIdx = u2snapidx;
                        tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                        deledge[u] = std::move(tmp);
                    }
                }
            }

            if (N * (pgid + 1) >= stADD[itADD]) {
                for (int j = 0; j < N; j++) {
                    int u = N * pgid + j;
                    auto lb = lower_bound(stADD.begin(), stADD.end(), u);
                    if (s == 278 && u == 586581)
                        s = 278;
                    // 这个点在插入边中
                    if (lb < stADD.end() && (*lb) == u) {

                        assert(pgid < lms.back().pages.size());
                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        std::vector<int> adj;
                        for (; ADDidx < ADD.size(); ADDidx++) {if (ADD[ADDidx].first != u) break; adj.emplace_back(ADD[ADDidx].second);}
                        // if (readd.find(u) != readd.end()) {
                        //     for (auto l : readd[u]) adj.emplace_back(l);
                        //     sort(adj.begin(), adj.end());
                        //     adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                        //     readd.erase(u);
                        // }
                        if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                        if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                        assert(u <  lms.back().ls.size());
                        auto ls_u = lms.back().ls[u];
                        assert(ls_u != -1);
                        if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                        // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                        if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                        lms.back().degree[ls_u] = adj.size();
                        // lms.back().trueD[ls_u] = adj.size();

                        // 这个点存在旧信息
                        EDGE tmp;
                        tmp.v = -1;
                        if (lms.back().pages[pgid].idx[j].first != -1) {
                            if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                            if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                            auto ls_u = lms.back().ls[u];
                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);
                            // 在删除时已经处理了这个点的con链
                            if (!readd.empty() && readd.find(u) != readd.end()) {
                                for (auto l : readd[u]) {
                                    adj.emplace_back(l);
                                }
                                readd.erase(u);
                            }
                            else {// 在删除时未处理这个点的con链
                                // 这个点之前存在con，处理con链
                                auto [u2snapidx, u2Idx] = lms[snapid].pages[pgid].idx[j];
                                auto u2snapls_u = lms[u2snapidx].ls[u];

                                if (egs[u2snapidx].edge_table[lms[u2snapidx].pre[lms[u2snapidx].ls[u] + 1] - 1].is_con) {
                                    /*int conSize = 0;
                                    for (int l = lms[u2snapidx].pre[lms[u2snapidx].ls[u] + 1] - 1; l >= lms[u2snapidx].pre[lms[u2snapidx].ls[u]]; l--) {
                                        if (egs[u2snapidx].edge_table[l].is_con) conSize++;
                                        else break;
                                    }

                                    for (int i = lms[u2snapidx].pre[u2snapls_u], idx = 0; i < lms[u2snapidx].pre[u2snapls_u + 1] - conSize; i++, idx++) {
                                        auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
                                        if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {

                                        }
                                        else
                                        {
                                            adj.emplace_back(egs[u2snapidx].edge_table[i].v);
                                        }
                                    }

                                    for (int i = lms[u2snapidx].pre[u2snapls_u + 1] - conSize; i < lms[u2snapidx].pre[u2snapls_u + 1]; i++) {
                                        auto eg = egs[u2snapidx].edge_table[i];
                                        auto& st = egs[u2snapidx].edge_table[i].st;
                                        assert(eg.is_con);
                                        for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                                            // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                                            if (egs[eg.egtableIdx].edge_table[l].is_con) {
                                                // std::cout << "ass453" << std::endl;
                                                continue;
                                            }
                                            auto lb = lower_bound(adj.begin(), adj.end(), egs[eg.egtableIdx].edge_table[l].v);
                                            if (st.empty() || !st[idx])
                                                if (lb == adj.end() || *lb != egs[eg.egtableIdx].edge_table[l].v) {
                                                    adj.emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                                            }
                                        }
                                    }
                                    */

                                    getAdj(u, u2snapidx, u2snapls_u, adj);

                                    sort(adj.begin(), adj.end());
                                    adj.erase(unique(adj.begin(), adj.end()), adj.end());
                                    if (!deledge.empty() && deledge.find(u) != deledge.end())
                                        deledge.erase(u);
                                }
                                else {
                                    // 这个点之前不存在con，创建一个con
                                    // assert(deledge.find(u) != deledge.end());
                                    if (deledge.find(u) != deledge.end()) {
                                        tmp = deledge[u];
                                        deledge.erase(u);

                                        if (tmp.len == 1) {
                                            tmp.is_con = false;
                                            adj.emplace_back(egs[tmp.egtableIdx].edge_table[tmp.v].v);
                                            tmp.v = -1;
                                        }
                                    }
                                    else {
                                        if (u2snapls_u == -1 || lms[u2snapidx].degree[u2snapls_u] == 0) continue;
                                        if (tmp.len == 1) {
                                            tmp.is_con = false;
                                            adj.emplace_back(egs[tmp.egtableIdx].edge_table[lms[u2snapidx].pre[u2snapls_u]].v);
                                            tmp.v = -1;
                                        }
                                        else {
                                            tmp.is_con = true;
                                            tmp.v = lms[u2snapidx].pre[u2snapls_u];
                                            tmp.egtableIdx = u2snapidx;
                                            tmp.len = lms[u2snapidx].pre[u2snapls_u + 1] - lms[u2snapidx].pre[u2snapls_u];
                                        }
                                    }
                                }
                            }
                        }

                        sort(adj.begin(), adj.end());
                        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());

                        lms.back().degree[ls_u] = adj.size();
                        for (auto l : adj) egs.back().edge_table.emplace_back(l);
                        if (tmp.v != -1) {
                            assert(tmp.is_con);
                            egs.back().edge_table.emplace_back(std::move(tmp));
                            lms.back().degree[ls_u]++;
                        }
                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                        lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                    }
                    // 这个点不在插入边中
                    else {
                        bool hasadd = false;
                        if (readd.find(u) != readd.end()) {
                            sort(readd[u].begin(), readd[u].end());
                            if (u >= lms.back().ls.size()) {
                                lms.back().ls.resize(u + 1, -1);
                                lms.back().ls[u] = vCount++;
                            }
                            auto ls_u = lms.back().ls[u];
                            if (lms.back().pages[pgid].elems.empty()) {
                                lms.back().pages[pgid].init(pgid);
                            }

                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                            lms.back().degree[ls_u] = readd[u].size();
                            // lms.back().trueD[ls_u] = readd[u].size();
                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                            for (auto l : readd[u]) egs.back().edge_table.emplace_back(l);
                            readd.erase(u);
                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                            hasadd = true;
                        }
                        if (deledge.find(u) != deledge.end()) {
                            if (u >= lms.back().ls.size()) {
                                lms.back().ls.resize(u + 1, -1);
                                lms.back().ls[u] = vCount++;
                            }
                            auto ls_u = lms.back().ls[u];

                            if (lms.back().pages[pgid].elems.empty()) {
                                lms.back().pages[pgid].init(pgid);
                            }

                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                            egs.back().edge_table.emplace_back(deledge[u]);
                            deledge.erase(u);
                            assert(ls_u >= 0);
                            assert(ls_u < lms.back().degree.size());
                            lms.back().degree[ls_u]++;

                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                            hasadd = true;
                        }
                        if (!hasadd) {
                            // 这个点不在插入边中 但是所在page涉及修改 所以需要维护该点在当前快照当前page中的信息
                            if (!lms[snapid].pages[pgid].idx.empty()) {
                                auto [u2snapidx, u2Idx] = lms[snapid].pages[pgid].idx[j];
                                if (u2snapidx != -1 && ! lms.back().pages[pgid].idx.empty() && lms.back().pages[pgid].idx[j].first != -1) {
                                    auto u2snapls_u = lms[u2snapidx].ls[u];

                                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                                    auto ls_u = lms.back().ls[u];

                                    if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                                    // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                                    if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                                    lms.back().degree[ls_u] = lms[u2snapidx].degree[u2snapls_u];

                                    for (int l = u2Idx; l < u2Idx + lms[u2snapidx].degree[u2snapls_u]; l++) {
                                        egs.back().edge_table.emplace_back(egs[u2snapidx].edge_table[l]);
                                    }
                                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                    lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                                }
                            }
                        }
                        /*else {
                            if (lms.back().pages[pgid].idx[j].first != -1) {
                                if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                                if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                                auto ls_u = lms.back().ls[u];
                                if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                                // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                                if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                                if (egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u] + 1] - 1].is_con) {
                                    int conSize = 0;
                                    while (egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u] + 1] - 1 - conSize].is_con) conSize++;
                                    // std::vector<int> adj;

                                    for (int i = lms[snapid].pre[lms[snapid].ls[u] + 1] - conSize; i < lms[snapid].pre[lms[snapid].ls[u] + 1]; i++) {
                                        auto eg = egs[snapid].edge_table[i];
                                        auto& st = egs[snapid].edge_table[i].st;
                                        // assert(eg.is_con);
                                        for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                                            // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                                            if (st.empty() || !st[idx]) {
                                                lms.back().degree[ls_u]++;
                                                // lms.back().trueD[ls_u]++;
                                                egs[snapid].edge_table.emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                                            }
                                        }
                                    }

                                    if (!egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].is_con) {
                                        EDGE tmp;
                                        tmp.is_con = true;
                                        tmp.v = lms[snapid].pre[lms[snapid].ls[u]];
                                        tmp.egtableIdx = snapid;
                                        tmp.len = lms[snapid].pre[lms[snapid].ls[u] + 1] - lms[snapid].pre[lms[snapid].ls[u]];
                                        lms.back().degree[ls_u]++;
                                        // lms.back().trueD[ls_u]++;
                                        egs.back().edge_table.emplace_back(std::move(tmp));
                                    }
                                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                }
                                else {

                                    // 把上一个快照中的相同点的边作为con存入这个快照
                                    EDGE tmp;

                                    tmp.len = lms[snapid].pre[lms[snapid].ls[u] + 1] - lms[snapid].pre[lms[snapid].ls[u]];
                                    assert(tmp.len > 0);
                                    if (tmp.len != 1) {
                                        tmp.is_con = true;
                                        tmp.v = lms[snapid].pre[lms[snapid].ls[u]];
                                        tmp.egtableIdx = snapid;
                                        lms.back().degree[ls_u]++;
                                        // lms.back().trueD[ls_u]++;
                                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                        egs.back().edge_table.emplace_back(std::move(tmp));

                                        tmp.st.resize(tmp.len);
                                        bool dont = true;

                                        for (int l = lms[snapid].pre[lms[snapid].ls[u]], idx = 0; l < lms[snapid].pre[lms[snapid].ls[u] + 1]; l++, idx++) {
                                            if (egs[snapid].edge_table[l].del) {
                                                dont = false;
                                                tmp.st[idx] = true;
                                            }
                                        }
                                        if (!dont) {
                                            tmp.st.clear();
                                            tmp.st.shrink_to_fit();
                                        }
                                    }
                                    else {
                                        tmp.is_con = false;
                                        assert(u < lms[snapid].ls.size());

                                        tmp.v = egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].v;
                                        if (!egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].del) {
                                            tmp.egtableIdx = egs.size() - 1;
                                            egs.back().edge_table.emplace_back(std::move(tmp));
                                            lms.back().degree[ls_u]++;
                                            // lms.back().trueD[ls_u]++;
                                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                                            if (lms.back().pages[pgid].elems.empty()) lms.back().pages[pgid].init(pgid);

                                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                                        }
                                    }

                                }
                            }
                        }*/
                    }
                }
            }
            else {
                auto p = deledge.begin();
                for (; p != deledge.end(); ) {
                    auto u = p->first;
                    if (u < N * (pgid + 1)) {
                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        if (u >= lms.back().ls.size()) {
                            lms.back().ls.resize(u + 1, -1);
                            lms.back().ls[u] = vCount++;
                        }
                        auto ls_u = lms.back().ls[u];
                        if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                        // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                        if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                        egs.back().edge_table.emplace_back(deledge[u]);
                        // deledge.erase(u);
                        lms.back().degree[ls_u]++;

                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                        auto lb = lower_bound(lms.back().pages[pgid].elems.begin(), lms.back().pages[pgid].elems.end(), u) - lms.back().pages[pgid].elems.begin();
                        lms.back().pages[pgid].idx[lb] = {lms.size() - 1, lms.back().pre[ls_u]};
                        p = deledge.erase(p);
                    }
                    else {
                        p++;
                    }
                }

                if (!readd.empty()) {
                    auto u_ = readd.begin();
                    for (; u_ != readd.end(); ) {
                        auto u = u_->first;
                        auto j = std::lower_bound(lms.back().pages[pgid].elems.begin(), lms.back().pages[pgid].elems.end(), u) - lms.back().pages[pgid].elems.begin();
                        auto lb = lower_bound(stADD.begin(), stADD.end(), u);
                        // 这个点在插入边中

                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        std::vector<int> adj;
                        if (lb < stADD.end() && (*lb) == u) {
                            for (; ADDidx < ADD.size(); ADDidx++) {if (ADD[ADDidx].first != u) break; adj.emplace_back(ADD[ADDidx].second);}
                        }
                        if (readd.find(u) != readd.end()) {
                            for (auto l : readd[u]) adj.emplace_back(l);
                            sort(adj.begin(), adj.end());
                            adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                            u_ = readd.erase(u_);
                        }
                        else {
                            u_++;
                        }
                        if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                        if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                        auto ls_u = lms.back().ls[u];
                        if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                        // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                        if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                        lms.back().degree[ls_u] = adj.size();
                        // lms.back().trueD[ls_u] = adj.size();

                        for (auto l : adj) egs.back().edge_table.emplace_back(l);

                        if (deledge.find(u) != deledge.end()) {
                            egs.back().edge_table.emplace_back(deledge[u]);
                            deledge.erase(u);
                            lms.back().degree[ls_u]++;
                        }

                        lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                    }
                    // 这个点不在插入边中

                        // if (u >= lms.back().ls.size()) {
                        //     lms.back().ls.resize(u + 1, -1);
                        // }
                        // if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                        //
                        // assert(readd.find(u) != readd.end());
                        // auto ls_u = lms.back().ls[u];
                        //
                        //
                        //
                        // if (deledge.find(u) != deledge.end()) {
                        //     if (lms.back().pages[pgid].elems.empty()) {
                        //         lms.back().pages[pgid].init(pgid);
                        //     }
                        //
                        //     if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                        //     if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                        //     if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);
                        //
                        //     egs.back().edge_table.emplace_back(deledge[u]);
                        //     deledge.erase(u);
                        //     lms.back().degree[ls_u]++;
                        //
                        //     lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                        //
                        //     lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                        // }


                }
            }

            snapid = lms.size() - 1;

            while (itDEL < DEL.size() && stDEL[itDEL] < N * (pgid + 1)) itDEL++;
            while (itADD < ADD.size() && stADD[itADD] < N * (pgid + 1)) itADD++;
        }
        for (; indirectionTablesIdx < lms.back().indirectionTables.size() && itDEL < DEL.size(); indirectionTablesIdx++) {
            auto [snapid, pgid] = lms.back().indirectionTables[indirectionTablesIdx];
            // std::cout << "del: " << indirectionTablesIdx << std::endl;
            // if (s == 23 && indirectionTablesIdx == 167)
            //     s = 23;
            if (snapid == -1) {
                snapid = lms.size() - 1;
                pgid = indirectionTablesIdx;
            }
            std::vector<int> adjDEL;

            if (N * (pgid + 1) - 1 < stDEL[itDEL]) continue;

            if (N * (pgid + 1) >= stDEL[itDEL]) {
                if (lms.back().pages[pgid].elems.empty())
                    lms.back().pages[pgid].init(pgid);
                if (N * (pgid + 1) >= stDEL[itDEL]) {
                    for (int j = 0; j < N; j++) {
                    int u = N * pgid + j;
                    if (s == 278 && u == 586581)
                        s = 278;
                    if (s == 279 && u == 586581)
                        s = 279;
                    if (s == 280 && u == 586581)
                        s = 280;
                    auto lb = lower_bound(stDEL.begin(), stDEL.end(), u);

                    // 这个点在删除边中
                    if (lb < stDEL.end() && (*lb) == u) {
                        auto [u2snapidx, u2idx] = lms[snapid].pages[pgid].idx[j];
                        auto ls_u = lms[u2snapidx].ls[u];

                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }

                        std::vector<int> adj;
                        for (; DELidx < DEL.size(); DELidx++) { if (DEL[DELidx].first != u) break; adj.emplace_back(DEL[DELidx].second); }
                        assert(lms[u2snapidx].pre.size() > ls_u + 1);
                        assert(ls_u + 1 < lms[u2snapidx].pre.size());
                        // 这个点有con链，处理con链
                        if (egs[u2snapidx].edge_table[lms[u2snapidx].pre[ls_u + 1] - 1].is_con) {


                            // 点的边删了一部分
                            // EDGE tmp;
                            // tmp.is_con = true;
                            // tmp.v = lms[u2snapidx].pre[ls_u];
                            // tmp.egtableIdx = u2snapidx;
                            // tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                            // tmp.st.resize(lms[u2snapidx].pre[ls_u + 1] - conSize - lms[u2snapidx].pre[ls_u]);
                            bool haveReadd = false;
                            getReadd(u, u2snapidx, ls_u, adj, readd, haveReadd);

                            // if (conSize != 1 && lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u] > 1)
                            //     deledge[u] = std::move(tmp);
                            if (haveReadd)
                                lms.back().pages[pgid].idx[j] = {lms.size() - 1, -1};
                            else
                                lms.back().pages[pgid].idx[j] = { - 1, -1};
                        }
                        //这个点无con链
                        else {
                            // 点的边全删了
                            assert(ls_u + 1 < lms[u2snapidx].pre.size());
                            if (adj.size() == lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u]) {
                                lms.back().pages[pgid].idx[j] = {-1, -1};
                            }
                            else {
                                // 点的边删了一部分
                                EDGE tmp;
                                tmp.is_con = true;

                                tmp.v = lms[u2snapidx].pre[ls_u];
                                tmp.egtableIdx = u2snapidx;

                                assert(ls_u != -1);
                                assert(ls_u + 1 < lms[u2snapidx].pre.size());
                                tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                                tmp.st.resize(lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u]);
                                for (int i = lms[u2snapidx].pre[ls_u], idx = 0; i < lms[u2snapidx].pre[ls_u + 1]; i++, idx++) {
                                    auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
                                    if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {
                                        // egs[u2snapidx].edge_table[i].del = true;
                                        tmp.st[idx] = true;
                                        tmp.stfalse++;
                                    }
                                }
                                tmp.stfalse = tmp.st.size() - tmp.stfalse;
                                deledge[u] = std::move(tmp);
                            }
                        }
                    }
                    // 这个点不在删除边中
                    else {
                        auto [u2snapidx, u2idx] = lms[snapid].pages[pgid].idx[j];
                        if (u2snapidx == -1 || u >= lms[u2snapidx].ls.size()) continue; // 这个点还未存在
                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        auto ls_u = lms[u2snapidx].ls[u];
                        if (ls_u == -1 || lms[u2snapidx].degree[ls_u] == 0) continue;
                        EDGE tmp;
                        tmp.is_con = true;
                        tmp.v = lms[u2snapidx].pre[ls_u];
                        tmp.egtableIdx = u2snapidx;
                        tmp.len = lms[u2snapidx].pre[ls_u + 1] - lms[u2snapidx].pre[ls_u];
                        deledge[u] = std::move(tmp);
                    }
                    }
                }
                lms.back().indirectionTables[indirectionTablesIdx] = {lms.size() - 1, indirectionTablesIdx};
            }

            while (itDEL < DEL.size() && stDEL[itDEL] < N * (pgid + 1)) itDEL++;
        }
        for (; indirectionTablesIdx < lms.back().indirectionTables.size() && itADD < ADD.size(); indirectionTablesIdx++) {
            auto [snapid, pgid] = lms.back().indirectionTables[indirectionTablesIdx];
            // std::cout << "add: " << indirectionTablesIdx << std::endl;
            if (s == 33 && indirectionTablesIdx == 167)
                s = 33;
            if (snapid == -1) {
                snapid = lms.size() - 1;
                pgid = indirectionTablesIdx;
            }

            std::vector<int> adjADD;

            if (N * (pgid + 1) - 1 < stADD[itADD]) continue;

            if (N * (pgid + 1) >= stADD[itADD]) {
                if (lms.back().pages[pgid].elems.empty())
                    lms.back().pages[pgid].init(pgid);

                for (int j = 0; j < N; j++) {
                    int u = N * pgid + j;
                    auto lb = lower_bound(stADD.begin(), stADD.end(), u);
                    if (s == 278 && u == 586581)
                        s = 278;
                    // 这个点在插入边中
                    if (lb < stADD.end() && (*lb) == u) {

                        assert(pgid < lms.back().pages.size());
                        if (lms.back().pages[pgid].elems.empty()) {
                            lms.back().pages[pgid].init(pgid);
                        }
                        std::vector<int> adj;
                        for (; ADDidx < ADD.size(); ADDidx++) {if (ADD[ADDidx].first != u) break; adj.emplace_back(ADD[ADDidx].second);}
                        // if (readd.find(u) != readd.end()) {
                        //     for (auto l : readd[u]) adj.emplace_back(l);
                        //     sort(adj.begin(), adj.end());
                        //     adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                        //     readd.erase(u);
                        // }
                        if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                        if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                        assert(u <  lms.back().ls.size());
                        auto ls_u = lms.back().ls[u];
                        assert(ls_u != -1);
                        if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                        // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                        if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                        lms.back().degree[ls_u] = adj.size();
                        // lms.back().trueD[ls_u] = adj.size();

                        // 这个点存在旧信息
                        EDGE tmp;
                        tmp.v = -1;
                        if (lms.back().pages[pgid].idx[j].first != -1) {
                            if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                            if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                            auto ls_u = lms.back().ls[u];
                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);
                            // 在删除时已经处理了这个点的con链
                            if (!readd.empty() && readd.find(u) != readd.end()) {
                                for (auto l : readd[u]) {
                                    adj.emplace_back(l);
                                }
                                readd.erase(u);
                            }
                            else {// 在删除时未处理这个点的con链
                                // 这个点之前存在con，处理con链
                                auto [u2snapidx, u2Idx] = lms[snapid].pages[pgid].idx[j];
                                auto u2snapls_u = lms[u2snapidx].ls[u];

                                if (egs[u2snapidx].edge_table[lms[u2snapidx].pre[lms[u2snapidx].ls[u] + 1] - 1].is_con) {
                                    /*int conSize = 0;
                                    for (int l = lms[u2snapidx].pre[lms[u2snapidx].ls[u] + 1] - 1; l >= lms[u2snapidx].pre[lms[u2snapidx].ls[u]]; l--) {
                                        if (egs[u2snapidx].edge_table[l].is_con) conSize++;
                                        else break;
                                    }

                                    for (int i = lms[u2snapidx].pre[u2snapls_u], idx = 0; i < lms[u2snapidx].pre[u2snapls_u + 1] - conSize; i++, idx++) {
                                        auto lb = lower_bound(adj.begin(), adj.end(), egs[u2snapidx].edge_table[i].v);
                                        if (lb < adj.end() && *lb == egs[u2snapidx].edge_table[i].v) {

                                        }
                                        else
                                        {
                                            adj.emplace_back(egs[u2snapidx].edge_table[i].v);
                                        }
                                    }

                                    for (int i = lms[u2snapidx].pre[u2snapls_u + 1] - conSize; i < lms[u2snapidx].pre[u2snapls_u + 1]; i++) {
                                        auto eg = egs[u2snapidx].edge_table[i];
                                        auto& st = egs[u2snapidx].edge_table[i].st;
                                        assert(eg.is_con);
                                        for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                                            // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                                            if (egs[eg.egtableIdx].edge_table[l].is_con) {
                                                // std::cout << "ass453" << std::endl;
                                                continue;
                                            }
                                            auto lb = lower_bound(adj.begin(), adj.end(), egs[eg.egtableIdx].edge_table[l].v);
                                            if (st.empty() || !st[idx])
                                                if (lb == adj.end() || *lb != egs[eg.egtableIdx].edge_table[l].v) {
                                                    adj.emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                                            }
                                        }
                                    }
                                    */

                                    getAdj(u, u2snapidx, u2snapls_u, adj);

                                    sort(adj.begin(), adj.end());
                                    adj.erase(unique(adj.begin(), adj.end()), adj.end());
                                    if (!deledge.empty() && deledge.find(u) != deledge.end())
                                        deledge.erase(u);
                                }
                                else {
                                    // 这个点之前不存在con，创建一个con
                                    // assert(deledge.find(u) != deledge.end());
                                    if (deledge.find(u) != deledge.end()) {
                                        tmp = deledge[u];
                                        deledge.erase(u);

                                        if (tmp.len == 1) {
                                            tmp.is_con = false;
                                            adj.emplace_back(egs[tmp.egtableIdx].edge_table[tmp.v].v);
                                            tmp.v = -1;
                                        }
                                    }
                                    else {
                                        if (u2snapls_u == -1 || lms[u2snapidx].degree[u2snapls_u] == 0) continue;
                                        if (tmp.len == 1) {
                                            tmp.is_con = false;
                                            adj.emplace_back(egs[tmp.egtableIdx].edge_table[lms[u2snapidx].pre[u2snapls_u]].v);
                                            tmp.v = -1;
                                        }
                                        else {
                                            tmp.is_con = true;
                                            tmp.v = lms[u2snapidx].pre[u2snapls_u];
                                            tmp.egtableIdx = u2snapidx;
                                            tmp.len = lms[u2snapidx].pre[u2snapls_u + 1] - lms[u2snapidx].pre[u2snapls_u];
                                        }
                                    }
                                }
                            }
                        }

                        sort(adj.begin(), adj.end());
                        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());

                        lms.back().degree[ls_u] = adj.size();
                        for (auto l : adj) egs.back().edge_table.emplace_back(l);
                        if (tmp.v != -1) {
                            assert(tmp.is_con);
                            egs.back().edge_table.emplace_back(std::move(tmp));
                            lms.back().degree[ls_u]++;
                        }
                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                        lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                    }
                    // 这个点不在插入边中
                    else {
                        bool hasadd = false;
                        if (readd.find(u) != readd.end()) {
                            sort(readd[u].begin(), readd[u].end());
                            if (u >= lms.back().ls.size()) {
                                lms.back().ls.resize(u + 1, -1);
                                lms.back().ls[u] = vCount++;
                            }
                            auto ls_u = lms.back().ls[u];
                            if (lms.back().pages[pgid].elems.empty()) {
                                lms.back().pages[pgid].init(pgid);
                            }

                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                            lms.back().degree[ls_u] = readd[u].size();
                            // lms.back().trueD[ls_u] = readd[u].size();
                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                            for (auto l : readd[u]) egs.back().edge_table.emplace_back(l);
                            readd.erase(u);
                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                            hasadd = true;
                        }
                        if (deledge.find(u) != deledge.end()) {
                            if (u >= lms.back().ls.size()) {
                                lms.back().ls.resize(u + 1, -1);
                                lms.back().ls[u] = vCount++;
                            }
                            auto ls_u = lms.back().ls[u];

                            if (lms.back().pages[pgid].elems.empty()) {
                                lms.back().pages[pgid].init(pgid);
                            }

                            if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                            // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                            if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                            egs.back().edge_table.emplace_back(deledge[u]);
                            deledge.erase(u);
                            assert(ls_u >= 0);
                            assert(ls_u < lms.back().degree.size());
                            lms.back().degree[ls_u]++;

                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                            hasadd = true;
                        }
                        if (!hasadd) {
                            // 这个点不在插入边中 但是所在page涉及修改 所以需要维护该点在当前快照当前page中的信息
                            if (!lms[snapid].pages[pgid].idx.empty()) {
                                auto [u2snapidx, u2Idx] = lms[snapid].pages[pgid].idx[j];
                                if (u2snapidx != -1 && ! lms.back().pages[pgid].idx.empty() && lms.back().pages[pgid].idx[j].first != -1) {
                                    auto u2snapls_u = lms[u2snapidx].ls[u];

                                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                                    auto ls_u = lms.back().ls[u];

                                    if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                                    // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                                    if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                                    lms.back().degree[ls_u] = lms[u2snapidx].degree[u2snapls_u];

                                    for (int l = u2Idx; l < u2Idx + lms[u2snapidx].degree[u2snapls_u]; l++) {
                                        egs.back().edge_table.emplace_back(egs[u2snapidx].edge_table[l]);
                                    }
                                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                    lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                                }
                            }
                        }
                        /*else {
                            if (lms.back().pages[pgid].idx[j].first != -1) {
                                if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                                if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                                auto ls_u = lms.back().ls[u];
                                if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                                // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                                if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                                if (egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u] + 1] - 1].is_con) {
                                    int conSize = 0;
                                    while (egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u] + 1] - 1 - conSize].is_con) conSize++;
                                    // std::vector<int> adj;

                                    for (int i = lms[snapid].pre[lms[snapid].ls[u] + 1] - conSize; i < lms[snapid].pre[lms[snapid].ls[u] + 1]; i++) {
                                        auto eg = egs[snapid].edge_table[i];
                                        auto& st = egs[snapid].edge_table[i].st;
                                        // assert(eg.is_con);
                                        for (int l = eg.v, idx = 0; idx < eg.len; l++, idx++) {
                                            // assert(!egs[eg.egtableIdx].edge_table[l].is_con);
                                            if (st.empty() || !st[idx]) {
                                                lms.back().degree[ls_u]++;
                                                // lms.back().trueD[ls_u]++;
                                                egs[snapid].edge_table.emplace_back(egs[eg.egtableIdx].edge_table[l].v);
                                            }
                                        }
                                    }

                                    if (!egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].is_con) {
                                        EDGE tmp;
                                        tmp.is_con = true;
                                        tmp.v = lms[snapid].pre[lms[snapid].ls[u]];
                                        tmp.egtableIdx = snapid;
                                        tmp.len = lms[snapid].pre[lms[snapid].ls[u] + 1] - lms[snapid].pre[lms[snapid].ls[u]];
                                        lms.back().degree[ls_u]++;
                                        // lms.back().trueD[ls_u]++;
                                        egs.back().edge_table.emplace_back(std::move(tmp));
                                    }
                                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                }
                                else {

                                    // 把上一个快照中的相同点的边作为con存入这个快照
                                    EDGE tmp;

                                    tmp.len = lms[snapid].pre[lms[snapid].ls[u] + 1] - lms[snapid].pre[lms[snapid].ls[u]];
                                    assert(tmp.len > 0);
                                    if (tmp.len != 1) {
                                        tmp.is_con = true;
                                        tmp.v = lms[snapid].pre[lms[snapid].ls[u]];
                                        tmp.egtableIdx = snapid;
                                        lms.back().degree[ls_u]++;
                                        // lms.back().trueD[ls_u]++;
                                        lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                                        egs.back().edge_table.emplace_back(std::move(tmp));

                                        tmp.st.resize(tmp.len);
                                        bool dont = true;

                                        for (int l = lms[snapid].pre[lms[snapid].ls[u]], idx = 0; l < lms[snapid].pre[lms[snapid].ls[u] + 1]; l++, idx++) {
                                            if (egs[snapid].edge_table[l].del) {
                                                dont = false;
                                                tmp.st[idx] = true;
                                            }
                                        }
                                        if (!dont) {
                                            tmp.st.clear();
                                            tmp.st.shrink_to_fit();
                                        }
                                    }
                                    else {
                                        tmp.is_con = false;
                                        assert(u < lms[snapid].ls.size());

                                        tmp.v = egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].v;
                                        if (!egs[snapid].edge_table[lms[snapid].pre[lms[snapid].ls[u]]].del) {
                                            tmp.egtableIdx = egs.size() - 1;
                                            egs.back().edge_table.emplace_back(std::move(tmp));
                                            lms.back().degree[ls_u]++;
                                            // lms.back().trueD[ls_u]++;
                                            lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                                            if (lms.back().pages[pgid].elems.empty()) lms.back().pages[pgid].init(pgid);

                                            lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                                        }
                                    }

                                }
                            }
                        }*/
                    }
                }
                lms.back().indirectionTables[indirectionTablesIdx] = {lms.size() - 1, pgid};
            }
            if (readd.size()) {
                assert(1);
                auto u_ = readd.begin();
                for (; u_ != readd.end(); ) {
                    auto u = u_->first;
                    auto j = std::lower_bound(lms.back().pages[pgid].elems.begin(), lms.back().pages[pgid].elems.end(), u) - lms.back().pages[pgid].elems.begin();
                    auto lb = lower_bound(stADD.begin(), stADD.end(), u);
                    // 这个点在插入边中

                    if (lms.back().pages[pgid].elems.empty()) {
                        lms.back().pages[pgid].init(pgid);
                    }
                    std::vector<int> adj;
                    if (lb < stADD.end() && (*lb) == u) {
                        for (; ADDidx < ADD.size(); ADDidx++) {if (ADD[ADDidx].first != u) break; adj.emplace_back(ADD[ADDidx].second);}
                    }
                    if (readd.find(u) != readd.end()) {
                        for (auto l : readd[u]) adj.emplace_back(l);
                        sort(adj.begin(), adj.end());
                        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                        u_ = readd.erase(u_);
                    }
                    else {
                        u_++;
                    }
                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];
                    if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                    // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                    if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                    lms.back().degree[ls_u] = adj.size();
                    // lms.back().trueD[ls_u] = adj.size();

                    for (auto l : adj) egs.back().edge_table.emplace_back(l);


                    if (deledge.find(u) != deledge.end()) {
                        egs.back().edge_table.emplace_back(deledge[u]);
                        deledge.erase(u);
                        lms.back().degree[ls_u]++;
                    }

                    lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                    // 这个点不在插入边中
                    // else {
                    //     if (u >= lms.back().ls.size()) continue;
                    //     assert(readd.find(u) != readd.end());
                    //     auto ls_u = lms.back().ls[u];
                    //     if (deledge.find(u) != deledge.end()) {
                    //         if (lms.back().pages[pgid].elems.empty()) {
                    //             lms.back().pages[pgid].init(pgid);
                    //         }
                    //
                    //         if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                    //         if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                    //         if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);
                    //
                    //         egs.back().edge_table.emplace_back(deledge[u]);
                    //         deledge.erase(u);
                    //         lms.back().degree[ls_u]++;
                    //
                    //         lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                    //
                    //         lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                    //     }
                    //
                    // }
                }
            }
            while (itADD < ADD.size() && stADD[itADD] < N * (pgid + 1)) itADD++;
        }

        if (!readd.empty()) {
            auto u_ = readd.begin();
            for (; u_ != readd.end();) {
                auto u = u_->first;
                auto pgid = u / N;
                auto j = std::lower_bound(lms.back().pages[pgid].elems.begin(), lms.back().pages[pgid].elems.end(), u) - lms.back().pages[pgid].elems.begin();
                auto lb = lower_bound(stADD.begin(), stADD.end(), u);
                // 这个点在插入边中

                    if (lms.back().pages[pgid].elems.empty()) {
                        lms.back().pages[pgid].init(pgid);
                    }
                    std::vector<int> adj;
                    if (lb < stADD.end() && (*lb) == u) {
                        for (; ADDidx < ADD.size(); ADDidx++) {if (ADD[ADDidx].first != u) break; adj.emplace_back(ADD[ADDidx].second);}
                    }
                    if (readd.find(u) != readd.end()) {
                        for (auto l : readd[u]) adj.emplace_back(l);
                        sort(adj.begin(), adj.end());
                        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                        u_ = readd.erase(u_);
                    }
                    else {
                        u_++;
                    }
                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];
                    if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                    // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                    if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                    lms.back().degree[ls_u] = adj.size();
                    // lms.back().trueD[ls_u] = adj.size();

                    for (auto l : adj) egs.back().edge_table.emplace_back(l);


                    if (deledge.find(u) != deledge.end()) {
                        egs.back().edge_table.emplace_back(deledge[u]);
                        deledge.erase(u);
                        lms.back().degree[ls_u]++;
                    }

                    lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                }
                // 这个点不在插入边中
                // else {
                //     // if (!readd.empty())
                //     //     assert(readd.find(u) != readd.end());
                //     // if (u >= lms.back().ls.size()) continue;
                //
                //     auto ls_u = lms.back().ls[u];
                //     if (deledge.find(u) != deledge.end()) {
                //         if (lms.back().pages[pgid].elems.empty()) {
                //             lms.back().pages[pgid].init(pgid);
                //         }
                //
                //         if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                //         if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                //         if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);
                //
                //         egs.back().edge_table.emplace_back(deledge[u]);
                //         deledge.erase(u);
                //         lms.back().degree[ls_u]++;
                //
                //         lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                //
                //         lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                //     }
                //
                // }
            // }
        }
        if (!deledge.empty()) {
            for (auto u_ : deledge) {
                auto u = u_.first;
                if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                auto ls_u = lms.back().ls[u];
                if (ls_u >= lms.back().degree.size()) lms.back().degree.resize(ls_u + 1);
                // if (ls_u >= lms.back().trueD.size()) lms.back().trueD.resize(ls_u + 1);
                if (ls_u + 1 >= lms.back().pre.size()) lms.back().pre.resize(ls_u + 2);

                egs.back().edge_table.emplace_back(deledge[u]);
                lms.back().degree[ls_u]++;
                auto pgid = u / N;
                if (lms.back().pages[pgid].elems.empty())
                    lms.back().pages[pgid].init(pgid);
                auto j = std::lower_bound(lms.back().pages[pgid].elems.begin(), lms.back().pages[pgid].elems.end(), u) - lms.back().pages[pgid].elems.begin();
                lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};

                lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
            }
        }
    }

    /*void ADDADDED(std::vector<std::pair<int, int>>& E, int maxV) {
        int la = lms.size() - 2;
        lms.back().indirectionTables = lms[la].indirectionTables;

        sort(E.begin(), E.end());
        std::vector<std::vector<int>> adj;
        lms.back().ls.resize(maxV + 1, -1);
        int vCount = 0;

        std::vector<int> st;
        st.reserve(E.size());
        for (auto [u, v] : E) st.emplace_back(u);
        std::sort(st.begin(), st.end());
        st.erase(std::unique(st.begin(), st.end()), st.end());
        if (maxV >= lms.back().indirectionTables.size()) {
            lms.back().indirectionTables.resize((maxV + lms.back().page_size_) / lms.back().page_size_, {-1, -1});
        }

        egs.back().delvec.resize(maxV + 1);

        auto it = 0;
        int Eidx = 0;
        int indirectionTablesIdx = 0;
        for (; indirectionTablesIdx < lms.back().indirectionTables.size() && it != st.size(); indirectionTablesIdx++) {
            auto& [snapid, pgid] = lms.back().indirectionTables[indirectionTablesIdx];
            if (snapid == -1) break;
            // 跳过未更新的页
            if (lms[snapid].pages[pgid].elems.back() < st[it]) continue;

            // 第snapid个快照的第pgid页需要更新
            // 当前快照拓展页
            if (pgid >= lms.back().pages.size()) lms.back().pages.resize(pgid + 1);

            // lms.back().pages.emplace_back(Page{});
            // 拷贝第snapid个快照的第pgid页到当前快照
            // lms.back().pages.back() = lms[snapid].pages[pgid];
            lms.back().pages[pgid].elems = lms[snapid].pages[pgid].elems;
            lms.back().pages[pgid].idx = lms[snapid].pages[pgid].idx;
            // lms.back().pages[pgid].loc = lms[snapid].pages[pgid].loc;
            lms.back().pages[pgid].vsize = lms[snapid].pages[pgid].vsize;

            if (lms.back().pages[pgid].elems.back() <= maxV) {
                lms.back().pages[pgid].vsize = N;
            }

            // 更新当前当前快照的页
            for (int j = 0; j < N; j++) {
                auto u = lms.back().pages[pgid].elems[j];

                // 当前点不用更新
                auto lb = lower_bound(st.begin(), st.end(), u);
                if (lb < st.end() && (*lb) == u) {
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];

                    adj.resize(vCount);
                    lms.back().degree.resize(vCount);
                    lms.back().trueD.resize(vCount);
                    lms.back().pre.resize(vCount + 1);
                    lms.back().fraSize.resize(vCount);
                    for (; Eidx < E.size(); Eidx++) {
                        if (E[Eidx].first != u) break;
                        adj[ls_u].emplace_back(E[Eidx].second);
                    }

                    lms.back().degree[ls_u] = adj[ls_u].size();
                    lms.back().trueD[ls_u] = adj[ls_u].size();
                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                    egs.back().edge_table.reserve(egs.back().edge_table.size() + adj[ls_u].size());
                    for (auto tmp : adj[ls_u]) egs.back().edge_table.emplace_back(-1, tmp, 0);
                    // std::move(adj[ls_u].begin(), adj[ls_u].end(), std::back_inserter(egs.back().edge_table));
                    // if (lms.back().pages[pgid].idx[j].first != -1) {
                    if (lms.back().pages[pgid].idx[j].first < lms[la].ls.size() &&
                        lms[la].ls[lms.back().pages[pgid].idx[j].first] != -1 &&
                        lms[la].degree[lms.back().pages[pgid].idx[j].first] > 0) {
                        lms.back().degree[ls_u]++;
                        // lms.back().trueD[ls_u]++;
                        lms.back().pre[ls_u + 1]++;
                        lms.back().fraSize[ls_u]++;

                        auto [a, b] = lms.back().pages[pgid].idx[j];
                        auto len = lms[a].pre[lms[a].ls[u] + 1] - lms[a].pre[lms[a].ls[u]];
                        lms.back().trueD[ls_u] += lms[a].trueD[lms[a].ls[u]];
                        egs.back().edge_table.emplace_back(a, b, len);
                    }
                    else {
                        // assert(lms.back().pages[pgid].idx[j].first == -1 && lms.back().pages[pgid].idx[j].second == -1);
                    }

                    lms.back().pages[pgid].idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                }
                // else if (lms.back().pages[pgid].idx[j].first != -1) {
                else if (lms.back().pages[pgid].idx[j].first < lms[la].ls.size() &&
                        lms[la].ls[lms.back().pages[pgid].idx[j].first] != -1 &&
                        lms[la].degree[lms.back().pages[pgid].idx[j].first] > 0) {
                    continuation t;
                    t.a = la;
                    t.b = lms[la].pre[lms[la].ls[u]];
                    t.len = lms[la].pre[lms[la].ls[u] + 1] - lms[la].pre[lms[la].ls[u]];
                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];
                    if (egs[la].edge_table[lms[la].pre[lms[la].ls[u]]].a != -1) {
                        assert(t.len == 1);
                        t = egs[la].edge_table[lms[la].pre[lms[la].ls[u]]];
                    }
                    lms.back().degree.resize(ls_u + 1);
                    lms.back().trueD.resize(ls_u + 1);
                    lms.back().pre.resize(ls_u + 2);
                    lms.back().fraSize.resize(ls_u + 1);

                    lms.back().degree[ls_u]++;
                    lms.back().trueD[ls_u] += lms[la].trueD[lms[la].ls[u]];
                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                    lms.back().fraSize[ls_u] = 1;
                    egs.back().edge_table.emplace_back(t);
                }
            }
            snapid = lms.size() - 1;

            assert(snapid >= 0 && pgid >= 0);
            // 跳到下一页
            while (it < st.size() && st[it] <= lms[snapid].pages[pgid].elems.back()) it++;
        }


        for (; indirectionTablesIdx < lms.back().indirectionTables.size(); indirectionTablesIdx++) {
            lms.back().pages.emplace_back();
            assert(N + 1 < lms.back().indirectionTables.size());
            for (int i = N * indirectionTablesIdx, j = 0; i <  N * (indirectionTablesIdx + 1); i++) {
                assert(j < lms.back().pages.back().elems.size());

                lms.back().pages.back().elems[j++] = i;
            }

            auto snapid = lms.size() - 1;
            auto pgid = lms.back().pages.size() - 1;
            lms.back().indirectionTables[indirectionTablesIdx] = {snapid, pgid};
            if (lms[snapid].pages[pgid].elems.back() < st[it]) continue;

            for (int j = 0; j < lms.back().pages.back().elems.size(); j++) {
                auto u = lms.back().pages.back().elems[j];
                // if (u > st.back()) {
                //     break;
                // }
                lms.back().pages.back().vsize = j + 1;
                // 当前点不用更新
                auto lb = lower_bound(st.begin(), st.end(), u);
                if (lb < st.end() && (*lb) == u) {
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];

                    adj.resize(vCount);
                    lms.back().degree.resize(vCount);
                    lms.back().trueD.resize(vCount);
                    lms.back().pre.resize(vCount + 1);
                    lms.back().fraSize.resize(vCount);

                    for (; Eidx < E.size(); Eidx++) {
                        if (E[Eidx].first != u) break;
                        adj[ls_u].emplace_back(E[Eidx].second);
                    }

                    lms.back().degree[ls_u] = adj[ls_u].size();
                    lms.back().trueD[ls_u] = adj[ls_u].size();
                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];

                    // egs.back().edge_table.reserve(egs.back().edge_table.size() + adj[ls_u].size());
                    for (auto tmp : adj[ls_u]) egs.back().edge_table.emplace_back(-1, tmp, 0);
                    // std::move(adj[ls_u].begin(), adj[ls_u].end(), std::back_inserter(egs.back().edge_table));

                    if (lms.back().pages.back().idx[j].first != -1) {
                        lms.back().degree[ls_u]++;
                        lms.back().pre[ls_u + 1]++;
                        lms.back().fraSize[ls_u]++;
                        assert(j < lms.back().pages.back().idx.size());
                        auto [a, b] = lms.back().pages.back().idx[j];
                        assert(a >= 0 && b >= 0 && a < lms.size() && u < lms[a].ls.size());
                        auto len = lms[a].degree[lms[a].ls[u]];
                        lms.back().trueD[ls_u] += lms[a].trueD[lms[a].ls[u]];
                        egs.back().edge_table.emplace_back(a, b, len);
                    }
                    else {
                        assert(lms.back().pages.back().idx[j].first == -1 && lms.back().pages.back().idx[j].second == -1);
                    }

                    lms.back().pages.back().idx[j] = {lms.size() - 1, lms.back().pre[ls_u]};
                }
                else if (u < lms[la].ls.size()) {
                    assert(lms[la].ls[u] != -1);
                    continuation t;
                    t.a = la;
                    t.b = lms[la].pre[lms[la].ls[u]];
                    t.len = lms[la].pre[lms[la].ls[u] + 1] - lms[la].pre[lms[la].ls[u]];
                    if (u >= lms.back().ls.size()) lms.back().ls.resize(u + 1, -1);
                    if (lms.back().ls[u] == -1) lms.back().ls[u] = vCount++;
                    auto ls_u = lms.back().ls[u];
                    if (egs[la].edge_table[lms[la].pre[lms[la].ls[u]]].a != -1) {
                        assert(t.len == 1);
                        t = egs[la].edge_table[lms[la].pre[lms[la].ls[u]]];
                    }
                    lms.back().degree.resize(ls_u + 1);
                    lms.back().trueD.resize(ls_u + 1);
                    lms.back().pre.resize(ls_u + 2);
                    lms.back().fraSize.resize(ls_u + 1);

                    lms.back().degree[ls_u]++;
                    lms.back().trueD[ls_u] += lms[la].trueD[lms[la].ls[u]];
                    lms.back().pre[ls_u + 1] = lms.back().pre[ls_u] + lms.back().degree[ls_u];
                    lms.back().fraSize[ls_u] = 1;
                    egs.back().edge_table.emplace_back(t);
                }
            }
        }

        lms.back().vertexCount = vCount;
    }*/

    std::vector<double> times;

    std::vector<std::pair<int, int>> tmpv;
    void searchBack(int snapIdx, int vIdx, int len, int u, int lasnapidx, int lai) {
        for (int i = vIdx, j = 0; i < vIdx + len; i++, j++) {
            auto a = egs[snapIdx].edge_table[i].is_con;
            auto b = egs[snapIdx].edge_table[i].v;

            // if (!a && !egs[snapIdx].edge_table[i].del) {
            if (!a) {
                if (lasnapidx == -1 && !egs[snapIdx].edge_table[i].del) {
                    tmpv.emplace_back(u, b);
                    // std::cout << u << " " << b << std::endl;
                    EdgeCount++;
                }
                else {
                    if (egs[lasnapidx].edge_table[lai].st.empty() || !egs[lasnapidx].edge_table[lai].st.empty() && !egs[lasnapidx].edge_table[lai].st[j]) {
                        tmpv.emplace_back(u, b);
                        // std::cout << u << " " << b << std::endl;
                        EdgeCount++;
                    }
                }
            }
            else if (a){
                searchBack(egs[snapIdx].edge_table[i].egtableIdx, b, egs[snapIdx].edge_table[i].len, u, snapIdx, i);
            }
        }
    }

    int PTeg(int idx) {
        // std::cout << idx << " snapshot's Vertex: " << std::endl;
        // int V = 0;
        // for (auto i : lms[idx].indirectionTables) {
        //     for (int j = 0; j < N; j++) {
        //         V++;
        //         if (lms[i.first].pages[i.second].elems.empty()) continue;
        //         auto u = lms[i.first].pages[i.second].elems[j];
        //         auto ls_u = lms[i.first].ls[u];
        //         if (lms[i.first].ls[u] != -1) {
        //             std::cout << u << " > " << lms[i.first].trueD[ls_u] << std::endl;
        //         }
        //     }
        // }
        // std::cout << V << std::endl;

        // std::cout << idx << " snapshot's Edge: ";
        tmpv.clear();
        tmpv.reserve(egs[idx].edge_table.size());
        EdgeCount = 0;
        for (auto i : lms[idx].indirectionTables) {
            if (i.first == -1) continue;
            for (int j = 0; j < lms[i.first].pages[i.second].elems.size(); j++) {
                auto u = lms[i.first].pages[i.second].elems[j];
                if (lms[i.first].pages[i.second].idx[j].first == -1) continue;
                if (u >= lms[lms[i.first].pages[i.second].idx[j].first].ls.size()) {
                    // std::cout << " break " << std::endl;
                    continue;
                }

                auto ls_u = lms[lms[i.first].pages[i.second].idx[j].first].ls[u];
                if (ls_u == -1) continue;
                searchBack(lms[i.first].pages[i.second].idx[j].first, lms[lms[i.first].pages[i.second].idx[j].first].pre[ls_u],
                    lms[lms[i.first].pages[i.second].idx[j].first].pre[ls_u + 1] - lms[lms[i.first].pages[i.second].idx[j].first].pre[ls_u], u, -1, -1);
            }
        }
        // std::cout << std::endl;
        // std::cout << EdgeCount << std::endl;
        return EdgeCount;
    }

    double getMemoryMB() {
        size_t sz = 0;
        for (auto i : lms) {
            sz += i.degree.size();
            sz += i.pre.size();
            sz += i.ls.size();
            sz += i.indirectionTables.size() * 2;
            for (auto j : i.pages) {
                sz += j.elems.size();
                sz += j.idx.size() * 2;
            }
        }
        for (const auto& i : egs) {
            for (const auto& j : i.edge_table) {
                sz += j.st.size();
                sz += 6;
            }
        }
        return 1.0 * sz * 4 / 1024 / 1024;
    }

    std::vector<std::pair<int, int>> res;
    int edges = 0;
    void query(int l, int r, bool fun) {
        edges = 0;
        for (int i = l; i <= r; i++) {
            if (i == 1159)
                std::cout << i << std::endl;
            auto eC = PTeg(i);
            std::sort(tmpv.begin(), tmpv.end());
            if (fun) {
                std::vector<std::pair<int, int>> tmp;
                std::set_union(res.begin(), res.end(),
                               tmpv.begin(), tmpv.end(),
                               std::back_inserter(tmp));
                res = std::move(tmp);
            }
            else {
                if (i == l) {
                    res = std::move(tmpv);
                } else {
                    std::vector<std::pair<int, int>> tmp;
                    std::set_intersection(res.begin(), res.end(),
                                                tmpv.begin(), tmpv.end(),
                                          std::back_inserter(tmp));
                    res = move(tmp);
                }
            }
        }
        for (const auto &i: res) {
            edges++;
        }
        res.clear();
        std::cout << l << " " << r << " " << edges << std::endl;
    }
};
