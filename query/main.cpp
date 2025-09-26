#include <chrono>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include "../utils/CLI11.hpp"
#include "../utils/types.h"
#include <boost/dynamic_bitset.hpp>
#include "deltagraph.h"
#include "ViLa.h"
#include <fstream>
#include <random>
#include <ctime>
#include <rocksdb/db.h>
#include <rocksdb/options.h>
#include <rocksdb/slice.h>
#include "staticcore.hpp"
#include "LLAMA.hpp"
#include <utility/utils.h>
#include "CECI.hpp"
#include "baseline.hpp"
#include "IIT.hpp"
#include "IITR.hpp"
#include "Delta.h"
using namespace rocksdb;
using namespace std;
int K = 18;
int level;
bool fun;
std::string query_path, dataset, input_file;
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

void DG_method(bool fun) {
    auto start = Get_Time();
    DeltaGraph DG(4, K);
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
            for (auto l : j.second) {
                intSz += l.csr.capacity();
            }
        }
    }
    auto query_start_time = Get_Time();
    DG.query(input_file, K, fun, level);
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_start_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << std::endl;

}

void ViLa_method(bool fun) {
    auto start = Get_Time();
    ViLa_Graph ViLa;
    ViLa.build(K, query_path);
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    auto query_time = Get_Time();
    ViLa.query(input_file, fun);

    size_t intSz = 0;
    for (auto& i : ViLa.ViLaEdge) {
        for (auto& j: i) {
            intSz ++;
            intSz += j.lifespan.num_blocks();
        }
    }
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << "Index Memory:" << intSz * 4 / 1024 / 1024 << "mb" << std::endl;
}

void IIT_method(int alpha, int beta, bool fun, string moreQuery) {
    IIT iit(K, alpha, beta, fun, query_path, dataset, input_file);
    auto start = Get_Time();
    iit.load();
    auto query_start_time = Get_Time();
    iit.query(moreQuery);
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_start_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << '\t' << "Index Memory:" << iit.getIndexMemory() << "mb" << std::endl;
}

void IITR_method(int alpha, int beta, bool fun) {
    IITR iitr(K, alpha, beta, fun, query_path, dataset, input_file);
    auto start = Get_Time();
    iitr.load();
    auto query_start_time = Get_Time();
    iitr.query();
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_start_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << '\t' << "Index Memory:" << iitr.getIndexMemory() << "mb" << std::endl;
}

void LLAMA_method(bool fun) {
    auto start = Get_Time();
    LLAMA llama(K);
    llama.load(query_path);
    auto query_time = Get_Time();
    llama.query(query_path, fun);
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000
                << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << '\t' << "Index Memory:" << llama.getMemoryMB() << "mb" << std::endl;
}

void Delta_method(bool fun) {
    GraphDelta dt;
    auto buildTimebegin = Get_Time();

    for (ui i = 0; i < K; i++) {
        std::string path = query_path + std::to_string(i);
        dt.lordDateset(path);
    }
    auto buildTime = Duration(buildTimebegin);

    auto query_time = Get_Time();
    auto IndexMemory = dt.getMemoryMB();
    int count = 0;
    std::fstream input(input_file, std::ios::in);
    int l, r;
    while (input >> l >> r) {
        // cout << l << " " << r << " start_query\n";
        vector<pair<ui, ui>> ans;
        if (fun == true) {
            ans = dt.queryAnd(l + 1, r + 1);
            // cout << "Query：" << l << " " << r << " ........And RES:" << ans.size() << "\n";
        } else {
            ans = dt.queryOr(l + 1, r + 1);
            // cout << "Query：" << l << " " << r << " ........OR RES:" << ans.size() << "\n";
        }
        count = 0;
        for (auto &it : ans) {
            count++;
        }
    }
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000
                << "s" << '\t' << "ALL Time: " << Duration(buildTimebegin) / 1000 << '\t' << "Index Memory:" << IndexMemory << "mb" << std::endl;
}

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

int main(int argc, char *argv[]) {
    auto memst = mem::getValue();
    cout << memst << endl;

    time_t timep;
    time(&timep);
    printf("%s", ctime(&timep));
    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << std::endl;
    CLI::App app{"App description"};
    string method = "IIT";
    dataset = "mo1M";
    level = 0;
    fun = 1;
    int alpha = 2;
    int beta = 64;
    int UPDATECOUNT;
    string moreQuery = "KCORE";
    app.add_option("-d,--dataset", dataset, "query graph path")->required();
    app.add_option("-m,--method", method, "method")->required();
    app.add_option("-f,--fun", fun, "method")->required();
    app.add_option("-a,--alpha", alpha, "alpha")->required();
    app.add_option("-b,--beta", beta, "beta")->required();
    app.add_option("-u,--update", UPDATECOUNT, "update");
    app.add_option("-q,--moreQuery", moreQuery, "moreQuery");
    CLI11_PARSE(app, argc, argv);

    getK();

    query_path = "/home/qsl/exp/IIT/dataset/" + dataset + "/q";
    input_file = "/home/qsl/exp/IIT/dataset/" + dataset + "/input.txt";

    if (moreQuery == "KCORE")
        input_file = "/home/qsl/exp/IIT/dataset/" + dataset + "/inputcore.txt";
    if (moreQuery == "CECI")
        input_file = "/home/qsl/exp/IIT/dataset/" + dataset + "/inputCECI.txt";

    cout << query_path << " " << input_file << " " << method << " K=" << K << endl;
    std::cout << "----------- Loading graphs ------------" << std::endl;
    if (method == "Base") {
        std::cout << "base" << std::endl;
        Base_method(K, query_path, dataset, input_file, fun);
        return 0;
    }
    if (method == "IIT") {
        std::cout << "IIT" << std::endl;
        IIT_method(alpha, beta, fun, moreQuery);
        return 0;
    }
    if (method == "IIT-R") {
        std::cout << "IIT-R" << std::endl;
        Base_method(K, query_path, dataset, input_file, fun);
        return 0;
    }
    if (method == "DG") {
        std::cout << "DG" << std::endl;
        DG_method(fun);
        return 0;
    }
    if (method == "ViLa") {
        std::cout << "ViLa" << std::endl;
        ViLa_method(fun);
        return 0;
    }
    if (method == "LLAMA") {
        std::cout << "LLAMA" << std::endl;
        LLAMA_method(fun);
        return 0;
    }
    if (method == "DELTA") {
        std::cout << "DELTA" << std::endl;
        Delta_method(fun);
        return 0;
    }

    return 0;
}
