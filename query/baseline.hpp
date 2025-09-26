#include <utility/utils.h>
#include "../utils/CLI11.hpp"
#include "../utils/types.h"

void Base_method(int K, string query_path, string dataset, string input_file, bool fun) {
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
        string ss;
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
    cout << "store in DB" << endl;
    return;
    cout << "Build Time: " << Duration(start) / 1000 << endl;
    std::cout << "Memory: " << mem::getValue() / 1024 << "mb" << endl;

    uint L, R;
    cout << input_file << endl;
    fstream input(input_file, ios::in);

    size_t intSz = 0;
    for (int snap_i = 0; snap_i < K; ++snap_i) {
        intSz += snaps[snap_i].size() * 2;
    }

    auto query_time = Get_Time();
    auto query_start_time = Get_Time();

    while (input >> L >> R) {
        cout << "Query: " << L << " " << R;
        int all_edges = 0;
        std::vector<std::pair<int, int>> res;
        for (int snap_i = L; snap_i <= R; ++snap_i) {
            if (fun) {
                std::vector<pair<int, int>> tmp;
                std::set_union(res.begin(), res.end(),
                               snaps[snap_i].begin(), snaps[snap_i].end(),
                               std::back_inserter(tmp));
                res = move(tmp);
            }
            else {
                if (snap_i == L) {
                    res = snaps[snap_i];
                } else {
                    vector<pair<int, int>> tmp;
                    std::set_intersection(res.begin(), res.end(),
                                          snaps[snap_i].begin(), snaps[snap_i].end(),
                                          std::back_inserter(tmp));
                    res = move(tmp);
                }
            }
        }

        for (const auto &i: res) {
            all_edges++;
        }

        std::cout << "  ----------------------------------------------Res: " << all_edges << std::endl;
        all_edges = 0;
    }
    input.close();

    std::cout << "Peak Mem: " << mem::getValue() / 1024 << "mb" << '\t' << "Query time:" << Duration(query_time) / 1000
            << "s" << '\t' << "ALL Time: " << Duration(start) / 1000 << '\t' << "Index Memory:" << intSz * 4 / 1024 / 1024 << "mb" << std::endl;
}