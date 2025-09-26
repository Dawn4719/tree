#include "Delta.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>
void findDeleted(const std::vector<ui> &a, const std::vector<ui> &b, std::vector<std::pair<ui, ui>> &ans, ui x) {
	if (a.empty())
		return;
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j]) {
			if (x < a[i]) {
				ans.emplace_back(x, a[i]);
			}
			++i;
		} else if (a[i] > b[j]) {
			++j;
		} else {
			++i;
			++j;
		}
	}
	while (i < a.size()) {
		if (x < a[i]) {
			ans.emplace_back(x, a[i]);
		}
		++i;
	}
}
void findAdded(const std::vector<ui> &a, const std::vector<ui> &b, std::vector<std::pair<ui, ui>> &ans, ui x) {
	if (b.empty())
		return;

	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j]) {
			++i;
		} else if (a[i] > b[j]) {
			if (x < b[j]) {
				ans.emplace_back(x, b[j]);
			}
			++j;
		} else {
			++i;
			++j;
		}
	}
	while (j < b.size()) {
		if (x < b[j]) {
			ans.emplace_back(x, b[j]);
		}
		++j;
	}
}

void GraphDelta::lordDateset(const std::string &file_path) {
	std::ifstream infile(file_path);
	if (!infile.is_open()) {
		std::cout << "Can not open the graph file " << file_path << std::endl;
		exit(-1);
	}
	sumT_++;
	ui u, v;
	std::vector<std::pair<ui, ui>> edge;
	ui maxID = 0;
	while (infile >> u >> v) {
		edge.emplace_back(u, v);
		maxID = std::max({maxID, u, v});
	}

	std::vector<std::vector<ui>> now_neighbor;
	vertices_count_ = std::max(maxID + 1, vertices_count_);
	edges_count_ = edge.size();

	now_neighbor.resize(vertices_count_);

	for (auto [u, v] : edge) {
		now_neighbor[u].push_back(v);
		now_neighbor[v].push_back(u);
	}

	for (ui i = 0; i < vertices_count_; i++) {
		std::sort(now_neighbor[i].begin(), now_neighbor[i].end());
	}

	if (sumT_ == 1) {
		neighbor = now_neighbor;
		last_neighbor = now_neighbor;
		addEdge_first.push_back({});
		deleteEdge_first.push_back({});
		addEdge_last.push_back({});
		deleteEdge_last.push_back({});

		addEdge_first.push_back({});
		deleteEdge_first.push_back({});
		addEdge_last.push_back({});
		deleteEdge_last.push_back({});
		first_vertices_count_ = vertices_count_;
		last_vertices_count_ = vertices_count_;
	} else {
		std::vector<std::pair<ui, ui>> first_add, first_delete;
		std::vector<std::pair<ui, ui>> last_add, last_delete;

		// 取所有可能点数的最大值，避免漏点
		ui max_vertices = std::max({(ui)neighbor.size(), (ui)last_neighbor.size(), vertices_count_});

		for (ui i = 0; i < max_vertices; i++) {
			const std::vector<ui> emptyVec; // 用来兜底

			const auto &cur = (i < now_neighbor.size()) ? now_neighbor[i] : emptyVec;
			const auto &pre_first = (i < neighbor.size()) ? neighbor[i] : emptyVec;
			const auto &pre_last = (i < last_neighbor.size()) ? last_neighbor[i] : emptyVec;

			// first 差分
			findDeleted(pre_first, cur, first_delete, i);
			findAdded(pre_first, cur, first_add, i);

			// last 差分
			findDeleted(pre_last, cur, last_delete, i);
			findAdded(pre_last, cur, last_add, i);
		}

		std::sort(first_add.begin(), first_add.end());
		std::sort(first_delete.begin(), first_delete.end());
		std::sort(last_add.begin(), last_add.end());
		std::sort(last_delete.begin(), last_delete.end());

		addEdge_first.push_back(first_add);
		deleteEdge_first.push_back(first_delete);
		addEdge_last.push_back(last_add);
		deleteEdge_last.push_back(last_delete);

		last_neighbor = now_neighbor;
		last_vertices_count_ = vertices_count_;
	}

	infile.close();
	// std::cout << "lord_over :" << file_path << "\n";
}

bool existsPair(const std::vector<std::pair<ui, ui>> &vec, const std::pair<ui, ui> &target) {
	return std::binary_search(vec.begin(), vec.end(), target);
}

std::vector<std::pair<ui, ui>> uniqueEdges(const std::vector<std::pair<ui, ui>> &edges) {
	std::set<std::pair<ui, ui>> s(edges.begin(), edges.end());
	return std::vector<std::pair<ui, ui>>(s.begin(), s.end());
}

void GraphDelta::getGraphL(ui L, std::vector<std::pair<ui, ui>> &ans) {
	if (L >= deleteEdge_first.size() || L >= addEdge_first.size())
		return;

	ui n = neighbor.size();
	for (ui i = 0; i < n; i++) {
		for (auto it : neighbor[i]) {
			if (it < i)
				continue;
			if (existsPair(deleteEdge_first[L], {i, it}))
				continue;
			ans.emplace_back(i, it);
		}
	}
	for (auto it : addEdge_first[L]) {
		ans.emplace_back(it);
	}
}

// queryOr 返回去重后的边
std::vector<std::pair<ui, ui>> GraphDelta::queryOr(ui l, ui r) {
	if (l > sumT_ || r > sumT_) {
		std::cout << "NO this range\n";
		exit(0);
	}

	std::vector<std::pair<ui, ui>> ans;
	getGraphL(l, ans);

	for (ui i = l + 1; i <= r; i++) {
		if (i < addEdge_last.size()) {
			ans.insert(ans.end(), addEdge_last[i].begin(), addEdge_last[i].end());
		}
	}

	return uniqueEdges(ans);
}

std::vector<std::pair<ui, ui>> GraphDelta::queryAnd(ui l, ui r) {
	if (l > sumT_ || r > sumT_) {
		std::cout << "NO this range\n";
		exit(0);
	}

	// 第 l 个图的边集
	std::vector<std::pair<ui, ui>> ans;
	getGraphL(l, ans);

	// alive[i] 表示 ans[i] 是否还存活
	std::vector<bool> alive(ans.size(), true);

	// 从 l+1 到 r，不断删除被删除的边
	for (ui i = l + 1; i <= r; i++) {
		for (size_t j = 0; j < ans.size(); j++) {
			if (alive[j] && existsPair(deleteEdge_last[i], ans[j])) {
				alive[j] = false;
			}
		}
	}

	// 输出存活边
	std::vector<std::pair<ui, ui>> res;
	for (size_t j = 0; j < ans.size(); j++) {
		if (alive[j])
			res.push_back(ans[j]);
	}

	return uniqueEdges(res); // 去重
}

double GraphDelta::getMemoryMB() {

	auto sizeVecVecUI = [](const std::vector<std::vector<ui>> &vv) {
		size_t sum = 0;
		for (const auto &v : vv) {
			sum += v.size() * sizeof(ui);
		}
		return sum;
	};

	auto sizeVecVecPair = [](const std::vector<std::vector<std::pair<ui, ui>>> &vv) {
		size_t sum = 0;
		for (const auto &v : vv) {
			sum += v.size() * sizeof(std::pair<ui, ui>);
		}
		return sum;
	};
	size_t totalBytes = 0;
	totalBytes += sizeVecVecUI(neighbor);
	totalBytes += sizeVecVecUI(last_neighbor);
	totalBytes += sizeVecVecPair(addEdge_first);
	totalBytes += sizeVecVecPair(deleteEdge_first);
	totalBytes += sizeVecVecPair(addEdge_last);
	totalBytes += sizeVecVecPair(deleteEdge_last);

	return (double)(totalBytes) / (1024.0 * 1024.0);
}
