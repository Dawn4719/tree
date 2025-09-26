#ifndef TYPE_H
#define TYPE_H

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdlib.h>
#include <vector>
typedef unsigned int ui;
typedef ui LabelID;
typedef uint32_t VertexID;
typedef uint64_t ul;
#include <cmath>
using namespace std;

class BloomFilter {
  private:
	size_t m; // bit 数组大小
	size_t k; // 哈希函数个数
	std::vector<uint64_t> bits;

	inline bool getBit(size_t idx) const {
		return (bits[idx >> 6] >> (idx & 63)) & 1ULL;
	}

	inline void setBit(size_t idx) {
		bits[idx >> 6] |= (1ULL << (idx & 63));
	}

	static uint64_t hash1(uint64_t x) {
		x ^= (x >> 33);
		x *= 0xff51afd7ed558ccdULL;
		x ^= (x >> 33);
		x *= 0xc4ceb9fe1a85ec53ULL;
		x ^= (x >> 33);
		return x;
	}

	static uint64_t hash2(uint64_t x) {
		x = (~x) + (x << 21);
		x ^= (x >> 24);
		x = (x + (x << 3)) + (x << 8);
		x ^= (x >> 14);
		x = (x + (x << 2)) + (x << 4);
		x ^= (x >> 28);
		x += (x << 31);
		return x;
	}

	std::pair<uint64_t, uint64_t> aggregate(const ui *arr, ui n) const {
		uint64_t h1 = 0, h2 = 0;
		for (ui i = 0; i < n; i++) {
			uint64_t v = arr[i];
			h1 ^= hash1(v);
			h2 += hash2(v);
		}
		return {h1, h2};
	}

  public:
	BloomFilter(size_t m_bits, size_t k_hash) : m(m_bits), k(k_hash) {
		bits.resize((m + 63) / 64);
	}

	BloomFilter(size_t expected_elements, double false_positive_rate = 0.01) {
		if (expected_elements == 0)
			expected_elements = 1000000; // 默认100万
		if (false_positive_rate <= 0 || false_positive_rate >= 1)
			false_positive_rate = 0.01;

		m = static_cast<size_t>(ceil(-(double)expected_elements * log(false_positive_rate) / (log(2) * log(2))));
		k = static_cast<size_t>(round((double)m / expected_elements * log(2)));
		if (k < 1)
			k = 1;
		if (k > 30)
			k = 30;

		bits.resize((m + 63) / 64);
	}

	BloomFilter() : BloomFilter(1000000, 0.01) {}

	bool possiblyContainsArray(const ui *arr, ui n) const {
		auto [h1, h2] = aggregate(arr, n);
		for (size_t i = 0; i < k; i++) {
			size_t pos = (h1 + i * h2) % m;
			if (!getBit(pos))
				return false;
		}
		return true;
	}

	void insertArray(const ui *arr, ui n) {
		auto [h1, h2] = aggregate(arr, n);
		for (size_t i = 0; i < k; i++) {
			size_t pos = (h1 + i * h2) % m;
			setBit(pos);
		}
	}

	bool checkAndInsertArray(const ui *arr, ui n) {
		auto [h1, h2] = aggregate(arr, n);
		bool exists = true;
		for (size_t i = 0; i < k; i++) {
			size_t pos = (h1 + i * h2) % m;
			if (!getBit(pos))
				exists = false;
		}
		if (!exists) {
			for (size_t i = 0; i < k; i++) {
				size_t pos = (h1 + i * h2) % m;
				setBit(pos);
			}
		}
		return exists;
	}

	void clear() {
		std::fill(bits.begin(), bits.end(), 0ULL);
	}
};

struct BlackHolePoint {
    std::vector<VertexID> PartBlackHole;
    std::vector<VertexID> endPoint;
    VertexID id;

    bool operator<(const BlackHolePoint &other) const {
        if (endPoint.size() != other.endPoint.size())
            return endPoint.size() > other.endPoint.size();  // 大的在前

        if (endPoint != other.endPoint)
            return endPoint > other.endPoint;  // 反过来，字典序大的在前

        if (PartBlackHole.size() != other.PartBlackHole.size())
            return PartBlackHole.size() > other.PartBlackHole.size();  // 大的在前

        return PartBlackHole > other.PartBlackHole;  // 字典序大的在前
    }
};


#endif