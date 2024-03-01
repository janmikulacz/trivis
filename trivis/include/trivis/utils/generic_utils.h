/**
 * File:   generic_utils.h
 *
 * Date:   25.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_UTILS_GENERIC_UTILS_H_
#define TRIVIS_UTILS_GENERIC_UTILS_H_

#include <vector>
#include <numeric>
#include <cassert>
#include <limits>
#include <algorithm>

namespace trivis::utils {

template<typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
[[nodiscard]] std::vector<Tp1, Alloc1> SortBy(
    const std::vector<Tp1, Alloc1> &vin,
    const std::vector<Tp2, Alloc2> &keys) {
    std::vector<std::size_t> is;
    is.reserve(vin.size());
    for (__attribute__((unused)) auto &&unused: keys) is.push_back(is.size());
    std::sort(begin(is), end(is), [&](std::size_t l, std::size_t r) { return keys[l] < keys[r]; });
    std::vector<Tp1, Alloc1> r;
    r.reserve(vin.size());
    for (std::size_t i: is) r.push_back(vin[i]);
    return r;
}

template<typename Tp, typename Alloc>
Tp MaxValue(const std::vector<std::vector<Tp, Alloc>> &mat) {
    auto max_val = std::numeric_limits<Tp>::lowest();
    for (const auto &vec : mat) max_val = std::max(max_val, *std::max_element(vec.begin(), vec.end()));
    return max_val;
}

template<typename Tp, typename Alloc>
Tp MinValue(const std::vector<std::vector<Tp, Alloc>> &mat) {
    auto min_val = std::numeric_limits<Tp>::max();
    for (const auto &vec : mat) min_val = std::min(min_val, *std::min_element(vec.begin(), vec.end()));
    return min_val;
}

template<typename Tp, typename Alloc>
Tp MaxValue(const std::vector<std::vector<std::vector<Tp, Alloc>>> &cub) {
    auto max_val = std::numeric_limits<Tp>::lowest();
    for (const auto &mat : cub) max_val = std::max(max_val, MaxValue(mat));
    return max_val;
}

template<typename Tp, typename Alloc>
Tp MinValue(const std::vector<std::vector<std::vector<Tp, Alloc>>> &cub) {
    auto min_val = std::numeric_limits<Tp>::max();
    for (const auto &mat : cub) min_val = std::min(min_val, MinValue(mat));
    return min_val;
}

template<typename Tp, typename Alloc>
void Append(
    std::vector<Tp, Alloc> &vec1,
    const std::vector<Tp, Alloc> &vec2) {
    vec1.insert(vec1.end(), vec2.begin(), vec2.end());
}

template<typename Tp, typename Alloc>
[[nodiscard]] std::vector<Tp, Alloc> Concatenate(
    const std::vector<Tp, Alloc> &vec1,
    const std::vector<Tp, Alloc> &vec2) {
    std::vector<Tp, Alloc> out(vec1.begin(), vec1.end());
    out.insert(out.end(), vec2.begin(), vec2.end());
    return out;
}

template<typename Tp, typename Alloc = std::allocator<Tp>>
[[nodiscard]] std::vector<Tp, Alloc> Range(Tp min, Tp max) {
    assert(std::numeric_limits<Tp>::is_integer);
    assert(min < max);
    std::vector<Tp, Alloc> out(max - min);
    std::iota(out.begin(), out.end(), min);
    return out;
}

template<typename Tp, typename Alloc = std::allocator<Tp>>
[[nodiscard]] std::vector<Tp, Alloc> Range(Tp max) {
    return Range(Tp(0), max);
}

template<typename Tp, typename Ti>
[[nodiscard]] std::vector<Tp> Select(const std::vector<Tp> &vec, const std::vector<Ti> &indices) {
    // const auto &m = vec.size();
    const auto &n = indices.size();
    std::vector<Tp> out;
    out.reserve(n);
    for (const auto &i : indices) {
        // assert(0 <= i < m);
        out.emplace_back(vec[i]);
    }
    return out;
}

template<typename Tp, typename Ti>
[[nodiscard]] std::vector<std::vector<Tp>> Select(
    const std::vector<std::vector<Tp>> &vec,
    const std::vector<Ti> &indices
) {
    // const auto &m = vec.size();
    const auto &n = indices.size();
    std::vector<std::vector<Tp>> out;
    out.reserve(n);
    for (const auto &i : indices) {
        // assert(0 <= i < m);
        out.emplace_back(vec[i].begin(), vec[i].end());
    }
    return out;
}

template<typename Tp, typename Ti>
[[nodiscard]] std::vector<Tp> Reorder(const std::vector<Tp> &vec, const std::vector<Ti> &indices) {
    assert(indices.size() == vec.size());
    return Select(vec, indices);
}

}

#endif //TRIVIS_UTILS_GENERIC_UTILS_H_
