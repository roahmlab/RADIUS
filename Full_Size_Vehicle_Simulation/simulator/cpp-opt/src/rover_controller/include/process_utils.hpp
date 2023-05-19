#ifndef ROAHM_PROCESS_UTILS_HPP_
#define ROAHM_PROCESS_UTILS_HPP_

#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <numeric>
#include <vector>

namespace roahm {

namespace {
template <typename T, typename V>
void SortRelative(T& container_to_sort,
                  const V& container_to_sort_relative_to) {
  if (container_to_sort.size() != container_to_sort_relative_to.size()) {
    throw std::runtime_error(
        "container_to_sort and container_to_sort_relative_to must be the same "
        "size");
  }
  std::vector<int> indices(container_to_sort.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(
      indices.begin(), indices.end(),
      [&container_to_sort_relative_to, indices](const int& a, const int& b) {
        return container_to_sort_relative_to.at(a) <
               container_to_sort_relative_to.at(b);
      });
  T container_to_sort_copy;
  for (const auto& idx : indices) {
    container_to_sort_copy.push_back(container_to_sort.at(idx));
  }
  container_to_sort = container_to_sort_copy;
}

template <typename RandomAccessContainer>
std::vector<int> ArgSort(RandomAccessContainer& vec) {
  std::vector<int> idx(vec.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&vec](int i1, int i2) { return vec.at(i1) < vec.at(i2); });
  return idx;
}
} // namespace

template <typename T>
void MoveColsToEnd(T& mat, std::vector<int> cols_to_move_idxs) {
  if (cols_to_move_idxs.size() > mat.cols()) {
    return;
  }
  auto sorted_move_indices = ArgSort(cols_to_move_idxs);
  std::reverse(sorted_move_indices.begin(), sorted_move_indices.end());
  std::vector<std::remove_cv_t<decltype(mat.col(0).eval())>> cols_to_move;
  for (const auto& col_move_idx_idx : sorted_move_indices) {
    const auto col_move_idx = cols_to_move_idxs.at(col_move_idx_idx);
    const auto col_move_val = mat.col(col_move_idx).eval();
    const auto num_cols_to_copy = mat.cols() - col_move_idx - 1;
    const auto copy_block = mat.rightCols(num_cols_to_copy).eval();
    mat.leftCols(mat.cols() - 1).rightCols(num_cols_to_copy) = copy_block;
    cols_to_move.push_back(col_move_val);
  }
  SortRelative(cols_to_move, sorted_move_indices);
  for (int i = 0; i < cols_to_move.size(); ++i) {
    const int final_col_idx = mat.cols() - cols_to_move.size() + i;
    mat.col(final_col_idx) = cols_to_move.at(i);
  }
}
} // namespace roahm
#endif // ROAHM_PROCESS_UTILS_HPP_