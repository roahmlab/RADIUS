#include "process_utils.hpp"

#include <gtest/gtest.h>

TEST(MoveColsToEnd, MoveMidEltOnly) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  const std::vector<int> cols_to_move{3};
  ::roahm::MoveColsToEnd(mat, cols_to_move);

  Eigen::Matrix<int, 1, 5> mat_expected_3;
  mat_expected_3 << 0, 1, 2, 4, 3;
  EXPECT_EQ(mat, mat_expected_3);
}
TEST(MoveColsToEnd, MoveLastColOnly) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  const std::vector<int> cols_to_move{4};
  ::roahm::MoveColsToEnd(mat, cols_to_move);

  Eigen::Matrix<int, 1, 5> mat_expected;
  mat_expected << 0, 1, 2, 3, 4;
  EXPECT_EQ(mat, mat_expected);
}

TEST(MoveColsToEnd, MoveFirstColOnly) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  const std::vector<int> cols_to_move{0};
  ::roahm::MoveColsToEnd(mat, cols_to_move);

  Eigen::Matrix<int, 1, 5> mat_expected;
  mat_expected << 1, 2, 3, 4, 0;
  EXPECT_EQ(mat, mat_expected);
}

TEST(MoveColsToEnd, MultipleColsInOrder) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  Eigen::Matrix<int, 1, 5> mat_expected;
  ::roahm::MoveColsToEnd(mat, {1, 2, 3});
  mat_expected << 0, 4, 1, 2, 3;
  EXPECT_EQ(mat, mat_expected);
}

TEST(MoveColsToEnd, MultipleColsOutOfOrder) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  Eigen::Matrix<int, 1, 5> mat_expected_32;
  ::roahm::MoveColsToEnd(mat, {3, 2});
  mat_expected_32 << 0, 1, 4, 3, 2;
  EXPECT_EQ(mat, mat_expected_32);
}

TEST(MoveColsToEnd, Subset) {
  Eigen::Matrix<int, 1, 5> mat;
  mat << 0, 1, 2, 3, 4;
  const std::vector<int> cols_to_move{0, 1};
  auto inner_portion = mat.rightCols(4).leftCols(3);
  ::roahm::MoveColsToEnd(inner_portion, cols_to_move);
  Eigen::Matrix<int, 1, 5> mat_expected;
  mat_expected << 0, 3, 1, 2, 4;
  EXPECT_EQ(mat, mat_expected);
}