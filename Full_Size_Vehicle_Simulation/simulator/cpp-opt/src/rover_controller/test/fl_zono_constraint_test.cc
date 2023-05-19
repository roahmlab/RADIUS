#include "fl_zono_constraint.hpp"

#include <gtest/gtest.h>

template <typename Container>
Container DivideContainer(Container container, double divide_by) {
  for (auto& val : container) {
    val /= divide_by;
  }
  return container;
}

TEST(FlZonoConstraintTest, Small) {
  const std::vector<double> a_mat{5.0, 10.0, 15.0, 20.0};
  const std::vector<double> b_mat{1.0, 3.0, 4.0, 2.0};
  const std::vector<int> constraint_start_indices{0, 2};
  const std::vector<int> constraint_sizes{2, 2};
  const double slice_center_val = 1.0;
  const double slice_generator_val = 2.0;
  std::array<double, 2> constraint_vals;
  std::array<double, 2> gradient_vals;

  {
    const double eval_pt = slice_center_val;
    // Evaluate at lambda = 0
    // -(A lambda - b) = b = {1, 3, 4, 2}
    // min(-(A_i lambda - b_i))
    // => constraints = {  1,   2 }
    //    gradients   = { -5, -20 } / slice_gen

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    const auto expected_grad{DivideContainer(std::array<double, 2>{-5.0, -20.0},
                                             slice_generator_val)};
    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 2);
    EXPECT_EQ(constraint_vals.at(0), 1.0);
    EXPECT_EQ(constraint_vals.at(1), 2.0);
    constraint.GetGradient(eval_pt, gradient_vals.data(), 2);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
  }

  {
    const double eval_pt = slice_center_val + slice_generator_val;
    // Evaluate at lambda = 1
    // -(A lambda - b) = {4, -7, -11, -18}
    // min(-(A_i lambda - b_i))
    // => constraints = {  -7, -18 }
    //    gradients   = { -10, -20 } / slice_gen

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    const auto expected_grad{DivideContainer(
        std::array<double, 2>{-10.0, -20.0}, slice_generator_val)};
    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 2);
    EXPECT_EQ(constraint_vals.at(0), -7.0);
    EXPECT_EQ(constraint_vals.at(1), -18.0);
    constraint.GetGradient(eval_pt, gradient_vals.data(), 2);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
  }

  {
    const double eval_pt = slice_center_val - slice_generator_val;
    // Evaluate at lambda = -1
    // -(A lambda - b) = {6, 13, 19, 22}
    // min(-(A_i lambda - b_i))
    // => constraints = {  6,  19 }
    //    gradients   = { -5, -15 } / slice_gen

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    const auto expected_grad{DivideContainer(std::array<double, 2>{-5.0, -15.0},
                                             slice_generator_val)};
    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 2);
    EXPECT_EQ(constraint_vals.at(0), 6.0);
    EXPECT_EQ(constraint_vals.at(1), 19.0);
    constraint.GetGradient(eval_pt, gradient_vals.data(), 2);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
  }
}

TEST(FlZonoConstraintTest, Medium) {
  const std::vector<double> a_mat{5.0,  10.1, -15.0, 20.0,
                                  18.3, 19.7, -16.8, 21.6};
  const std::vector<double> b_mat{1.0,   -3.0, 4.0,  2.0,
                                  -14.0, 25.3, 94.5, -12.0};
  const std::vector<int> constraint_start_indices{0, 1, 2, 5};
  const std::vector<int> constraint_sizes{1, 1, 3, 3};
  const double slice_center_val = -1.5;
  const double slice_generator_val = 2.5;

  std::array<double, 4> constraint_vals;
  std::array<double, 4> gradient_vals;

  {
    const double eval_pt = slice_center_val;
    // Evaluate at lambda = 0

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 4);
    EXPECT_EQ(constraint_vals.at(0), 1.0);
    EXPECT_EQ(constraint_vals.at(1), -3.0);
    EXPECT_EQ(constraint_vals.at(2), -14.0);
    EXPECT_EQ(constraint_vals.at(3), -12.0);
    const std::array<double, 4> expected_grad{DivideContainer(
        std::array<double, 4>{-5.0, -10.1, -18.3, -21.6}, slice_generator_val)};
    constraint.GetGradient(eval_pt, gradient_vals.data(), 4);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
    EXPECT_DOUBLE_EQ(gradient_vals.at(2), expected_grad.at(2));
    EXPECT_DOUBLE_EQ(gradient_vals.at(3), expected_grad.at(3));
  }

  {
    const double eval_pt = slice_center_val + slice_generator_val;
    // Evaluate at lambda = 1

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    const std::array<double, 4> expected_cons{-4.0, -13.1, -32.3, -33.6};
    const std::array<double, 4> expected_grad{DivideContainer(
        std::array<double, 4>{-5.0, -10.1, -18.3, -21.6}, slice_generator_val)};
    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 4);
    EXPECT_EQ(constraint_vals.at(0), expected_cons.at(0));
    EXPECT_EQ(constraint_vals.at(1), expected_cons.at(1));
    EXPECT_EQ(constraint_vals.at(2), expected_cons.at(2));
    EXPECT_EQ(constraint_vals.at(3), expected_cons.at(3));
    constraint.GetGradient(eval_pt, gradient_vals.data(), 4);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
    EXPECT_DOUBLE_EQ(gradient_vals.at(2), expected_grad.at(2));
    EXPECT_DOUBLE_EQ(gradient_vals.at(3), expected_grad.at(3));
  }

  {
    const double eval_pt = slice_center_val - slice_generator_val;
    // Evaluate at lambda = -1

    ::roahm::fl_zono_constraint::FlZonoConstraint constraint(
        a_mat, b_mat, constraint_start_indices, constraint_sizes,
        slice_center_val, slice_generator_val);
    const std::array<double, 4> expected_cons{6.0, 7.1, -11.0, 9.6};
    const std::array<double, 4> expected_grad{DivideContainer(
        std::array<double, 4>{-5.0, -10.1, 15.0, -21.6}, slice_generator_val)};

    constraint.GetConstraintEvaluations(eval_pt, constraint_vals.data(), 4);
    EXPECT_DOUBLE_EQ(constraint_vals.at(0), expected_cons.at(0));
    EXPECT_DOUBLE_EQ(constraint_vals.at(1), expected_cons.at(1));
    EXPECT_DOUBLE_EQ(constraint_vals.at(2), expected_cons.at(2));
    EXPECT_DOUBLE_EQ(constraint_vals.at(3), expected_cons.at(3));
    constraint.GetGradient(eval_pt, gradient_vals.data(), 4);
    EXPECT_DOUBLE_EQ(gradient_vals.at(0), expected_grad.at(0));
    EXPECT_DOUBLE_EQ(gradient_vals.at(1), expected_grad.at(1));
    EXPECT_DOUBLE_EQ(gradient_vals.at(2), expected_grad.at(2));
    EXPECT_DOUBLE_EQ(gradient_vals.at(3), expected_grad.at(3));
  }
}
