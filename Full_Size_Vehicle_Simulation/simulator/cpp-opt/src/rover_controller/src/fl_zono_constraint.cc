#include "fl_zono_constraint.hpp"

#include <cassert>   // for assert
#include <stdexcept> // for runtime_error
#include <string>    // for operator+, char_traits, to_string

namespace roahm::fl_zono_constraint {
FlZonoConstraint::FlZonoConstraint(std::vector<double> a_full,
                                   std::vector<double> b_full,
                                   std::vector<int> constraint_start_indices,
                                   std::vector<int> constraint_sizes,
                                   double slice_center_val,
                                   double slice_generator_val)
    : a_full_{std::move(a_full)}, b_full_{std::move(b_full)},
      constraint_start_indices_{std::move(constraint_start_indices)},
      constraint_sizes_{std::move(constraint_sizes)},
      slice_center_val_{slice_center_val},
      slice_generator_val_{slice_generator_val},
      prev_evaluation_value_{
          std::numeric_limits<decltype(prev_evaluation_value_)>::quiet_NaN()},
      constraint_evals_{}, gradient_evals_{} {}

[[nodiscard]] int FlZonoConstraint::NumConstraints() const {
  return constraint_sizes_.size();
}

void FlZonoConstraint::EvaluateAt(const double val) {
  if (prev_evaluation_value_ == val) {
    assert(constraint_sizes_.size() == constraint_start_indices_.size());
    assert(gradient_evals_.size() == constraint_start_indices_.size());
    assert(constraint_evals_.size() == constraint_start_indices_.size());
    return;
  }
  const double lambda = (val - slice_center_val_) / slice_generator_val_;

  if (constraint_evals_.size() != constraint_start_indices_.size()) {
    constraint_evals_.resize(constraint_start_indices_.size());
  }
  if (gradient_evals_.size() != constraint_start_indices_.size()) {
    gradient_evals_.resize(constraint_start_indices_.size());
  }

  // Check preconditions
  assert(constraint_sizes_.size() == constraint_start_indices_.size());
  assert(gradient_evals_.size() == constraint_start_indices_.size());
  assert(constraint_evals_.size() == constraint_start_indices_.size());

  for (int i = 0; i < constraint_start_indices_.size(); ++i) {
    const auto row_start = constraint_start_indices_.at(i);
    const auto num_rows = constraint_sizes_.at(i);
    // Evaluate min(Ax - b)
    int min_row = -1;
    double min_val = std::numeric_limits<decltype(min_val)>::max();
    assert(num_rows > 0);
    assert(row_start >= 0);
    assert(row_start + num_rows <= a_full_.size());
    for (int row = row_start; row < row_start + num_rows; ++row) {
      const double row_eval = -(a_full_.at(row) * lambda - b_full_.at(row));
      if (row_eval < min_val) {
        min_val = row_eval;
        min_row = row;
      }
    }

    constraint_evals_.at(i) = min_val;
    gradient_evals_.at(i) = -a_full_.at(min_row) / slice_generator_val_;
  }
}

void FlZonoConstraint::GetConstraintEvaluations(double val, double* const out,
                                                const int num_expected_vals) {
  EvaluateAt(val);
  if (num_expected_vals != constraint_evals_.size()) {
    throw std::runtime_error("Expected number of values [" +
                             std::to_string(num_expected_vals) +
                             "] != actual number of values [" +
                             std::to_string(constraint_evals_.size()) + "]");
  }
  for (int i = 0; i < constraint_evals_.size(); ++i) {
    out[i] = constraint_evals_[i];
  }
}

void FlZonoConstraint::GetGradient(double val, double* const out,
                                   const int num_expected_vals) {
  EvaluateAt(val);
  if (num_expected_vals != gradient_evals_.size()) {
    throw std::runtime_error("Expected number of values [" +
                             std::to_string(num_expected_vals) +
                             "] != actual number of values [" +
                             std::to_string(gradient_evals_.size()) + "]");
  }
  for (int i = 0; i < gradient_evals_.size(); ++i) {
    out[i] = gradient_evals_[i];
  }
}

} // namespace roahm::fl_zono_constraint
