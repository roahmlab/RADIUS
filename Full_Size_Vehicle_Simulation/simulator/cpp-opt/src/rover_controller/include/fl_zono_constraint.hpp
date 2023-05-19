#ifndef ROAHM_FL_ZONO_CONSTRAINT_HPP_
#define ROAHM_FL_ZONO_CONSTRAINT_HPP_

#include "fl_zono_obs_set.hpp"
#include "gencon.hpp"
#include "vehrs.hpp"
#include <fmt/format.h>
#include <limits>  // for numeric_limits
#include <utility> // for move
#include <vector>  // for vector

namespace roahm::fl_zono_constraint {

class FlZonoConstraint {
private:
  // For generation
  const std::vector<double> a_full_;
  const std::vector<double> b_full_;
  const std::vector<int> constraint_start_indices_;
  const std::vector<int> constraint_sizes_;
  const double slice_center_val_;
  const double slice_generator_val_;

  // Computed / cached values
  double prev_evaluation_value_;
  std::vector<double> constraint_evals_;
  std::vector<double> gradient_evals_;

  void EvaluateAt(double val);

public:
  FlZonoConstraint(std::vector<double> a_full, std::vector<double> b_full,
                   std::vector<int> constraint_start_indices,
                   std::vector<int> constraint_sizes, double slice_center_val,
                   double slice_generator_val);

  [[nodiscard]] int NumConstraints() const;

  void GetConstraintEvaluations(double val, double* const out,
                                const int num_expected_vals);

  void GetGradient(double val, double* const out, const int num_expected_vals);

  [[nodiscard]] inline bool ConstraintsSatisfied(double val) {
    std::vector<double> curr_cons;
    curr_cons.resize(NumConstraints(), 0.0);
    GetConstraintEvaluations(val, curr_cons.data(), NumConstraints());
    for (const auto& val : curr_cons) {
      if (val > 0.0) {
        return false;
      }
    }
    return true;
  }
};

static inline ::roahm::fl_zono_constraint::FlZonoConstraint
GetFlZonoConstraint(const Vehrs& vehrs, const FlZonoObsSet& fl_zono_obs) {
  auto old_classical_constraints_ = GenerateConstraints(vehrs, fl_zono_obs);
  std::vector<double> cons_a_mat;
  std::vector<double> cons_b_mat;
  std::vector<int> cons_start_idxs;
  std::vector<int> cons_sizes;
  double slc_center{};
  double slc_gen{};
  if (fl_zono_obs.GetNumObs() > 0) {
    for (int i = 0; i < old_classical_constraints_.zono_startpoints_.size();
         ++i) {
      cons_start_idxs.push_back(
          old_classical_constraints_.zono_startpoints_.at(i));
      cons_sizes.push_back(old_classical_constraints_.zono_obs_sizes_.at(i));
    }
    for (int i = 0; i < old_classical_constraints_.num_b_elts_; ++i) {
      cons_a_mat.push_back(old_classical_constraints_.a_con_arr_[i]);
      cons_b_mat.push_back(old_classical_constraints_.b_con_arr_[i]);
    }
    const auto slc_cen_gen = vehrs.GetCenterGenTrajParam();
    slc_center = slc_cen_gen.first;
    slc_gen = slc_cen_gen.second;
  }
  return ::roahm::fl_zono_constraint::FlZonoConstraint{
      cons_a_mat, cons_b_mat, cons_start_idxs, cons_sizes, slc_center, slc_gen};
}

}; // namespace roahm::fl_zono_constraint

#endif // ROAHM_FL_ZONO_CONSTRAINT_HPP_