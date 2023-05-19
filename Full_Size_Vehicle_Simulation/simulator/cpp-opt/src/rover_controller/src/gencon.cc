#include "gencon.hpp"

#include <algorithm>        // for max, minmax_element
#include <array>            // for array, array<>::value...
#include <cassert>          // for assert
#include <cmath>            // for atan2, cos, sin, sqrt
#include <cstdlib>          // for abs, size_t
#include <fstream>          // for ifstream, ofstream
#include <initializer_list> // for initializer_list
#include <iostream>         // for operator<<, basic_ost...
#include <utility>          // for pair, make_pair

#include "simple_util.hpp" // for ToAngle, Square, Angl...
#include <fmt/format.h>

namespace roahm {

namespace {

/// Takes a pair \f$ (a, b) \f$ and returns a normalized pair,
/// \f$ (\frac{a}{\sqrt{a^2 + b^2}}, \frac{b}{\sqrt{a^2 + b^2}}) \f$
/// \param a the first number of the set
/// \param b the second number of the set
/// \return the inputs normalized
/// \f$ (\frac{a}{\sqrt{a^2 + b^2}}, \frac{b}{\sqrt{a^2 + b^2}}) \f$
inline std::pair<double, double> NormPair(double a, double b) {
  const double norm_recip = 1.0 / Norm(a, b);
  return std::make_pair(a * norm_recip, b * norm_recip);
}

// TODO make sure that this is necessary, compare with constraint generation
/// Fills arrays with zonotope information from FRSes, leaving room for obstacle
/// generators to be appended for constraint generation
/// \param vehrs the forward reachable set
/// \return filled arrays containing the reachable sets' zonotope information
/// and space to append obstacle generators for constraint generation
ZonoInfo FillZonosFromVehrs(const Vehrs& vehrs) {
  // part of gen constraints (maybe)

  // SUS

  ZonoInfo ret;
  // Number of total input zonotopes (i.e. the "size" of the FRS)
  IndexT num_zonos = vehrs.GetNumZonos();
  ret.num_zonos_ = num_zonos;

  // Number of output constraints for each input zonotope, per obstacle
  ret.num_out_zono_gens_arr_ = std::make_unique<IndexT[]>(num_zonos);
  auto& num_out_zono_gens_arr = ret.num_out_zono_gens_arr_;
  ret.cum_size_arr_ = std::make_unique<IndexT[]>(num_zonos);
  auto& cum_size_arr = ret.cum_size_arr_;

  for (IndexT i = 0; i < num_zonos; ++i) {
    // Since we always have two generators per obstacles, we leave room for them
    num_out_zono_gens_arr[i] = vehrs.zono_sizes_.at(i) + 2;
  }

  // Compute total sizes
  IndexT cumsum_num_out_zono_gens = 0;
  for (IndexT i = 0; i < num_zonos; ++i) {
    cum_size_arr[i] = cumsum_num_out_zono_gens;
    cumsum_num_out_zono_gens += num_out_zono_gens_arr[i];
  }
  ret.cumsum_num_out_zono_gens_ = cumsum_num_out_zono_gens;

  // Each generator has x, y so multiply by 2 for size
  IndexT zono_arr_size = 2 * cumsum_num_out_zono_gens;

  // Create arrays
  ret.zono_arr_ = std::make_unique<double[]>(zono_arr_size);
  auto& zono_arr = ret.zono_arr_;

  for (IndexT zono_idx = 0; zono_idx < num_zonos; ++zono_idx) {
    const IndexT out_start_idx = 2 * cum_size_arr[zono_idx];
    const IndexT num_out_gens = num_out_zono_gens_arr[zono_idx];
    const IndexT num_out_elts = num_out_gens * 2;
    const IndexT num_in_gens = num_out_gens - 2;
    const IndexT num_in_elts = num_in_gens * 2;
    const IndexT in_start_idx = out_start_idx - 4 * zono_idx;

    // Fill the forward reachable set zonotopes
    for (IndexT i = 0; i < num_in_elts; ++i) {
      zono_arr[out_start_idx + i] = vehrs.zono_xy_.at(in_start_idx + i);
    }
    double x_sum = 0.0;
    double y_sum = 0.0;
    for (IndexT i = 0; i < (num_in_elts / 2); ++i) {
      const double x = vehrs.zono_xy_.at(in_start_idx + (2 * i) + 0);
      const double y = vehrs.zono_xy_.at(in_start_idx + (2 * i) + 1);
      x_sum += x;
      y_sum += y;
    }

    // Leave the obstacle generators empty for now
    for (IndexT i = num_in_elts; i < num_out_elts; ++i) {
      zono_arr[out_start_idx + i] = 0;
    }
  }
  return ret;
}

} // namespace

class PHOut {
public:
  Eigen::Matrix<double, Eigen::Dynamic, 2> a_mat_;
  Eigen::VectorXd b_mat_;
  Eigen::VectorXd idx_mat_;
};

PHOut PolytopeHalfspace(Eigen::Vector2d center,
                        Eigen::Matrix<double, 2, Eigen::Dynamic> gens) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> c;
  c.conservativeResize(2, gens.cols());
  c.row(0) = -gens.row(1);
  c.row(1) = gens.row(0);
  c.colwise().normalize();
  c.transposeInPlace();
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_d =
      ((((c * gens).transpose()).array().abs()).colwise().sum()).transpose();
  Eigen::Matrix<double, Eigen::Dynamic, 1> d = c * center;

  PHOut cons;
  auto& a_mat = cons.a_mat_;
  a_mat.conservativeResize(gens.cols() * 2, 2);
  a_mat << c, -c;
  auto& b_mat = cons.b_mat_;
  b_mat.conservativeResize(gens.cols() * 2, 1);
  auto b_mat1 = delta_d + d;
  auto b_mat2 = delta_d - d;
  b_mat << b_mat1, b_mat2;
  return cons;
}
inline int PolytopeHalfspaceInPlace(
    const Eigen::Vector2d& center,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& zono_gens,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& zono_gen_t,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& delta_d_ul_submat,
    const Eigen::Matrix<double, 2, 2>& obs_gens,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& a_mat_out,
    Eigen::VectorXd& b_mat_out, int row_out_start,
    const Eigen::Vector2d& slc_xy) {
  const double obs_x1_un = obs_gens(0, 0);
  const double obs_x2_un = obs_gens(0, 1);
  const double obs_y1_un = obs_gens(1, 0);
  const double obs_y2_un = obs_gens(1, 1);
  const double obs_norm_1 =
      std::sqrt(obs_x1_un * obs_x1_un + obs_y1_un * obs_y1_un);
  const double obs_norm_2 =
      std::sqrt(obs_x2_un * obs_x2_un + obs_y2_un * obs_y2_un);
  const double obs_x1_n = obs_x1_un / obs_norm_1;
  const double obs_y1_n = obs_y1_un / obs_norm_1;
  const double obs_x2_n = obs_x2_un / obs_norm_2;
  const double obs_y2_n = obs_y2_un / obs_norm_2;
  zono_gen_t(zono_gen_t.rows() - 2, 0) = -obs_y1_n;
  zono_gen_t(zono_gen_t.rows() - 2, 1) = obs_x1_n;
  zono_gen_t(zono_gen_t.rows() - 1, 0) = -obs_y2_n;
  zono_gen_t(zono_gen_t.rows() - 1, 1) = obs_x2_n;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_d_bl_submat =
      ((zono_gen_t.block(zono_gen_t.rows() - 2, 0, 2, 2) * zono_gens)
           .array()
           .abs())
          .rowwise()
          .sum();
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& c = zono_gen_t;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_d =
      ((zono_gen_t * obs_gens).array().abs()).rowwise().sum();
  delta_d.block(0, 0, delta_d.rows() - 2, 1) += delta_d_ul_submat;
  delta_d.block(delta_d.rows() - 2, 0, 2, 1) += delta_d_bl_submat;
  Eigen::Matrix<double, Eigen::Dynamic, 1> d = c * center;
  auto c_times_slc_xy = c * slc_xy;
  const int new_rows = 2 * c.rows();
  a_mat_out.block(row_out_start, 0, new_rows, 1) << c_times_slc_xy,
      -c_times_slc_xy;
  b_mat_out.block(row_out_start, 0, new_rows, 1) << delta_d + d, delta_d - d;
  return new_rows;
}

struct Gencon2Zono {
  Eigen::Vector2d center_xy_;
  Eigen::Matrix<double, 2, Eigen::Dynamic> generators_xy_;
};

// Just a (currently very slow), but correct Eigen implementation for REFINE
// constraints. May come in handy as a ground truth or if we do REFINE/RISK
// constraint selection on a (zono, obs) pair basis instead of max vs.
// the full FRS.
Gencon2Out GenerateConstraints2(const Vehrs& reach_set,
                                const FlZonoObsSet& local_obs) {
  Gencon2Out cons;
  const int n_obs = local_obs.GetNumObs();
  const int n_reach = reach_set.GetNumZonos();
  {
    int prev_start = 0;
    int prev_size = 0;
    for (int reach_idx = 0; reach_idx < n_reach; ++reach_idx) {
      for (int obs_idx = 0; obs_idx < n_obs; ++obs_idx) {
        const int num_reach_zono_gens = 2;
        const int next_size = (reach_set.zono_sizes_.at(reach_idx) + 2) * 2;
        const int next_start = prev_start + prev_size;

        cons.start_idx_mat_.conservativeResize(cons.start_idx_mat_.rows() + 1,
                                               1);
        cons.start_idx_mat_(cons.start_idx_mat_.rows() - 1, 0) = next_start;
        cons.size_mat_.conservativeResize(cons.size_mat_.rows() + 1, 1);
        cons.size_mat_(cons.size_mat_.rows() - 1, 0) = next_size;

        prev_size = next_size;
        prev_start = next_start;
      }
    }
    cons.a_mat_.conservativeResize(prev_start + prev_size, 1);
    cons.b_mat_.conservativeResize(prev_start + prev_size, 1);
  }

  int prev_start = 0;
  for (int reach_idx = 0; reach_idx < n_reach; ++reach_idx) {
    // Conversion
    const double zono_center_x = reach_set.xy_centers_.at(reach_idx * 2);
    const double zono_center_y = reach_set.xy_centers_.at(reach_idx * 2 + 1);
    const Eigen::Vector2d zono_center_xy{zono_center_x, zono_center_y};
    Eigen::Matrix<double, 2, Eigen::Dynamic> zono_gens;
    int zono_gen_start = 0;
    for (int i = 0; i < reach_idx; ++i) {
      zono_gen_start += 2 * reach_set.zono_sizes_.at(i);
    }
    const int n_gens_curr = reach_set.zono_sizes_.at(reach_idx);
    zono_gens.conservativeResize(2, n_gens_curr);
    for (int gen_idx = 0; gen_idx < n_gens_curr; ++gen_idx) {
      zono_gens(0, gen_idx) =
          reach_set.zono_xy_.at(zono_gen_start + (gen_idx * 2));
      zono_gens(1, gen_idx) =
          reach_set.zono_xy_.at(zono_gen_start + (gen_idx * 2) + 1);
    }

    Gencon2Zono zono{zono_center_xy, zono_gens};
    const double zono_slc_x = reach_set.slc_x_.at(reach_idx);
    const double zono_slc_y = reach_set.slc_y_.at(reach_idx);
    const Eigen::Vector2d slc_xy{zono_slc_x, zono_slc_y};

    const auto zono_time_interval =
        reach_set.zono_time_intervals_.at(reach_idx);

    Eigen::Matrix<double, Eigen::Dynamic, 2> zono_gen_t;
    zono_gen_t.conservativeResize(zono_gens.cols() + 2, 2);
    zono_gen_t.col(0) << -zono_gens.row(1).transpose(), 0, 0;
    zono_gen_t.col(1) << zono_gens.row(0).transpose(), 0, 0;
    zono_gen_t.rowwise().normalize();
    Eigen::Matrix<double, Eigen::Dynamic, 1> delta_d_ul_submat =
        ((zono_gen_t.block(0, 0, zono_gen_t.rows() - 2, 2) * zono_gens)
             .array()
             .abs())
            .rowwise()
            .sum();

    for (int obs_idx = 0; obs_idx < n_obs; ++obs_idx) {
      // Conversion
      const auto obs_center_pre =
          local_obs.GetCenter(obs_idx, zono_time_interval);
      const Eigen::Vector2d obs_center{obs_center_pre.x_, obs_center_pre.y_};
      const auto obs_gen_pre =
          local_obs.GetGenerators(obs_idx, zono_time_interval);

      Eigen::Matrix<double, 2, 2> obs_gens;
      obs_gens(0, 0) = obs_gen_pre.x1_;
      obs_gens(0, 1) = obs_gen_pre.x2_;
      obs_gens(1, 0) = obs_gen_pre.y1_;
      obs_gens(1, 1) = obs_gen_pre.y2_;
      const Eigen::Vector2d center_ph = obs_center - zono_center_xy;

      // Cons gen
      // TODO should be num_decision_vars instead of 1
      const int row_start = prev_start;
      const int new_rows = PolytopeHalfspaceInPlace(
          center_ph, zono_gens, zono_gen_t, delta_d_ul_submat, obs_gens,
          cons.a_mat_, cons.b_mat_, row_start, slc_xy);
      prev_start = row_start + new_rows;
    }
  }
  return cons;
}

Constraints GenerateConstraints(const Vehrs& vehrs,
                                const FlZonoObsSet& obs_info) {
  // SUS
  // TODO: get rid of shared and unique pointers and replace them with vectors

  const IndexT num_obs_gens = 2;
  const auto filled_test_zonos = FillZonosFromVehrs(vehrs);
  const auto num_obs = obs_info.GetNumObs();

  const auto num_zonos = filled_test_zonos.num_zonos_;
  const auto& cumsum_num_out_zono_gens =
      filled_test_zonos.cumsum_num_out_zono_gens_;
  const auto& num_out_zono_gens_arr = filled_test_zonos.num_out_zono_gens_arr_;
  const auto& cum_size_arr = filled_test_zonos.cum_size_arr_;
  const auto& zono_arr = filled_test_zonos.zono_arr_;
  if (num_obs < 1 || num_zonos < 1) {
    return Constraints{};
  }

  const IndexT total_d_size = cumsum_num_out_zono_gens;
  const IndexT total_c_size = 2 * cumsum_num_out_zono_gens;
  double* const delta_d_arr = new double[total_d_size];
  double* const c_arr = new double[total_c_size];
  for (IndexT i = 0; i < total_d_size; ++i) {
    delta_d_arr[i] = 0;
  }

  // Compute initial values.
  for (IndexT zono_idx = 0; zono_idx < num_zonos; ++zono_idx) {
    // Number of output generators for current zonotope
    const IndexT curr_n = num_out_zono_gens_arr[zono_idx];

    // Offsets for aggregated array output
    const IndexT n_start_offset = cum_size_arr[zono_idx];
    const IndexT d_start_idx = n_start_offset;
    const IndexT zono_start_idx = 2 * n_start_offset;
    const IndexT c_start_idx = 2 * n_start_offset;

    double* const curr_delta_d_arr = delta_d_arr + d_start_idx;
    double* const curr_c_arr = c_arr + c_start_idx;
    double* const curr_zono_arr = zono_arr.get() + zono_start_idx;

    const IndexT max_r_without_obs = curr_n - num_obs_gens;
    // delta_d(r, c) = abs(C(r,:) * G(:,c))
    for (IndexT r = 0; r < max_r_without_obs; ++r) {
      curr_delta_d_arr[r] = 0;

      const auto g_r0 = curr_zono_arr[r * 2];
      const auto g_r1 = curr_zono_arr[r * 2 + 1];
      // C(r,:) = Normalize([-G(1, r); G(2, r)]), one-indexed
      const auto [c_r0, c_r1] =
          NormPair(-curr_zono_arr[r * 2 + 1], curr_zono_arr[r * 2]);

      // Save to C array
      curr_c_arr[r * 2] = c_r0;
      curr_c_arr[r * 2 + 1] = c_r1;

      for (IndexT c = 0; c < curr_n - num_obs_gens; ++c) {
        const double g_c0 = curr_zono_arr[c * 2];
        const double g_c1 = curr_zono_arr[c * 2 + 1];
        curr_delta_d_arr[r] += std::abs((c_r0 * g_c0) + (c_r1 * g_c1));
      }
    }
  }

  // TODO allow for multiple slice
  const IndexT pa_size_per_obs = (2 * total_c_size); // PA = [-C; C]
  // Pb = [d + deltaD; -d + deltaD]
  const IndexT pb_size_per_obs = (2 * total_d_size);
  const IndexT total_pa_size = pa_size_per_obs * num_obs;
  const IndexT total_pb_size = pb_size_per_obs * num_obs;

  constexpr int kNumSlcGens = 1;
  static_assert(kNumSlcGens == 1); // Not implemented for other values yet
  const IndexT total_a_con_size = total_pb_size * kNumSlcGens;
  double* const d_arr = new double[total_d_size];
  double* const pa_arr = new double[total_pa_size];
  Constraints ret;
  ret.zono_startpoints_.push_back(0);
  for (IndexT zono_idx = 0; zono_idx < num_zonos; ++zono_idx) {
    for (IndexT obs_idx = 0; obs_idx < num_obs; ++obs_idx) {
      const auto prev_start = ret.zono_startpoints_.back();
      const auto additional_gens = 2 * num_out_zono_gens_arr[zono_idx];
      if (!(zono_idx == num_zonos - 1 && obs_idx == num_obs - 1)) {
        ret.zono_startpoints_.push_back(prev_start + additional_gens);
      }
      ret.zono_obs_sizes_.push_back(additional_gens);
    }
  }
  assert(ret.zono_obs_sizes_.size() == num_zonos * num_obs);
  assert(ret.zono_startpoints_.size() == num_zonos * num_obs);
  ret.num_a_elts_ = total_a_con_size;
  ret.num_b_elts_ = total_pb_size;
  ret.a_con_arr_ = std::make_unique<double[]>(total_a_con_size); // SUS
  ret.b_con_arr_ = std::make_unique<double[]>(total_pb_size);
  double* const pb_arr = ret.b_con_arr_.get();
  double* const a_con_arr = ret.a_con_arr_.get();

  // Initalize elements to zero
  for (IndexT i = 0; i < total_pa_size; ++i) {
    pa_arr[i] = 0;
  }
  for (IndexT i = 0; i < total_pb_size; ++i) {
    pb_arr[i] = 0;
  }

  for (IndexT zono_idx = 0; zono_idx < num_zonos; ++zono_idx) {
    const IndexT curr_n = num_out_zono_gens_arr[zono_idx];
    const IndexT n_start_offset = cum_size_arr[zono_idx];
    const IndexT c_start_idx = 2 * n_start_offset;
    const IndexT zono_start_idx = 2 * n_start_offset;
    const IndexT d_start_idx = n_start_offset;
    const IndexT max_r_without_obs = curr_n - num_obs_gens;
    assert(curr_n >= num_obs_gens);
    double* const curr_c_arr = c_arr + c_start_idx;
    double* const curr_zono_arr = zono_arr.get() + zono_start_idx;
    double* const curr_delta_d_arr = delta_d_arr + d_start_idx;

    const auto t_interval = vehrs.zono_time_intervals_.at(zono_idx);
    const double zono_center_x = vehrs.xy_centers_.at(zono_idx * 2);
    const double zono_center_y = vehrs.xy_centers_.at(zono_idx * 2 + 1);
    const double slc_x = vehrs.slc_x_.at(zono_idx);
    const double slc_y = vehrs.slc_y_.at(zono_idx);

    for (IndexT obs_idx = 0; obs_idx < num_obs; ++obs_idx) {
      const IndexT curr_row_out_start =
          (cum_size_arr[zono_idx] * num_obs) + (curr_n * obs_idx);

      // 2 for plus/minus components, 2 for obs dimensionality
      const IndexT pa_start_idx = 2 * 2 * curr_row_out_start;
      const IndexT pb_start_idx = 2 * curr_row_out_start;
      const auto [obs_center_x, obs_center_y] =
          obs_info.GetCenter(obs_idx, t_interval);
      double* const curr_pa_arr = pa_arr + pa_start_idx;
      double* const curr_pb_arr = pb_arr + pb_start_idx;

      // Copy the non-obstacle C portion
      for (IndexT r = 0; r < max_r_without_obs; ++r) {
        curr_pa_arr[r * 2] = curr_c_arr[r * 2];
        curr_pa_arr[r * 2 + 1] = curr_c_arr[r * 2 + 1];
        curr_pa_arr[(r + curr_n) * 2] = -curr_c_arr[r * 2];
        curr_pa_arr[(r + curr_n) * 2 + 1] = -curr_c_arr[r * 2 + 1];
      }

      // Copy the precomputed delta_d values
      for (IndexT r = 0; r < curr_n; ++r) {
        const double delta_d_val = curr_delta_d_arr[r];
        curr_pb_arr[r] = delta_d_val;
        curr_pb_arr[r + curr_n] = delta_d_val;
      }

      for (IndexT obs_gen_idx = 0; obs_gen_idx < num_obs_gens; ++obs_gen_idx) {
        const auto [obs_x, obs_y] =
            obs_info.GetSingleGenerator(obs_idx, obs_gen_idx, t_interval);
        const auto [c0_val, c1_val] = NormPair(-obs_y, obs_x);

        // Copy the normalized obstacle into C
        const IndexT pa_obs_start = (max_r_without_obs + obs_gen_idx) * 2;
        curr_pa_arr[pa_obs_start + 0] = c0_val;
        curr_pa_arr[pa_obs_start + 1] = c1_val;

        // Set the -C portion of PA
        curr_pa_arr[pa_obs_start + (curr_n * 2) + 0] = -c0_val;
        curr_pa_arr[pa_obs_start + (curr_n * 2) + 1] = -c1_val;
      }

      for (IndexT r = 0; r < curr_n - num_obs_gens; ++r) {
        for (IndexT obs_gen_idx = 0; obs_gen_idx < num_obs_gens;
             ++obs_gen_idx) {
          // In this, obs_gen_idx corresponds to the last num_obs_gens columns
          // TODO REMOVE const IndexT obs_gen_start_idx =
          // TODO REMOVE     (obs_idx * num_obs_gens * 2) + obs_gen_idx * 2;
          // TODO REMOVE const double obs_x = obs_arr.at(obs_gen_start_idx + 0);
          // TODO REMOVE const double obs_y = obs_arr.at(obs_gen_start_idx + 1);
          const auto [obs_x, obs_y] =
              obs_info.GetSingleGenerator(obs_idx, obs_gen_idx, t_interval);
          const double pa_val_x = curr_pa_arr[r * 2 + 0];
          const double pa_val_y = curr_pa_arr[r * 2 + 1];
          const double abs_cg_val =
              std::abs(pa_val_x * obs_x + pa_val_y * obs_y);
          curr_pb_arr[r] += abs_cg_val;
          curr_pb_arr[r + curr_n] += abs_cg_val;
        }
      }

      // Sum from the bottom num_obs_gen x num_obs_gen corner of abs(CG) matrix
      for (IndexT c_obs_gen_idx = 0; c_obs_gen_idx < num_obs_gens;
           ++c_obs_gen_idx) {
        const IndexT r = curr_n - (num_obs_gens - c_obs_gen_idx);
        const IndexT c_obs_gen_start_idx =
            (obs_idx * num_obs_gens * 2) + c_obs_gen_idx * 2;
        const auto [c_obs_x_unnorm, c_obs_y_unnorm] =
            obs_info.GetSingleGenerator(obs_idx, c_obs_gen_idx, t_interval);
        const auto [c_obs_x, c_obs_y] =
            NormPair(c_obs_x_unnorm, c_obs_y_unnorm);
        // TODO REMOVE const auto [c_obs_x, c_obs_y] = NormPair(
        // TODO REMOVE     obs_arr.at(c_obs_gen_start_idx),
        // obs_arr.at(c_obs_gen_start_idx + 1));
        for (IndexT g_obs_gen_idx = 0; g_obs_gen_idx < num_obs_gens;
             ++g_obs_gen_idx) {
          // TODO REMOVE const IndexT g_obs_gen_start_idx =
          // TODO REMOVE     (obs_idx * num_obs_gens * 2) + g_obs_gen_idx * 2;
          const auto [g_obs_x, g_obs_y] =
              obs_info.GetSingleGenerator(obs_idx, g_obs_gen_idx, t_interval);
          // TODO REMOVE const double g_obs_x = obs_arr.at(g_obs_gen_start_idx +
          // 0);
          // TODO REMOVE const double g_obs_y = obs_arr.at(g_obs_gen_start_idx +
          // 1);
          curr_pb_arr[r] += std::abs(-c_obs_y * g_obs_x + c_obs_x * g_obs_y);
        }
      }

      for (IndexT c_obs_gen_idx = 0; c_obs_gen_idx < num_obs_gens;
           ++c_obs_gen_idx) {
        const auto [c_obs_x_unnorm, c_obs_y_unnorm] =
            obs_info.GetSingleGenerator(obs_idx, c_obs_gen_idx, t_interval);
        const auto [c_obs_x, c_obs_y] =
            NormPair(c_obs_x_unnorm, c_obs_y_unnorm);
        const IndexT r = curr_n - (num_obs_gens - c_obs_gen_idx);
        // const IndexT c_obs_gen_start_idx =
        //     (obs_idx * num_obs_gens * 2) + c_obs_gen_idx * 2;
        //  TODO REMOVE const auto [c_obs_x, c_obs_y] = NormPair(
        //  TODO REMOVE     obs_arr.at(c_obs_gen_start_idx),
        //  obs_arr.at(c_obs_gen_start_idx + 1));
        for (IndexT c = 0; c < curr_n - num_obs_gens; ++c) {
          const double g_x = curr_zono_arr[c * 2 + 0];
          const double g_y = curr_zono_arr[c * 2 + 1];
          curr_pb_arr[r] += std::abs(-c_obs_y * g_x + c_obs_x * g_y);
        }
      }

      for (IndexT r = 0; r < curr_n; ++r) {
        curr_pb_arr[r + curr_n] = curr_pb_arr[r];
      }

      const double center_delta_x = obs_center_x - zono_center_x;
      const double center_delta_y = obs_center_y - zono_center_y;

      for (IndexT r = 0; r < curr_n; ++r) {
        const double d_val = curr_pa_arr[r * 2 + 0] * center_delta_x +
                             curr_pa_arr[r * 2 + 1] * center_delta_y;
        curr_pb_arr[r] += d_val;
        curr_pb_arr[r + curr_n] -= d_val;
        /// TODO make sure the offset by curr_n matches with constraint eval
        const double a_con_val =
            (curr_pa_arr[r * 2 + 0] * slc_x) + (curr_pa_arr[r * 2 + 1] * slc_y);
        a_con_arr[pb_start_idx + r] = a_con_val;
        a_con_arr[pb_start_idx + r + curr_n] = -a_con_val;
      }
    }
  }

  // Cleanup
  delete[] c_arr;
  delete[] delta_d_arr;

  // Cleanup [obstacle portion specific]
  delete[] d_arr;
  delete[] pa_arr;

  /*
  auto gc2_ret = GenerateConstraints2(vehrs, obs_info);
  std::cout << "GC2: " << gc2_ret.a_mat_.rows() << ", " << gc2_ret.a_mat_.cols()
  << std::endl; std::cout << "GC1: " << total_a_con_size << std::endl; if
  (total_a_con_size == gc2_ret.a_mat_.rows()) { for (int i = 0; i <
  total_a_con_size; ++i) { const double new_val = gc2_ret.a_mat_(i, 0); const
  double orig_val = a_con_arr[i]; const double delta = std::abs(new_val -
  orig_val); if (delta > 1.0e-9) { std::cout << "[DBBG] A " << orig_val << " !=
  " << new_val << " [Delta: " << delta << "]" << std::endl;
    }
  }
  } else {
    std::cout << "[DBBG] NEW A != OLD A Rows" << std::endl;
  }

  if (ret.num_b_elts_ == gc2_ret.b_mat_.rows()) {
  for (int i = 0; i < total_a_con_size; ++i) {
    const double new_val = gc2_ret.b_mat_(i, 0);
    const double orig_val = ret.b_con_arr_[i];
    const double delta = std::abs(new_val - orig_val);
    if (delta > 1.0e-9) {
      std::cout << "[DBBG] B " << orig_val << " != " << new_val << " [Delta: "
  << delta << "]" << std::endl;
    }
  }
  } else {
    std::cout << "[DBBG] NEW B != OLD B Rows" << std::endl;

  }
  */

  return ret;
}

void PrintConstraints(const Constraints& constraints, std::ostream& o_stream) {
  const auto& a_con_arr = constraints.a_con_arr_;
  const auto& pb_arr = constraints.b_con_arr_;
  o_stream << "bvals_in = [";
  IndexT endpoint = constraints.num_b_elts_;
  for (IndexT i = 0; i < endpoint; ++i) {
    // o_stream << "pb_arr[" << i << "]: " << pb_arr[i] << std::endl;
    o_stream << pb_arr[i];
    if (i + 1 != endpoint) {
      o_stream << ", ";
    }
  }
  o_stream << "];\n";

  o_stream << "avals_in = [";
  endpoint = constraints.num_a_elts_;
  for (IndexT i = 0; i < endpoint; ++i) {
    o_stream << a_con_arr[i];
    if (i + 1 != endpoint) {
      o_stream << ", ";
    }
  }
  o_stream << "];\n";
}

/*
void WriteObsTestInfo(const FlZonoObsSet& obs_info, std::string fname) {
  std::ofstream f_out;
  f_out.open(fname, std::ofstream::out);

  f_out << "function [num_obs_in, obs_arr_in] = load_obs_info_latest()\n";
  f_out << "num_obs_in = " << obs_info.num_obs_ << ";\n";
  f_out << "obs_arr_in = {};\n";
  for (std::size_t i = 0; i < obs_info.num_obs_; ++i) {
    f_out << "obs_arr_in{end+1} = [" << obs_info.obs_centers_[2 * i + 0] << ", "
          << obs_info.obs_arr_[4 * i + 0] << ", "
          << obs_info.obs_arr_[4 * i + 2] << ";\n"
          << obs_info.obs_centers_[2 * i + 1] << ", "
          << obs_info.obs_arr_[4 * i + 1] << ", "
          << obs_info.obs_arr_[4 * i + 3] << "];\n";
  }
*/
/*
 * TODO CLEANUP
f_out << "obs_arr_in = [";
for (int i = 0; i < obs_info.obs_arr_.size(); ++i) {
  f_out << obs_info.obs_arr_.at(i);
  if (i +1 != obs_info.obs_arr_.size()) {
    f_out << ",";
  }
}
f_out << "];\n";
f_out << "obs_centers_in = [";
for (int i = 0; i < obs_info.obs_centers_.size(); ++i) {
  f_out << obs_info.obs_centers_.at(i);
  if (i +1 != obs_info.obs_centers_.size()) {
    f_out << ",";
  }
}
f_out << "];";
*/
/*
  f_out << "\nend";
  f_out.close();
}
*/

void WriteConstraints(const Constraints& constraints) {
  std::ofstream f_out;
  f_out.open("/home/TODO/constraint_vals.m", std::ofstream::out);
  PrintConstraints(constraints, f_out);
  f_out.close();
}

void PrintConstraints(const Constraints& constraints) {
  PrintConstraints(constraints, std::cout);
}

} // namespace roahm
