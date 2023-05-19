#ifndef ROAHM_GENCON_HPP_
#define ROAHM_GENCON_HPP_

#include <eigen3/Eigen/Dense>
#include <iosfwd> // for ostream
#include <memory> // for unique_ptr, shared_ptr
#include <string> // for string
#include <vector> // for vector

#include "common.hpp" // for IndexT
#include "dyn_obs.hpp"
#include "fl_zono_obs_set.hpp" // for FlZonoObsSet
#include "point_xyh.hpp"       // for PointXYH
#include "vehrs.hpp"           // for Vehrs

/// @file gencon.hpp Contains constraint generation code and related utilities

namespace roahm {

/// TODO
struct ZonoInfo {
  /// The number of zonotopes, \f$ n \f$
  IndexT num_zonos_;
  /// TODO
  IndexT cumsum_num_out_zono_gens_;
  /// TODO
  std::unique_ptr<IndexT[]> num_out_zono_gens_arr_;
  /// TODO
  std::unique_ptr<IndexT[]> cum_size_arr_;
  /// TODO
  std::unique_ptr<double[]> zono_arr_;

  /// Default constructor
  ZonoInfo()
      : num_zonos_{0}, cumsum_num_out_zono_gens_{0}, num_out_zono_gens_arr_{},
        cum_size_arr_{}, zono_arr_{} {}
};

/// TODO
struct Constraints {
  /// Number of elements in the a matrix
  IndexT num_a_elts_;
  /// Number of elements in the b vector
  IndexT num_b_elts_;
  /// TODO
  std::shared_ptr<double[]> a_con_arr_;
  /// TODO
  std::shared_ptr<double[]> b_con_arr_;
  /// TODO
  std::vector<IndexT> zono_startpoints_;
  /// TODO
  std::vector<IndexT> zono_obs_sizes_;
  /// Default constructor
  Constraints()
      : num_a_elts_{0}, num_b_elts_{0}, a_con_arr_{}, b_con_arr_{},
        zono_startpoints_{}, zono_obs_sizes_{} {}
};

/// Generates constraints given a reachable set, predicted \f$ (u, v, r) \f$
/// which the reachable set can be sliced at, and obstacles
/// \param vehrs the reachable set, sliced with initial conditions
/// \param obs_info the obstacles in the local frame relative to the predicted
/// pose
/// \return constraints to be used in evaluating the optimization problem
Constraints GenerateConstraints(const Vehrs& vehrs,
                                const FlZonoObsSet& obs_info);

class Gencon2Out {
public:
  // WRT XY
  Eigen::Matrix<double, Eigen::Dynamic, 1> a_mat_;
  Eigen::VectorXd b_mat_;
  Eigen::VectorXd start_idx_mat_;
  Eigen::VectorXd size_mat_;
};

// Just a (currently very slow), but correct Eigen implementation for REFINE
// constraints. May come in handy as a ground truth or if we do REFINE/RISK
// constraint selection on a (zono, obs) pair basis instead of max vs.
// the full FRS.
Gencon2Out GenerateConstraints2(const Vehrs& reach_set,
                                const FlZonoObsSet& local_obs);

/// Writes constraints in a form that is MATLAB-parseable to a std::ostream
/// \param constraints the constraints to write
/// \param o_stream the std::ostream to write to
void PrintConstraints(const Constraints& constraints, std::ostream& o_stream);

/// Writes obstacle information to a file that is MATALB-parseable
/// \param obs_info the obstacle information to write
/// \param fname the file name to write the output to
void WriteObsTestInfo(const FlZonoObsSet& obs_info, std::string fname);

/// Prints constraints to std::cout in a MATLAB readable form
/// \param constraints the constraints to write to the file
void WriteConstraints(const Constraints& constraints);

/// Writes constraints to file that is readable by MATLAB
/// \param constraints the constraints to write to the file
void PrintConstraints(const Constraints& constraints);
} // namespace roahm

#endif // ROAHM_GENCON_HPP_
