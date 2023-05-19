#ifndef ROAHM_FL_ZONO_OBS_SET_HPP_
#define ROAHM_FL_ZONO_OBS_SET_HPP_

#include <vector>

#include "dyn_obs.hpp"

namespace roahm {

/// Contains zonotope representations of a set of obstacles, which have two
/// generators per obstacle
struct FlZonoObsSet {
private:
  std::vector<DynObs> dyn_obs_;

public:
  FlZonoObsSet RelativeTo(PointXYH pt, bool mirror) const;

  /// Default constructor
  FlZonoObsSet();

  struct CenterXY {
    double x_;
    double y_;
  };

  struct GenXY {
    double x1_;
    double y1_;
    double x2_;
    double y2_;
  };

  void ClearAll();

  std::size_t GetNumObs() const;

  void PushObs(DynObs obs);

  void PrintObsConv(const int obs_idx,
                    const ::roahm::Interval t_interval) const;
  CenterXY ComputeCenter(const int obs_idx,
                         const ::roahm::Interval t_interval) const;

  CenterXY GetCenter(const int obs_idx,
                     const ::roahm::Interval t_interval) const {
    // PrintObsConv(obs_idx, t_interval);
    return ComputeCenter(obs_idx, t_interval);
  }

  GenXY GetGenerators(int obs_idx, ::roahm::Interval t_interval) const;

  CenterXY GetSingleGenerator(int obs_idx, int gen_idx,
                              ::roahm::Interval t_interval) const;

  // TODO move to impl
  [[nodiscard]] DynObs GetSingleDynObs(int obs_idx) const;
};

} // namespace roahm
#endif // ROAHM_FL_ZONO_OBS_SET_HPP_