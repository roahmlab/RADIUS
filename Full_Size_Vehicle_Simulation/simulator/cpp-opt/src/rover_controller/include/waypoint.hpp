#ifndef ROAHM_WAYPOINT_HPP_
#define ROAHM_WAYPOINT_HPP_

#include "manu_type.hpp"
#include "point_xyh.hpp"
namespace roahm {

constexpr bool Implies(const bool a, const bool b) {
  // returns true iff a => b
  return (not a) or b;
}

template <bool EgoFrame, bool MirrorTakenIntoAccount>
  requires(Implies(MirrorTakenIntoAccount, EgoFrame))
struct Waypoint : public PointXYH {
  [[nodiscard]] auto RelativeToEgoFrameNoMirror(const PointXYH& ego_frame) const
    requires((not EgoFrame) and (not MirrorTakenIntoAccount))
  {
    return Waypoint<true, false>{ToLocalFrame(ego_frame)};
  }

  [[nodiscard]] auto TakeMirrorIntoAccount(bool mirror) const
    requires(not MirrorTakenIntoAccount)
  {
    if (mirror) {
      return Waypoint<EgoFrame, true>{PointXYH::Mirror()};
    } else {
      return Waypoint<EgoFrame, true>{*this};
    }
  }

  [[nodiscard]] inline Waypoint<true, true>
  HeuristicAdjust(const ::roahm::ManuType& manu_type, const double u0) const
    requires(EgoFrame and MirrorTakenIntoAccount)
  {
    const double h{h_};
    double x{x_};
    double y{y_};
    if (IsSpd(manu_type)) {
      x = 1.5 * (std::min(x / 3.0, u0 + 4.0) + u0);
    } else if (IsLan(manu_type)) {
      x = (u0)*6.0;
      y *= 1.0;
    } else if (IsDir(manu_type)) {
      x = (u0 + 1.0) * 3.0;
      y *= 0.7;
    }
    return Waypoint<true, true>{PointXYH{x, y, h}};
  }
};

using WaypointLocalNoMirrorTakenIntoAccount = Waypoint<true, false>;
using WaypointLocalMirrorTakenIntoAccount = Waypoint<true, true>;
using WaypointGlobalNoMirror = Waypoint<false, false>;

} // namespace roahm

#endif // ROAHM_WAYPOINT_HPP_