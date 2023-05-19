#include "rover_state.hpp"

namespace roahm {

double RoverState::GetHeading() const { return ToAngle(h_); }

void RoverState::SetHeading(double heading) { h_ = ToAngle(heading); }

std::string RoverState::ToStringXYHUVR() const {
  return ::roahm::VariadicListToStr(x_, y_, GetHeading(), u_, v_, r_);
}

PointXYH RoverState::GetXYH() const { return {x_, y_, GetHeading()}; }
} // namespace roahm
