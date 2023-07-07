#include "manu_type.hpp"

#include <cassert> // for assert
#include <ostream> // for operator<<, basic_ostream

#include "ros/console.h" // for LogLocation, ROS_WARN, ROS_WARN_STREAM

namespace roahm {

ManuType FromString(const std::string& manu_type_str) {
  if (manu_type_str == "DIR") {
    return ManuType::kDirChange;
  } else if (manu_type_str == "LAN") {
    return ManuType::kLanChange;
  } else if (manu_type_str == "SPD") {
    return ManuType::kSpdChange;
  } else if (manu_type_str == "NON") {
    return ManuType::kNone;
  } else {
    ROS_WARN_STREAM("ManuType " << manu_type_str << " does not exist");
    return ManuType::kSpdChange;
  }
}

std::string ToString(const ManuType& manu_type) {
  switch (manu_type) {
  case ManuType::kDirChange:
    return "DIR";
  case ManuType::kLanChange:
    return "LAN";
  case ManuType::kSpdChange:
    return "SPD";
  case ManuType::kNone:
    return "NON";
  }
  ROS_WARN("Manu Type Error in ToString");
  assert(false);
  return "UKN";
}

double GetGenParamAu(const ManuType& manu_type, double u0, double k_param) {
  // Au = K for speed change, Au = u0 otherwise
  if (manu_type == ManuType::kSpdChange) {
    return k_param;
  }
  return u0;
}

double GetGenParamAy(const ManuType& manu_type, double u0, double k_param) {
  if (manu_type == ManuType::kSpdChange) {
    return 0.0;
  }
  return k_param;
}

double GetPreBrakeTime(ManuType manu_type) {
  // TODO from consts
  if (IsSpdDir(manu_type)) {
    return 1.5;
  } else if (IsLan(manu_type)) {
    return 3.0;
  } else if (IsNone(manu_type)) {
    return 0.0;
  }
  // TODO warn
  return 0.0;
}

std::uint8_t ManuToUint8(ManuType manu_type) {
  if (IsSpd(manu_type)) {
    return '0';
  } else if (IsDir(manu_type)) {
    return '1';
  } else if (IsLan(manu_type)) {
    return '2';
  } else if (IsNone(manu_type)) {
    return 'N';
  }
  ROS_WARN("Unknown manu type in conversion to uint8");
  assert(false);
  return 'U';
}

ManuType ManuFromUint8(std::uint8_t val) {
  if (val == '0') {
    return ManuType::kSpdChange;
  } else if (val == '1') {
    return ManuType::kDirChange;
  } else if (val == '2') {
    return ManuType::kLanChange;
  } else if (val == 'N') {
    return ManuType::kNone;
  }
  ROS_WARN_STREAM("Manu Type Character '" << val << "' not recognized");
  assert(false);
  return ManuType::kNone;
}

} // namespace roahm
