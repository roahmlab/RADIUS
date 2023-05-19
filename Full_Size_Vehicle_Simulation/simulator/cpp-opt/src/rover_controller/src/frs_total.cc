#include "frs_total.hpp"

#include "frs_io.hpp"
namespace roahm {
namespace {

/// Finds the minimum element in a container, given a comparison operator
/// \tparam ContainerType the type of container to search over, must provide
/// begin and end methods that provide forward iterators
/// \tparam Op The comparison operator type, must provide a valid comparison to
/// be provided to std::min_element
/// \param container the list of values to search
/// \param op the instance of the comparison operator
/// \return an iterator to the minimum element of the provided container,
/// according to the provided comparison operator
template <typename ContainerType, typename Op>
auto MinElement(const ContainerType& container, const Op& op) {
  return std::min_element(container.begin(), container.end(), op);
}

/// Finds the index of the minimum element in a container, given a
/// comparison operator
/// \tparam ContainerType the type of container to search over, must provide
/// begin and end methods that provide forward iterators
/// \tparam Op The comparison operator type, must provide a valid comparison to
/// be provided to std::min_element
/// \param container the list of values to search
/// \param op the instance of the comparison operator
/// \return the index of the minimum element according to the provided
/// comparison operator \p op
template <typename ContainerType, typename Op>
auto MinElementIdx(const ContainerType& container, const Op& op) {
  return std::distance(container.begin(), MinElement(container, op));
}
} // namespace

const Vehrs& FrsTotal::GetVehrs(const FrsSelectInfo& info) const {
  const auto& curr_mega = megas_.at(info.idxu0_);
  return (info.manu_type_ == ManuType::kSpdChange)
             ? (curr_mega.au_.at(info.idx0_))
             : ((info.manu_type_ == ManuType::kDirChange)
                    ? curr_mega.dir_.at(info.idx0_).at(info.idx1_)
                    : curr_mega.lan_.at(info.idx0_).at(info.idx1_));
}
bool FrsTotal::AlternatingAu() const { return alternating_au_; }
void FrsTotal::WriteToBinFileDirect(std::ofstream& file) const {
  WriteToBinFile(file, u0_intervals_);
  WriteToBinFile(file, megas_);
  WriteToBinFile(file, u0_min_);
  WriteToBinFile(file, u0_max_);
  WriteToBinFile(file, successful_);
  WriteToBinFile(file, alternating_au_);
}

FrsTotal FrsTotal::ReadFromBinFileDirect(std::ifstream& file) {
  FrsTotal frs_total;
  ReadToVarFromBinFile(frs_total.u0_intervals_, file);
  ReadToVarFromBinFile(frs_total.megas_, file);
  ReadToVarFromBinFile(frs_total.u0_min_, file);
  ReadToVarFromBinFile(frs_total.u0_max_, file);
  ReadToVarFromBinFile(frs_total.successful_, file);
  ReadToVarFromBinFile(frs_total.alternating_au_, file);
  return frs_total;
}
bool FrsTotal::NearEqual(const FrsTotal& rhs, double tol) const {
  const bool u0_intervals_near_eq = VecsNearEq(u0_intervals_, rhs.u0_intervals_,
                                               tol, "FrsTotal::u0_intervals");
  const bool megas_near_eq =
      VecsNearEq(megas_, rhs.megas_, tol, "FrsTotal::megas");
  const bool u0_min_near_eq = NearEq(u0_min_, rhs.u0_min_, tol);
  const bool u0_max_near_eq = NearEq(u0_max_, rhs.u0_max_, tol);
  const bool successful_eq = successful_ == rhs.successful_;
  const bool alternating_au_eq = alternating_au_ == rhs.alternating_au_;
  return u0_intervals_near_eq and megas_near_eq and u0_min_near_eq and
         u0_max_near_eq and successful_eq and alternating_au_eq;
}
bool FrsTotal::operator==(const FrsTotal& other) const {
  return VecsEq(u0_intervals_, other.u0_intervals_, "Total::u0_intervals") &&
         VecsEq(megas_, other.megas_, "Total::megas") &&
         u0_min_ == other.u0_min_ && u0_max_ == other.u0_max_ &&
         successful_ == other.successful_ &&
         alternating_au_ == other.alternating_au_;
}
double FrsTotal::MinU0() const { return u0_min_; }
double FrsTotal::MaxU0() const { return u0_max_; }

long FrsTotal::SelectU0IdxNoClamp(const double u0) const {
  const auto comp_dist_to_u0 = [u0](const Interval& a, const Interval& b) {
    return a.DistanceTo(u0) < b.DistanceTo(u0);
  };
  return static_cast<long>(MinElementIdx(u0_intervals_, comp_dist_to_u0));
}

long FrsTotal::SelectU0Idx(double& u0, bool can_clamp_u) const {
  if (can_clamp_u) {
    u0 = ClampWithWarn(u0, MinU0(), MaxU0(), "U0");
  }
  return SelectU0IdxNoClamp(u0);
}

} // namespace roahm