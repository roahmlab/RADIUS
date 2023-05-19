#include "vehrs_plot.hpp"

#include <gtest/gtest.h>

#include <chrono>

#include "frs_loader.hpp"
#include "frs_total.hpp"
#include "plot_interface.hpp"
#include "timing_util.hpp"

TEST(VehrsPlotTest, LoadFromFileTest) {
  const ::roahm::FrsTotal frs{
      ::roahm::LoadFrsBinary("/data/frs_bin_processed.txt")};
  roahm::plot_utils::OpenFigure(19);
  const auto& frs_unsliced = frs.megas_.at(5).dir_.at(0).at(0);
  const auto frs_sliced =
      frs_unsliced.SliceAtParam(frs_unsliced.GetUCenterGen().first, 0.0, 0.0,
                                frs_unsliced.GetTrajParamMax());
  const ::roahm::plot_utils::PlotInfo plot_info_unsliced{
      ::roahm::plot_utils::PlotColor{0.01, 0.0, 0.0, 1.0}, 0.5, false,
      ::roahm::plot_utils::PlotLineType::Continuous()};
  roahm::PlotVehrs(frs_unsliced, {0, 0, 0}, false, plot_info_unsliced);
  roahm::PlotVehrs(frs_unsliced, {0, 0, 0}, true, plot_info_unsliced);
  const ::roahm::plot_utils::PlotInfo plot_info_sliced{
      ::roahm::plot_utils::PlotColor{0.01, 0.0, 1.0, 0.0}, 0.5, false,
      ::roahm::plot_utils::PlotLineType::Continuous()};
  roahm::PlotVehrs(frs_sliced, {0.0, 0.0, 0.0}, false, plot_info_sliced);
  roahm::plot_utils::SetAxisEqual();
  roahm::plot_utils::Show();
}