#ifndef ROAHM_PRE_SLICE_HPP_
#define ROAHM_PRE_SLICE_HPP_

#include <array>
#include <memory>

#include "cuda_info.hpp"

namespace roahm::pre_slice {
struct PreSliceOutputs {
  std::shared_ptr<double> x0_;
  std::shared_ptr<double> y0_;
  std::shared_ptr<double> H1_;
  std::shared_ptr<double> H2_;
  std::shared_ptr<double> H4_;
  std::shared_ptr<float> gpu_time_;
};

PreSliceOutputs
PreSliceImpl(const int grid_size, const int num_zono, const double* host_x0,
             const double* host_y0, const double* host_dx,
             const double* host_dy, const double* host_u0v0r0_slice_beta,
             const double* host_g_u0_x, const double* host_g_u0_y,
             const double* host_g_v0_x, const double* host_g_v0_y,
             const double* host_g_r0_x, const double* host_g_r0_y,
             const double* host_block_inzono_list, const double* host_rot_angle,
             const double* host_mu_sigma, const double* host_g_p_x);

PreSliceOutputs PreSlice(const ::roahm::CudaInfo& cuda_info,
                         const std::array<double, 3>& u0v0r0_slice_beta_in,
                         const double* const mu_sigma);
} // namespace roahm::pre_slice

#endif // ROAHM_PRE_SLICE_HPP_