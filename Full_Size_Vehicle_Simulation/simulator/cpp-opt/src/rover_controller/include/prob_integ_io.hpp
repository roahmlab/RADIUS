#ifndef ROAHM_PROB_INTEG_IO_H_
#define ROAHM_PROB_INTEG_IO_H_

#include <memory>

namespace roahm {
struct ProbIntegrationInputs {
  // TODO move documentation here
  int num_zonos_;
  const double* x0_;
  const double* y0_;
  const double* dx_;
  const double* dy_;
  const double* block_inzono_list_;
  const double* rot_angle_;
  const double* mu_sigma_;
  const double* cg_p_;
  const double* g_p_x_;
  int grid_size_;
  double p_;
  double* H1_;
  double* H2_;
  double* H4_;
};

struct ProbIntegrationOutputs {
  // TODO move documentation here
  double constraint_val_;
  double d_constraint_val_;
  double d2_constraint_val_;
  float computation_time_;
};
} // namespace roahm

#endif // ROAHM_PROB_INTEG_IO_H_