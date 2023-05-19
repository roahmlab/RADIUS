#include "individual_mu_sigma.hpp"

#include <cmath>

#include "point_xyh.hpp"

namespace roahm {
IndividualMuSigma IndividualMuSigma::RelativeTo(const ::roahm::PointXYH& frame,
                                                const bool mirror) const {
  const double x0 = frame.x_;
  const double y0 = frame.y_;
  const double h0 = frame.h_;
  const double cos_h = std::cos(h0);
  const double sin_h = std::sin(h0);
  const double dx = mu_x_ - x0;
  const double dy = mu_y_ - y0;

  // R(-h) * [mu_x - x0; mu_y - y0]
  const double out_mu_x = cos_h * dx + sin_h * dy;
  const double out_mu_y = -sin_h * dx + cos_h * dy;
  const double A = sigma_1_;
  const double B = sigma_2_;
  const double D = sigma_4_;
  const double sigma_1 =
      cos_h * (A * cos_h + B * sin_h) + sin_h * (B * cos_h + D * sin_h);
  const double sigma_2 =
      sin_h * (D * cos_h - B * sin_h) + cos_h * (B * cos_h - A * sin_h);
  const double sigma_4 =
      -sin_h * (B * cos_h - A * sin_h) + cos_h * (D * cos_h - B * sin_h);
  const double mirror_mult = mirror ? -1 : 1;
  return {out_mu_x, out_mu_y * mirror_mult, sigma_1, sigma_2 * mirror_mult,
          sigma_4};
}
} // namespace roahm