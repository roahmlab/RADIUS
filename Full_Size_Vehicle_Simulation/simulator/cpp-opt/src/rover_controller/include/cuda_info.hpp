#ifndef ROAHM_CUDA_INFO_HPP_
#define ROAHM_CUDA_INFO_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace roahm {

struct CudaInfo {
  int grid_size_;
  int num_zono_;
  int manu_type_;
  std::vector<double> grid_x0_;
  std::vector<double> grid_y0_;
  std::vector<double> grid_dx_;
  std::vector<double> grid_dy_;
  std::vector<double> g_u0_x_;
  std::vector<double> g_u0_y_;
  std::vector<double> g_v0_x_;
  std::vector<double> g_v0_y_;
  std::vector<double> g_r0_x_;
  std::vector<double> g_r0_y_;
  /// Whether each block is in the zonotope
  /// -1: not in zonotope
  /// other: index, in zonotope
  /// Starts in lower left, goes right, then up (row major but flipped)
  std::vector<double> block_inzono_list_;
  std::vector<double> rot_angle_;
  std::vector<double> cg_p_;
  std::vector<double> g_p_x_;
  std::vector<double> cg_u0v0r0_;
  std::vector<double> t_range_;

  bool operator==(const CudaInfo& oth) const;

  void WriteToBinFileDirect(std::ofstream& file) const;
  static CudaInfo ReadFromBinFileDirect(std::ifstream& file);

  void Print() const;
};

} // namespace roahm
#endif