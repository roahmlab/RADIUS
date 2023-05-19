#include "cuda_info.hpp"

#include "frs_io.hpp"
#include "frs_loader.hpp"

namespace roahm {
namespace {
template <typename T> static void PrintVar(const T& var, std::string name) {
  std::cout << name << ": " << var << std::endl;
}

template <typename T>
static void PrintVec(const std::vector<T>& vec, std::string name) {
  std::cout << name << " [N: " << vec.size() << "]: ";
  if (not vec.empty()) {
    for (int i = 0; (i + 1) < vec.size(); ++i) {
      std::cout << vec.at(i) << ", ";
    }
    std::cout << vec.back();
  }
  std::cout << std::endl;
}
} // namespace

void CudaInfo::Print() const {
  PrintVar(grid_size_, "grid_size_");
  PrintVar(num_zono_, "num_zono_");
  PrintVar(manu_type_, "manu_type_");
  PrintVec(grid_x0_, "grid_x0_");
  PrintVec(grid_y0_, "grid_y0_");
  PrintVec(grid_dx_, "grid_dx_");
  PrintVec(grid_dy_, "grid_dy_");
  PrintVec(g_u0_x_, "g_u0_x_");
  PrintVec(g_u0_y_, "g_u0_y_");
  PrintVec(g_v0_x_, "g_v0_x_");
  PrintVec(g_v0_y_, "g_v0_y_");
  PrintVec(g_r0_x_, "g_r0_x_");
  PrintVec(g_r0_y_, "g_r0_y_");
  PrintVec(block_inzono_list_, "block_inzono_list_");
  PrintVec(rot_angle_, "rot_angle_");
  PrintVec(cg_p_, "cg_p_");
  PrintVec(g_p_x_, "g_p_x_");
  PrintVec(cg_u0v0r0_, "cg_u0v0r0_");
  PrintVec(t_range_, "t_range_");
}

bool CudaInfo::operator==(const CudaInfo& oth) const {
  return grid_size_ == oth.grid_size_ and num_zono_ == oth.num_zono_ and
         manu_type_ == oth.manu_type_ and VecsEq(grid_x0_, oth.grid_x0_) and
         VecsEq(grid_y0_, oth.grid_y0_) and VecsEq(grid_dx_, oth.grid_dx_) and
         VecsEq(grid_dy_, oth.grid_dy_) and VecsEq(g_u0_x_, oth.g_u0_x_) and
         VecsEq(g_u0_y_, oth.g_u0_y_) and VecsEq(g_v0_x_, oth.g_v0_x_) and
         VecsEq(g_v0_y_, oth.g_v0_y_) and VecsEq(g_r0_x_, oth.g_r0_x_) and
         VecsEq(g_r0_y_, oth.g_r0_y_) and
         VecsEq(block_inzono_list_, oth.block_inzono_list_) and
         VecsEq(rot_angle_, oth.rot_angle_) and VecsEq(cg_p_, oth.cg_p_) and
         VecsEq(g_p_x_, oth.g_p_x_) and VecsEq(cg_u0v0r0_, oth.cg_u0v0r0_) and
         VecsEq(t_range_, oth.t_range_);
}

void CudaInfo::WriteToBinFileDirect(std::ofstream& file) const {
  WriteToBinFile(file, grid_size_);
  WriteToBinFile(file, num_zono_);
  WriteToBinFile(file, manu_type_);
  WriteToBinFile(file, grid_x0_);
  WriteToBinFile(file, grid_y0_);
  WriteToBinFile(file, grid_dx_);
  WriteToBinFile(file, grid_dy_);
  WriteToBinFile(file, g_u0_x_);
  WriteToBinFile(file, g_u0_y_);
  WriteToBinFile(file, g_v0_x_);
  WriteToBinFile(file, g_v0_y_);
  WriteToBinFile(file, g_r0_x_);
  WriteToBinFile(file, g_r0_y_);
  WriteToBinFile(file, block_inzono_list_);
  WriteToBinFile(file, rot_angle_);
  WriteToBinFile(file, cg_p_);
  WriteToBinFile(file, g_p_x_);
  WriteToBinFile(file, cg_u0v0r0_);
  WriteToBinFile(file, t_range_);
}

CudaInfo CudaInfo::ReadFromBinFileDirect(std::ifstream& file) {
  CudaInfo ret;
  ret.grid_size_ = ReadFromBinFile<int>(file);
  ret.num_zono_ = ReadFromBinFile<int>(file);
  ret.manu_type_ = ReadFromBinFile<int>(file);
  ret.grid_x0_ = ReadFromBinFile<std::vector<double>>(file);
  ret.grid_y0_ = ReadFromBinFile<std::vector<double>>(file);
  ret.grid_dx_ = ReadFromBinFile<std::vector<double>>(file);
  ret.grid_dy_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_u0_x_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_u0_y_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_v0_x_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_v0_y_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_r0_x_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_r0_y_ = ReadFromBinFile<std::vector<double>>(file);
  ret.block_inzono_list_ = ReadFromBinFile<std::vector<double>>(file);
  ret.rot_angle_ = ReadFromBinFile<std::vector<double>>(file);
  ret.cg_p_ = ReadFromBinFile<std::vector<double>>(file);
  ret.g_p_x_ = ReadFromBinFile<std::vector<double>>(file);
  ret.cg_u0v0r0_ = ReadFromBinFile<std::vector<double>>(file);
  ret.t_range_ = ReadFromBinFile<std::vector<double>>(file);
  return ret;
}

} // namespace roahm