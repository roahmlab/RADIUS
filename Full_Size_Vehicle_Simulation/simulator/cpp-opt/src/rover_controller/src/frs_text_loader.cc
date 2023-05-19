#include "frs_text_loader.hpp"
namespace roahm {

namespace {

/// Prints an error message when two values are not equal
/// \tparam T the values' type
/// \param a_name the name of the first value
/// \param b_name the name of the second value
/// \param a the first value
/// \param b the second value
template <typename T>
void CheckEqMsg(const std::string& a_name, const std::string& b_name,
                const T& a, const T& b) {
  std::cerr << a_name << " != " << b_name << " (" << a << " != " << b << ")\n";
}

/// Loads an individual zonotope as part of loading an FRS
/// \param f_in the input stream to read from
/// \param vehrs out param FRS to add the zonotope to
/// \return true iff successful
bool LoadZono(std::istream& f_in, Vehrs& vehrs) {
  std::string in_str;
  IndexT num_dims = 0;
  IndexT num_gens = 0;
  f_in >> in_str >> num_dims >> num_gens;
  if (num_dims != 3) { // XYH
    std::cerr << "num_dims != 3\n";
    return false;
  }
  vehrs.zono_sizes_.push_back(num_gens);
  for (IndexT gen_idx = 0; gen_idx < num_gens; ++gen_idx) {
    double x, y, h;
    f_in >> x >> y >> h;
    vehrs.zono_xy_.push_back(x);
    vehrs.zono_xy_.push_back(y);
    vehrs.zono_h_.push_back(h);
  }

  double c_x, c_y, c_h;
  f_in >> in_str >> c_x >> c_y >> c_h;
  if (in_str != "CENTER") {
    std::cerr << "in_str != \"CENTER\" (\"" << in_str << "\" != CENTER)\n";
    return false;
  }
  vehrs.xy_centers_.push_back(c_x);
  vehrs.xy_centers_.push_back(c_y);
  vehrs.h_centers_.push_back(c_h);
  IndexT num_slc_vals;
  f_in >> in_str >> num_slc_vals;
  if (in_str != "SLICE_VALS:") {
    std::cerr << "in_str != \"SLICE_VALS:\"\n";
    return false;
  }
  vehrs.slc_infos_.emplace_back();
  for (IndexT i = 0; i < num_slc_vals; ++i) {
    int dim;
    double center_val, slc_val, x, y, h;
    f_in >> in_str >> dim >> center_val >> slc_val >> x >> y >> h;
    if (in_str != "SLC") {
      std::cerr << "in_str != \"SLC\"\n";
      return false;
    }
    vehrs.slc_infos_.back().slc_vals_.emplace_back(dim, center_val, slc_val, x,
                                                   y, h);
  }
  f_in >> in_str;
  if (in_str != "TIME_RNG:") {
    std::cerr << "in_str != \"TIME_RNG:\"\n";
    return false;
  }
  double t_lb{};
  double t_ub{};
  if (not(f_in >> t_lb >> t_ub)) {
    std::cerr << "Could not read [t_lb, t_ub]\n";
    return false;
  }
  vehrs.zono_time_intervals_.emplace_back(t_lb, t_ub);
  return true;
}

/// Loads a single FRS (forward reachable set) from an existing input stream
/// \param f_in the input stream to read from
/// \param vehrs out param to write the FRS to
/// \return true iff successful
bool LoadVehrs(std::istream& f_in, Vehrs& vehrs) {
  std::string in_str;
  IndexT num_zonos;
  f_in >> in_str; // ({MANU_TYPE}_{IDX})
  f_in >> in_str; // vehRS_save
  if (in_str != "vehRS_save") {
    std::cerr << "in_str != \"vehRS_save\"\n";
    return false;
  }
  int t_eval_idx = -1;
  f_in >> num_zonos >> t_eval_idx;
  if (t_eval_idx < 0) {
    std::cerr << "t_eval_idx <= 0\n";
    return false;
  }

  f_in >> in_str;
  if (in_str != "K_RNG") {
    std::cout << "in_str != \"K_RNG\"\n";
    return false;
  }
  double k_rng = 0.0;
  f_in >> k_rng;
  if (k_rng < 1.0e-2) {
    std::cout << "k_rng < 1.0e-2\n";
    return false;
  }
  std::cout << "k_rng: " << k_rng << std::endl;
  vehrs.k_rng_ = k_rng;

  vehrs.t_eval_idx_ = t_eval_idx;
  for (IndexT i = 0; i < num_zonos; ++i) {
    if (!LoadZono(f_in, vehrs)) {
      std::cerr << "LoadZono failed\n";
      return false;
    }
  }
  return true;
}

/// Load one set of FRSes from an existing input stream
/// \param f_in the stream to read from
/// \param mega out param to write the data to
/// \return true iff successful
bool LoadMega(std::istream& f_in, FrsMega& mega) {
  std::string in_str;
  f_in >> in_str;
  if (in_str != "M_MEGA") {
    std::cerr << "in_str != \"M_MEGA\"\n";
    return false;
  }
  f_in >> in_str; // mega idx
  f_in >> in_str; // (MEGA_{N})

  // Load AU
  std::string manu_type;
  IndexT au_num = 0;
  f_in >> manu_type >> au_num;
  if (manu_type != "AU") {
    std::cerr << "manu_type != \"AU\"\n";
    return false;
  }
  mega.au_.resize(au_num);
  for (IndexT i = 0; i < au_num; ++i) {
    if (!LoadVehrs(f_in, mega.au_.at(i))) {
      std::cerr << "AU LoadVehrs failed\n";
      return false;
    }
  }

  // Load AUTB
  IndexT autb_rows, autb_cols;
  f_in >> in_str >> in_str >> autb_rows >> autb_cols;
  if (in_str != "AUTB") {
    std::cerr << "in_str != \"AUTB\"\n";
    return false;
  }
  // TODO if (autb_rows != 14) {
  // TODO   std::cerr << "autb_rows != 14\n";
  // TODO   return false;
  // TODO }
  if (autb_cols != au_num) {
    std::cerr << "autb_cols [" << autb_cols << "] != au_num [" << au_num
              << "]\n";
    return false;
  }
  // NOTE: this temporary is just handle the old format, we no longer store autb
  // anywhere
  std::vector<std::vector<double>> autb_tmp;
  autb_tmp.resize(autb_rows);
  for (IndexT autb_r = 0; autb_r < autb_rows; ++autb_r) {
    for (IndexT autb_c = 0; autb_c < autb_cols; ++autb_c) {
      double tmp_val;
      if (!(f_in >> tmp_val)) {
        std::cerr << "!(f_in >> tmp_val)\n";
        return false;
      }
      autb_tmp.at(autb_r).push_back(tmp_val);
    }
  }
  if (autb_tmp.size() != autb_rows) {
    CheckEqMsg("autb size", "autb rows", autb_tmp.size(), autb_rows);
    // std::cerr << "autb_tmp.size() != autb_rows * autb_cols\n";
    return false;
  }
  for (auto& autb_row : autb_tmp) {
    if (autb_row.size() != autb_cols) {
      CheckEqMsg("autb row size", "autb cols", autb_row.size(), autb_cols);
      return false;
    }
  }

  // Load DIR
  f_in >> in_str >> in_str;
  if (in_str != "DIR") {
    std::cerr << "in_str != \"DIR\"\n";
    return false;
  }
  IndexT dir_rows, dir_cols;
  f_in >> dir_rows >> dir_cols;
  if (not(dir_rows > 0)) {
    std::cerr << "check (dir_rows > 0) failed\n";
    return false;
  }

  mega.dir_.resize(dir_rows);
  for (auto& dir_column : mega.dir_) {
    dir_column.resize(dir_cols);
  }
  for (IndexT dir_r = 0; dir_r < dir_rows; ++dir_r) {
    for (IndexT dir_c = 0; dir_c < dir_cols; ++dir_c) {
      if (!LoadVehrs(f_in, mega.dir_.at(dir_r).at(dir_c))) {
        std::cerr << "DIR LoadVehrs failed\n";
        return false;
      }
    }
  }

  // Load DIRTB
  IndexT dirtb_dim0, dirtb_dim1, dirtb_dim2;
  f_in >> in_str >> in_str >> dirtb_dim0 >> dirtb_dim1 >> dirtb_dim2;

  // Verify things
  if (in_str != "DIRTB") {
    std::cerr << "in_str != \"DIRTB\"\n";
    return false;
  }
  if (dirtb_dim2 > 0) {
    // if (dirtb_dim0 != 14) {
    //   std::cerr << "dirtb_dim1 != 14\n";
    //   return false;
    // }
    if (dirtb_dim1 != dir_rows) {
      std::cerr << "dirtb_dim1 != dir_rows\n";
      return false;
    }
  }
  if (dirtb_dim2 != dir_cols) {
    std::cerr << "dirtb_dim2 != dir_cols\n";
    return false;
  }

  // NOTE: this temporary is just handle the old format, we no longer store
  // dirtb anywhere
  std::vector<std::vector<std::vector<double>>> dirtb_tmp;
  dirtb_tmp.resize(dirtb_dim0);
  for (IndexT dirtb_i = 0; dirtb_i < dirtb_dim0; ++dirtb_i) {
    dirtb_tmp.at(dirtb_i).resize(dirtb_dim1);
    for (IndexT dirtb_j = 0; dirtb_j < dirtb_dim1; ++dirtb_j) {
      for (IndexT dirtb_k = 0; dirtb_k < dirtb_dim2; ++dirtb_k) {
        double tmp_val;
        f_in >> tmp_val;
        dirtb_tmp.at(dirtb_i).at(dirtb_j).push_back(tmp_val);
      }
    }
  }
  if (dirtb_tmp.size() != dirtb_dim0) {
    std::cerr << "dirtb_tmp.size() != dirtb_dim0\n";
    return false;
  }
  for (const auto& dirtb_col : dirtb_tmp) {
    if (dirtb_col.size() != dirtb_dim1) {
      std::cerr << "dirtb_tmp.at(i).size() != dirtb_dim1\n";
      return false;
    }
    for (const auto& dirtb_t0 : dirtb_col) {
      if (dirtb_t0.size() != dirtb_dim2) {
        std::cerr << "dirtb_tmp.at(i).at(j).size() != dirtb_dim2\n";
        return false;
      }
    }
  }

  // Load LAN
  f_in >> in_str >> in_str;
  if (in_str != "LAN") {
    std::cerr << "in_str != \"LAN\"\n";
    return false;
  }
  IndexT lan_rows, lan_cols;
  f_in >> lan_rows >> lan_cols;

  // if (not(lan_rows > 0)) {
  //   std::cerr << "check (lan_rows > 0) failed\n";
  //   return false;
  // }

  mega.lan_.resize(lan_rows);
  for (auto& lan_column : mega.lan_) {
    lan_column.resize(lan_cols);
  }
  for (IndexT lan_r = 0; lan_r < lan_rows; ++lan_r) {
    for (IndexT lan_c = 0; lan_c < lan_cols; ++lan_c) {
      if (!LoadVehrs(f_in, mega.lan_.at(lan_r).at(lan_c))) {
        std::cerr << "LAN LoadVehrs failed\n";
        return false;
      }
    }
  }

  // Load LANTB
  IndexT lantb_dim0, lantb_dim1, lantb_dim2;
  f_in >> in_str >> in_str >> lantb_dim0 >> lantb_dim1 >> lantb_dim2;
  if (in_str != "LANTB") {
    std::cerr << "in_str != \"LANTB\"\n";
    return false;
  }
  if (lantb_dim2 > 0) {
    if (lantb_dim0 != 14) {
      std::cerr << "lantb_dim0 != 14\n";
      return false;
    }
    if (lantb_dim1 != lan_rows) {
      std::cerr << "lantb_dim1 != lan_rows\n";
      return false;
    }
  }
  if (lantb_dim2 != lan_cols) {
    std::cerr << "lantb_dim2 != lan_cols\n";
    return false;
  }
  // NOTE: this temporary is just handle the old format, we no longer store
  // lantb anywhere
  std::vector<std::vector<std::vector<double>>> lantb_tmp;
  lantb_tmp.resize(lantb_dim0);
  for (IndexT lantb_i = 0; lantb_i < lantb_dim0; ++lantb_i) {
    lantb_tmp.at(lantb_i).resize(lantb_dim1);
    for (IndexT lantb_j = 0; lantb_j < lantb_dim1; ++lantb_j) {
      for (IndexT lantb_k = 0; lantb_k < lantb_dim2; ++lantb_k) {
        double tmp_val;
        f_in >> tmp_val;
        lantb_tmp.at(lantb_i).at(lantb_j).push_back(tmp_val);
      }
    }
  }
  if (lantb_tmp.size() != lantb_dim0) {
    std::cerr << "lantb_tmp.size() != lantb_dim0 * lantb_dim1 * lantb_dim2\n";
    return false;
  }
  for (auto& lantb_row : lantb_tmp) {
    if (lantb_row.size() != lantb_dim1) {
      std::cerr << "lantb_row.size() != lantb_dim1\n";
      return false;
    }
    for (auto& lantb_t0 : lantb_row) {
      if (lantb_t0.size() != lantb_dim2) {
        std::cerr << "lantb_t0.size() != lantb_dim2\n";
        return false;
      }
    }
  }
  return true;
}

/// Prints a failure message when trying to load the FRS
/// \param f_name the FRS file name that failed
void FrsFailureMsg(const std::string& f_name) {
  std::cerr << "Failed to load FRS '" << f_name << "'\n";
}

bool LoadU0Intervals(std::istream& f_in, std::vector<Interval>& u0_intervals) {
  std::string in_str = "";
  int u0_rows = 0;
  int u0_cols = 0;
  f_in >> in_str >> u0_rows >> u0_cols;
  if (in_str != "u0_intervals") {
    std::cerr << "in_str != \"u0_intervals\"\n";
    return false;
  }
  if (u0_rows != 2) {
    std::cerr << "u0_rows != 2\n";
    return false;
  }
  u0_intervals.resize(u0_cols);
  for (int i = 0; i < u0_cols; ++i) {
    double tmp_val_min;
    f_in >> tmp_val_min;
    u0_intervals.at(i).SetMin(tmp_val_min);
  }
  for (int i = 0; i < u0_cols; ++i) {
    double tmp_val_max;
    f_in >> tmp_val_max;
    u0_intervals.at(i).SetMax(tmp_val_max);
  }
  return true;
}

} // namespace

FrsTotal LoadFrs(const std::string& f_name) {
  std::cout << "Loading FRS '" << f_name << "'...\n";
  std::ifstream f_in(f_name);
  std::string in_str;
  IndexT num_megas = 0;
  FrsTotal frs_total;
  frs_total.successful_ = true;
  f_in >> in_str >> frs_total.u0_min_ >> frs_total.u0_max_;
  if (in_str != "MINMAX_U0") {
    std::cerr << "in_str != \"MINMAX_U0\"\n";
    frs_total.successful_ = false;
    return frs_total;
  }
  f_in >> in_str >> frs_total.alternating_au_;
  if (in_str != "ALTERNATING_AU") {
    std::cerr << "in_str != \"ALTERNATING_AU\"\n";
    frs_total.successful_ = false;
    return frs_total;
  }
  // Set u0's
  f_in >> in_str >> num_megas;
  if (in_str != "NUM_MEGAS") {
    std::cerr << "in_str != \"NUM_MEGAS\"\n";
    frs_total.successful_ = false;
    return frs_total;
  }
  std::cout << "NUM MEGAS: " << num_megas << std::endl;
  if (num_megas <= 0) {
    std::cerr << "num_megas <= 0" << std::endl;
    frs_total.successful_ = false;
    FrsFailureMsg(f_name);
    return frs_total;
  }
  frs_total.megas_.resize(num_megas);
  for (IndexT i = 0; i < num_megas; ++i) {
    if (!LoadMega(f_in, frs_total.megas_.at(i))) {
      std::cerr << "LoadMega failed\n";
      frs_total.successful_ = false;
      FrsFailureMsg(f_name);
      return frs_total;
    }
  }
  if (!LoadU0Intervals(f_in, frs_total.u0_intervals_)) {
    std::cerr << "Load u0 intervals failed.\n";
    frs_total.successful_ = false;
    return frs_total;
  }
  f_in.close();
  if (frs_total.u0_intervals_.size() != frs_total.megas_.size()) {
    std::cerr << "frs_total.u0_intervals_.size() ["
              << frs_total.u0_intervals_.size()
              << "] != frs_total.megas_.size() [" << frs_total.megas_.size()
              << "]\n";
    FrsFailureMsg(f_name);
    frs_total.successful_ = false;
    return frs_total;
  }
  if (frs_total.successful_) {
    std::cout << "Successfully loaded '" << f_name << "'\n";
    std::cout << "FRS MEGAS: " << frs_total.megas_.size() << std::endl;
  } else {
    FrsFailureMsg(f_name);
  }
  return frs_total;
}

} // namespace roahm
