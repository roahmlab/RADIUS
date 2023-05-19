#include "frs_io.hpp"

namespace roahm {

//
// bool
//
template <> void WriteToBinFile(std::ofstream& file, const bool& data) {
  file.write(reinterpret_cast<const char*>(&data), sizeof(bool));
}

template <> bool ReadFromBinFile(std::ifstream& file) {
  bool data;
  file.read(reinterpret_cast<char*>(&data), sizeof(bool));
  return data;
}

//
// int
//
template <> void WriteToBinFile(std::ofstream& file, const int& data) {
  file.write(reinterpret_cast<const char*>(&data), sizeof(int));
}

template <> int ReadFromBinFile(std::ifstream& file) {
  int data;
  file.read(reinterpret_cast<char*>(&data), sizeof(int));
  return data;
}

//
// std::size_t
//
template <> void WriteToBinFile(std::ofstream& file, const std::size_t& data) {
  file.write(reinterpret_cast<const char*>(&data), sizeof(std::size_t));
}

template <> std::size_t ReadFromBinFile(std::ifstream& file) {
  std::size_t data;
  file.read(reinterpret_cast<char*>(&data), sizeof(std::size_t));
  return data;
}

//
// double
//
template <> void WriteToBinFile(std::ofstream& file, const double& data) {
  file.write(reinterpret_cast<const char*>(&data), sizeof(double));
}

template <> double ReadFromBinFile(std::ifstream& file) {
  double data;
  file.read(reinterpret_cast<char*>(&data), sizeof(double));
  return data;
}

//
// ::roahm::Interval
//
template <> void WriteToBinFile(std::ofstream& file, const Interval& data) {
  WriteToBinFile(file, data.Min());
  WriteToBinFile(file, data.Max());
}

template <> Interval ReadFromBinFile(std::ifstream& file) {
  const double data_min = ReadFromBinFile<double>(file);
  const double data_max = ReadFromBinFile<double>(file);
  return Interval{data_min, data_max};
}

//
// ::roahm::Sliceable
//
template <> void WriteToBinFile(std::ofstream& file, const Sliceable& data) {
  data.WriteToBinFileDirect(file);
}

template <> Sliceable ReadFromBinFile(std::ifstream& file) {
  return Sliceable::ReadFromBinFileDirect(file);
}

// ZonoSliceInfo
template <>
void WriteToBinFile(std::ofstream& file, const ZonoSliceInfo& data) {
  WriteToBinFile(file, data.slc_vals_);
}

template <> ZonoSliceInfo ReadFromBinFile(std::ifstream& file) {
  ZonoSliceInfo data;
  data.slc_vals_ = ReadFromBinFile<std::vector<Sliceable>>(file);
  return data;
}

//
// std::vector<bool>
//
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<bool>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<bool> ReadFromBinFile(std::ifstream& file) {
  std::vector<bool> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<bool>(file));
  }
  return data;
}

//
// std::vector<int>
//
template <> std::vector<int> ReadFromBinFile(std::ifstream& file) {
  std::vector<int> data;
  int size = ReadFromBinFile<int>(file);
  data.resize(size);
  file.read(reinterpret_cast<char*>(data.data()), size * sizeof(int));
  return data;
}

template <>
void WriteToBinFile(std::ofstream& file, const std::vector<int>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  file.write(reinterpret_cast<const char*>(data.data()),
             data.size() * sizeof(int));
}

//
// std::vector<std::size_t>
//
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<std::size_t>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  file.write(reinterpret_cast<const char*>(data.data()),
             data.size() * sizeof(std::size_t));
}

template <> std::vector<std::size_t> ReadFromBinFile(std::ifstream& file) {
  std::vector<std::size_t> data;
  int size = ReadFromBinFile<int>(file);
  data.resize(size);
  file.read(reinterpret_cast<char*>(data.data()), size * sizeof(std::size_t));
  return data;
}

//
// std::vector<double>
//
template <> std::vector<double> ReadFromBinFile(std::ifstream& file) {
  std::vector<double> data;
  int size = ReadFromBinFile<int>(file);
  data.resize(size);
  file.read(reinterpret_cast<char*>(data.data()), size * sizeof(double));
  return data;
}

template <>
void WriteToBinFile(std::ofstream& file, const std::vector<double>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  file.write(reinterpret_cast<const char*>(data.data()),
             data.size() * sizeof(double));
}

//
// std::vector<Interval>
//
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Interval>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<Interval> ReadFromBinFile(std::ifstream& file) {
  std::vector<Interval> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<Interval>(file));
  }
  return data;
}

//
// std::vector<std::vector<double>>
//
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<double>>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <>
std::vector<std::vector<double>> ReadFromBinFile(std::ifstream& file) {
  std::vector<std::vector<double>> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<std::vector<double>>(file));
  }
  return data;
}

//
// std::vector<std::vector<std::vector<double>>>
//
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<std::vector<double>>>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <>
std::vector<std::vector<std::vector<double>>>
ReadFromBinFile(std::ifstream& file) {
  std::vector<std::vector<std::vector<double>>> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<std::vector<std::vector<double>>>(file));
  }
  return data;
}

//
// std::vector<Sliceable>
//
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Sliceable>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<Sliceable> ReadFromBinFile(std::ifstream& file) {
  std::vector<Sliceable> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<Sliceable>(file));
  }
  return data;
}

//==============================================================================
// UNTESTED
//==============================================================================

template <> void WriteToBinFile(std::ofstream& file, const CudaInfo& data) {
  data.WriteToBinFileDirect(file);
}

template <> CudaInfo ReadFromBinFile(std::ifstream& file) {
  return CudaInfo::ReadFromBinFileDirect(file);
}

// std::vector<ZonoSliceInfo>
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<ZonoSliceInfo>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<ZonoSliceInfo> ReadFromBinFile(std::ifstream& file) {
  std::vector<ZonoSliceInfo> data;
  int size = ReadFromBinFile<int>(file);
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<ZonoSliceInfo>(file));
  }
  return data;
}

template <> void WriteToBinFile(std::ofstream& file, const Vehrs& data) {
  data.WriteToBinFileDirect(file);
}

template <> Vehrs ReadFromBinFile(std::ifstream& file) {
  return Vehrs::ReadFromBinFileDirect(file);
}

template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Vehrs>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<Vehrs> ReadFromBinFile(std::ifstream& file) {
  int size = ReadFromBinFile<int>(file);
  std::vector<Vehrs> data;
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<Vehrs>(file));
  }
  return data;
}

template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<Vehrs>>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <>
std::vector<std::vector<Vehrs>> ReadFromBinFile(std::ifstream& file) {
  int size = ReadFromBinFile<int>(file);
  std::vector<std::vector<Vehrs>> data;
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<std::vector<Vehrs>>(file));
  }
  return data;
}

template <> void WriteToBinFile(std::ofstream& file, const FrsMega& data) {
  data.WriteToBinFileDirect(file);
}

template <> FrsMega ReadFromBinFile(std::ifstream& file) {
  return FrsMega::ReadFromBinFileDirect(file);
}

// std::vector<FrsMega>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<FrsMega>& data) {
  WriteToBinFile(file, static_cast<int>(data.size()));
  for (const auto& d : data) {
    WriteToBinFile(file, d);
  }
}

template <> std::vector<FrsMega> ReadFromBinFile(std::ifstream& file) {
  int size = ReadFromBinFile<int>(file);
  std::vector<FrsMega> data;
  data.reserve(size);
  for (int i = 0; i < size; ++i) {
    data.push_back(ReadFromBinFile<FrsMega>(file));
  }
  return data;
}

// FrsTotal
template <> void WriteToBinFile(std::ofstream& file, const FrsTotal& data) {
  data.WriteToBinFileDirect(file);
}

template <> FrsTotal ReadFromBinFile(std::ifstream& file) {
  return FrsTotal::ReadFromBinFileDirect(file);
}
} // namespace roahm