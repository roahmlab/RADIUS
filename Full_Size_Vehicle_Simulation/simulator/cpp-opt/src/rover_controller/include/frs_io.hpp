#ifndef ROAHM_FRS_IO_HPP_
#define ROAHM_FRS_IO_HPP_
#include <fstream>
#include <vector>

#include "frs_loader.hpp"

namespace roahm {
// Generic
template <typename T> void WriteToBinFile(std::ofstream& file, const T& data);

template <typename T> T ReadFromBinFile(std::ifstream& file);

template <typename T>
void ReadToVarFromBinFile(T& out_param, std::ifstream& file) {
  out_param = ReadFromBinFile<T>(file);
}

// double
template <> void WriteToBinFile(std::ofstream& file, const double& data);

template <> double ReadFromBinFile(std::ifstream& file);

// int
template <> void WriteToBinFile(std::ofstream& file, const int& data);

template <> int ReadFromBinFile(std::ifstream& file);

// bool
template <> void WriteToBinFile(std::ofstream& file, const bool& data);

template <> bool ReadFromBinFile(std::ifstream& file);

// Interval
template <> void WriteToBinFile(std::ofstream& file, const Interval& data);

template <> Interval ReadFromBinFile(std::ifstream& file);

// std::vector<Interval>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Interval>& data);

template <> std::vector<Interval> ReadFromBinFile(std::ifstream& file);

// std::vector<double>
template <> std::vector<double> ReadFromBinFile(std::ifstream& file);

template <>
void WriteToBinFile(std::ofstream& file, const std::vector<double>& data);

// std::vector<int>
template <> std::vector<int> ReadFromBinFile(std::ifstream& file);

template <>
void WriteToBinFile(std::ofstream& file, const std::vector<int>& data);

// std::vector<std::vector<double>> (for Autb)
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<double>>& data);

template <>
std::vector<std::vector<double>> ReadFromBinFile(std::ifstream& file);

// std::vector<std::vector<std::vector<double>>> (for dirtb/lantb)
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<std::vector<double>>>& data);

template <>
std::vector<std::vector<std::vector<double>>>
ReadFromBinFile(std::ifstream& file);

// Sliceable
template <> void WriteToBinFile(std::ofstream& file, const Sliceable& data);

template <> Sliceable ReadFromBinFile(std::ifstream& file);

// std::vector<Sliceable>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Sliceable>& data);

template <> std::vector<Sliceable> ReadFromBinFile(std::ifstream& file);

// ZonoSliceInfo
template <> void WriteToBinFile(std::ofstream& file, const ZonoSliceInfo& data);

template <> ZonoSliceInfo ReadFromBinFile(std::ifstream& file);

// std::vector<ZonoSliceInfo>
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<ZonoSliceInfo>& data);

template <> std::vector<ZonoSliceInfo> ReadFromBinFile(std::ifstream& file);

// std::vector<bool>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<bool>& data);

template <> std::vector<bool> ReadFromBinFile(std::ifstream& file);

// std::vector<Interval>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Interval>& data);

template <> std::vector<Interval> ReadFromBinFile(std::ifstream& file);

// std::size_t
template <> void WriteToBinFile(std::ofstream& file, const std::size_t& data);

template <> std::size_t ReadFromBinFile(std::ifstream& file);

// std::vector<std::size_t>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<std::size_t>& data);

template <> std::vector<std::size_t> ReadFromBinFile(std::ifstream& file);

// Vehrs
template <> void WriteToBinFile(std::ofstream& file, const Vehrs& data);

template <> Vehrs ReadFromBinFile(std::ifstream& file);

// std::vector<Vehrs>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<Vehrs>& data);

template <> std::vector<Vehrs> ReadFromBinFile(std::ifstream& file);

// std::vector<std::vector<Vehrs>>
template <>
void WriteToBinFile(std::ofstream& file,
                    const std::vector<std::vector<Vehrs>>& data);

template <>
std::vector<std::vector<Vehrs>> ReadFromBinFile(std::ifstream& file);

// FrsMega
template <> void WriteToBinFile(std::ofstream& file, const FrsMega& data);

template <> FrsMega ReadFromBinFile(std::ifstream& file);

// std::vector<FrsMega>
template <>
void WriteToBinFile(std::ofstream& file, const std::vector<FrsMega>& data);

template <> std::vector<FrsMega> ReadFromBinFile(std::ifstream& file);

// FrsTotal
template <> void WriteToBinFile(std::ofstream& file, const FrsTotal& data);

template <> FrsTotal ReadFromBinFile(std::ifstream& file);

// CudaInfo
template <> void WriteToBinFile(std::ofstream& file, const CudaInfo& data);

template <> CudaInfo ReadFromBinFile(std::ifstream& file);

} // namespace roahm
#endif // ROAHM_FRS_IO_HPP_