#include "frs_mega.hpp"

#include "frs_io.hpp"

namespace roahm {

const std::vector<std::vector<Vehrs>>&
FrsMega::GetDirLanSet(const ManuType manu_type) const {
  return manu_type == ManuType::kDirChange ? dir_ : lan_;
}

bool FrsMega::operator==(const FrsMega& rhs) const {
  return VecsEq(au_, rhs.au_, "Mega::Au") &
         VecsEq(dir_, rhs.dir_, "Mega::Dir") &
         VecsEq(lan_, rhs.lan_, "Mega::Lan");
}

// FrsMega
void FrsMega::WriteToBinFileDirect(std::ofstream& file) const {
  WriteToBinFile(file, au_);
  WriteToBinFile(file, dir_);
  WriteToBinFile(file, lan_);
}
FrsMega FrsMega::ReadFromBinFileDirect(std::ifstream& file) {
  FrsMega ret;
  ReadToVarFromBinFile(ret.au_, file);
  ReadToVarFromBinFile(ret.dir_, file);
  ReadToVarFromBinFile(ret.lan_, file);
  return ret;
}

} // namespace roahm