#include "frs_io.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>

#include "frs_text_loader.hpp"
#include "timing_util.hpp"

TEST(FrsLoaderTest, LoadFrs) {
  using ::roahm::GetDeltaS;
  using ::roahm::Tick;

  const auto t0 = Tick();
  const ::roahm::FrsTotal frs{::roahm::LoadFrs("/data/car_frs.txt")};
  std::cout << "FRS Num Megas: " << frs.megas_.size() << std::endl;
  const auto t1 = Tick();
  const std::string file_name{"/data/frs_bin_no_cuda.txt"};
  std::ofstream file_out(file_name, std::ios::binary);
  const auto t2 = Tick();
  WriteToBinFile(file_out, frs);
  const auto t3 = Tick();
  file_out.close();
  const auto t4 = Tick();

  std::ifstream file_in(file_name, std::ios::binary);
  const auto t5 = Tick();
  ::roahm::FrsTotal frs2{::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)};
  const auto t6 = Tick();
  file_in.close();
  const auto t7 = Tick();
  std::cout << "Load TXT:  " << GetDeltaS(t1, t0) << " s" << std::endl;
  std::cout << "Write BIN: " << GetDeltaS(t3, t2) << " s" << std::endl;
  std::cout << "Load BIN:  " << GetDeltaS(t6, t5) << " s" << std::endl;
  std::cout << "Equality: " << (frs == frs2) << std::endl;
  EXPECT_EQ(frs, frs2);
}

class BinaryReadWrite : public ::testing::Test {
private:
  std::string file_name_;
  std::optional<std::ofstream> file_out_;
  std::optional<std::ifstream> file_in_;

  void CloseFiles() {
    if (file_in_.has_value()) {
      file_in_.value().close();
      file_in_ = std::nullopt;
    }
    if (file_out_.has_value()) {
      file_out_.value().close();
      file_out_ = std::nullopt;
    }
  }

protected:
  void SetUp(std::string file_name) {
    file_name_ = file_name;
    file_out_ = std::nullopt;
    file_in_ = std::nullopt;
  }

  void TearDown() { CloseFiles(); }

  std::ofstream& GetFileOut() {
    CloseFiles();
    file_out_ = std::ofstream(file_name_, std::ios::binary);
    return file_out_.value();
  }

  std::ifstream& GetFileIn() {
    CloseFiles();
    file_in_ = std::ifstream(file_name_, std::ios::binary);
    return file_in_.value();
  }

  template <typename T> void BinaryWrite(const std::vector<T>& data) {
    auto& file_out = GetFileOut();
    for (const auto& d : data) {
      ::roahm::WriteToBinFile<T>(file_out, d);
    }
  }
  template <typename T> std::vector<T> BinaryRead(int num) {
    std::vector<T> data;
    auto& file_in = GetFileIn();
    for (int i = 0; i < num; ++i) {
      data.push_back(::roahm::ReadFromBinFile<T>(file_in));
    }
    return data;
  }
};

TEST_F(BinaryReadWrite, Bool) {
  SetUp("/data/test_out_bool.bin");
  using T = bool;
  const std::vector<T> vals_in{false, true, true, false, false};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, Int) {
  SetUp("/data/test_out_int.bin");
  using T = int;
  const std::vector<T> vals_in{4, -13};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, SizeT) {
  SetUp("/data/test_out_size_t.bin");
  using T = std::size_t;
  constexpr auto kSizeTMax = std::numeric_limits<std::size_t>::max();
  const std::vector<T> vals_in{4, 13, 193, kSizeTMax, 18};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, Double) {
  SetUp("/data/test_out_double.bin");
  using T = double;
  const std::vector<T> vals_in{4, -13, 15.03, -19.05};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, Interval) {
  SetUp("/data/test_out_interval.bin");
  using T = ::roahm::Interval;
  const std::vector<T> vals_in{{0.0, 1.0}, {-13, 15.03}, {-19.05, 18.0}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, Sliceable) {
  SetUp("/data/test_out_sliceable.bin");
  using T = ::roahm::Sliceable;
  const std::vector<T> vals_in{
      {0, 14.0, 0.01, -18.0, 13.0, 19.0},
      {18, -1.0, 5.0, 17.3, 13.9, 19.1},
  };
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, ZonoSliceInfo) {
  SetUp("/data/test_out_zono_slice_info.bin");
  using T = ::roahm::ZonoSliceInfo;
  const std::vector<T> vals_in{
      T{{
          ::roahm::Sliceable{0, 14.0, 0.01, -18.0, 13.0, 19.0},
          ::roahm::Sliceable{18, -1.0, 5.0, 17.3, 13.9, 19.1},
      }},
      T{},
      T{{
          ::roahm::Sliceable{1, 14.0, 0.03, -18.0, 13.0, 78.0},
          ::roahm::Sliceable{3, -19.0, 3.0, 19.71, -193.9, 70.1},
          ::roahm::Sliceable{3, -1.0, -2.0, -17.3, 21.9, 56.1},
      }}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecBool) {
  SetUp("/data/test_out_vec_bool.bin");
  using T = std::vector<bool>;
  const std::vector<T> vals_in{
      {false}, {true, false, true}, {true}, {true, true, false, true}, {}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecInt) {
  SetUp("/data/test_out_vec_int.bin");
  using T = std::vector<int>;
  const std::vector<T> vals_in{
      {0, 1}, {-13, 15, 19}, {-21}, {14, 15, 16, 17}, {}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecSizeT) {
  SetUp("/data/test_out_vec_size_t.bin");
  using T = std::vector<std::size_t>;
  constexpr auto kSizeTMax = std::numeric_limits<std::size_t>::max();
  const std::vector<T> vals_in{
      {0, kSizeTMax}, {13, 15, 19}, {21}, {14, 15, 16, 17}, {}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecDouble) {
  SetUp("/data/test_out_vec_double.bin");
  using T = std::vector<double>;
  const std::vector<T> vals_in{
      {0.5, 1.1}, {-13.2, 15.1991, 19.43}, {-21.0}, {14.0, 15.3, 16.4, 17.2}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecInterval) {
  SetUp("/data/test_out_vec_interval.bin");
  using T = std::vector<::roahm::Interval>;
  const std::vector<T> vals_in{{{0.5, 1.1}, {19.0, 34.0}},
                               {{-13.2, 15.1991}},
                               {{19.43, 156.7}, {-175.0, 14.3}, {-2.0, 0.0}},
                               {},
                               {{-21.0, 14.0}, {15.3, 16.4}, {17.2, 14.0}}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecSliceable) {
  SetUp("/data/test_out_vec_sliceable.bin");
  using T = std::vector<::roahm::Sliceable>;
  const std::vector<T> vals_in{{
                                   {0, 14.0, 0.01, -18.0, 13.0, 19.0},
                                   {18, -1.0, 5.0, 17.3, 13.9, 19.1},
                               },
                               {},
                               {
                                   {1, 14.0, 0.03, -18.0, 13.0, 78.0},
                                   {3, -19.0, 3.0, 19.71, -193.9, 70.1},
                                   {3, -1.0, -2.0, -17.3, 21.9, 56.1},
                               }};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecVecDouble) {
  SetUp("/data/test_out_vec_vec_double.bin");
  using T = std::vector<std::vector<double>>;
  const std::vector<T> vals_in{{{0.5, 1.1}, {-13.2, 15.1991, 19.43}},
                               {{}, {-21.0}, {}, {14.0, 15.3, 16.4, 17.2}},
                               {}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

TEST_F(BinaryReadWrite, VecVecVecDouble) {
  SetUp("/data/test_out_vec_vec_vec_double.bin");
  using T = std::vector<std::vector<std::vector<double>>>;
  const std::vector<T> vals_in{{{{0.5, 1.1}, {-13.2, 15.1991, 19.43}},
                                {{}, {-21.0}, {}, {14.0, 15.3, 16.4, 17.2}},
                                {}},
                               {{{0.1}, {-1300.2, 1.12, -16.43}, {}, {1.4}},
                                {{}, {-21.0}, {}, {14.0, 15.3, 16.4, 17.2}},
                                {{1.67, 1.91}},
                                {},
                                {{1.31}, {1.43, 1.5, 1.9}, {18.4}}}};
  BinaryWrite<T>(vals_in);
  const auto vals_out = BinaryRead<T>(vals_in.size());
  EXPECT_EQ(vals_in, vals_out);
}

/*

TEST(FrsLoaderTest, ImportExportTest) {
  const std::string frs_in_fname{"/data/frs_bin_processed.txt"};
  const std::string frs_out_fname{"/data/frs_bin_processed_out_test.txt"};

  std::ifstream file_in(frs_in_fname, std::ios::binary);
  const auto t0 = Tick();
  ::roahm::FrsTotal
frs_in{::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in)}; const auto t1 =
Tick(); file_in.close(); std::cout << "Load IN:   " << GetDeltaS(t1, t0) << " s"
<< std::endl;

  std::ofstream file_out(frs_out_fname, std::ios::binary);
  const auto t2 = Tick();
  ::roahm::WriteToBinFile(file_out, frs_in);
  const auto t3 = Tick();
  file_out.close();
  std::cout << "Write OUT: " << GetDeltaS(t3, t2) << " s" << std::endl;

  std::ifstream file_out_in(frs_out_fname, std::ios::binary);
  const auto t4 = Tick();
  ::roahm::FrsTotal
frs_out{::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_out_in)}; const auto t5
= Tick(); file_out_in.close(); std::cout << "Load OUT:  " << GetDeltaS(t5, t4)
<< " s" << std::endl;

  const auto t6 = Tick();
  EXPECT_EQ(frs_in, frs_out);
  const auto t7 = Tick();
  std::cout << "Compare IN/OUT:  " << GetDeltaS(t7, t6) << " s" << std::endl;
}
*/

/*
TEST(FrsLoaderTest, CompareFRS) {
  const auto t_begin = Tick();
  const std::string matlab_exported_frs_fname{"/data/frs_bin_processed.txt"};
  const std::string cpp_exported_frs_fname{"/data/frs_bin_no_cuda.txt"};

  std::ifstream file_in_mat(matlab_exported_frs_fname, std::ios::binary);
  ::roahm::FrsTotal
frs_mat{::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in_mat)};
  file_in_mat.close();

  const auto t0 = Tick();
  std::ifstream file_in_cpp(cpp_exported_frs_fname, std::ios::binary);
  ::roahm::FrsTotal
frs_cpp{::roahm::ReadFromBinFile<::roahm::FrsTotal>(file_in_cpp)};
  file_in_cpp.close();
  const auto t1 = Tick();

  constexpr double kTol = 1.0e-10;
  EXPECT_TRUE(frs_mat.NearEqual(frs_cpp, kTol));
  const auto t_end = Tick();
  std::cout << "Time to read FRS: " << GetDeltaS(t1, t0) << " s" << std::endl;
  std::cout << "Time TOTAL: " << GetDeltaS(t_end, t_begin) << " s" << std::endl;
}
*/
