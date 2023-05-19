#ifndef ROAHM_RISK_PDF_H_
#define ROAHM_RISK_PDF_H_

#include <ros/console.h>
#include <unistd.h>

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "rover_state.hpp" // for ComplexHeadin...

namespace roahm {
struct PDF {
  double t_0;
  double t_1;
  double mu_0;
  double mu_1;
  double sigma_00;
  double sigma_01;
  double sigma_11;
};

class Risk_PDF {
  std::vector<std::vector<PDF>> pdf_for_all_maneuvers;
  int scenario;

public:
  bool readPDFDataFromCSV(std::string filename) {
    std::ifstream readFile(filename);

    std::vector<PDF> pdfs;
    PDF pdf;

    if (!readFile) {
      ROS_ERROR("PDF Trajectory File %s Read Error!", filename.c_str());
      return false;
    }

    std::string line, word;
    std::vector<double> dataLine;

    // stringstream str(line);

    while (getline(readFile, line)) {
      std::stringstream str(line);
      dataLine.clear();

      while (getline(str, word, ','))
        dataLine.push_back(stod(word));

      pdf.t_0 = dataLine[0];
      pdf.t_1 = dataLine[1];
      pdf.mu_0 = dataLine[2];
      pdf.mu_1 = dataLine[3];
      pdf.sigma_00 = dataLine[4];
      pdf.sigma_01 = dataLine[5];
      pdf.sigma_11 = dataLine[6];
      pdfs.push_back(pdf);
    }

    pdf_for_all_maneuvers.push_back(pdfs);
    ROS_INFO("%ld pdf data points loaded from %s.", pdfs.size(),
             filename.c_str());

    readFile.close();
    return true;
  }

  void GetStoredPDFs(const std::string& filepath, int scenario) {
    this->scenario = scenario;
    std::string full_file_name;
    int delta_t_start = 10;
    int delta_t_end = 160;
    int delta_t = 10;

    // PDF Files are stored in [1, .., n] indexed manner
    scenario += 1;

    for (int delta = delta_t_start; delta <= delta_t_end; delta += delta_t) {
      full_file_name = filepath + "/" + "scenario_" + std::to_string(scenario) +
                       "_delta_" + std::to_string(delta) + "ms.csv";
      readPDFDataFromCSV(full_file_name);
    }
    ROS_INFO("%ld saved pdf trajectories loaded for scenario %d",
             pdf_for_all_maneuvers.size(), scenario);
  }

  void matrixMultiply(std::vector<std::vector<double>>& A,
                      std::vector<std::vector<double>>& B,
                      std::vector<std::vector<double>>& AB) {
    for (int i = 0; i < A.size(); ++i) {
      for (int k = 0; k < A[0].size(); ++k) {
        for (int j = 0; j < B[0].size(); ++j) {
          AB[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }

  void pushBackTransformedPdfToDoubleVector(
      std::vector<double>& mu_sigma, const PDF& pdf,
      const RoverState& tracked_rover_state, const bool first_pdf,
      double& determinant_normalization_factor) {
    double heading_global_frame;
    double pdf_mu_x_global_frame;
    double pdf_mu_y_global_frame;
    std::vector<std::vector<double>> pdf_sigma(2, std::vector<double>(2));
    std::vector<std::vector<double>> rotation_matrix(2, std::vector<double>(2));

    // Transform Mu in the global frame
    heading_global_frame = tracked_rover_state.GetHeading();

    pdf_mu_x_global_frame = std::cos(heading_global_frame) * pdf.mu_0 -
                            std::sin(heading_global_frame) * pdf.mu_1 +
                            tracked_rover_state.x_;

    pdf_mu_y_global_frame = std::sin(heading_global_frame) * pdf.mu_0 +
                            std::cos(heading_global_frame) * pdf.mu_1 +
                            tracked_rover_state.y_;

    // Normalize PDF before Rotating
    double default_x_variance = 0.01;
    double default_y_variance = 0.01;

    double determinant_pdf_sigma =
        (pdf.sigma_00 * pdf.sigma_11) - (pdf.sigma_01 * pdf.sigma_01);
    double determinant_default_pdf_sigma =
        default_x_variance * default_y_variance;

    if (first_pdf) {
      pdf_sigma = {{default_x_variance, 0.0}, {0.0, default_y_variance}};

      determinant_normalization_factor =
          determinant_pdf_sigma / determinant_default_pdf_sigma;
    } else {
      if (determinant_normalization_factor == 0.0) {
        ROS_ERROR("Determinant Normalization Factor is 0.0");
      }
      pdf_sigma[0][0] = pdf.sigma_00 / determinant_normalization_factor;
      pdf_sigma[0][1] = pdf.sigma_01 / determinant_normalization_factor;
      pdf_sigma[1][0] = pdf_sigma[0][1];
      pdf_sigma[1][1] = pdf.sigma_11 / determinant_normalization_factor;
    }

    // Transform Sigma in the global frame
    // Rotation Matrix ^ T
    rotation_matrix = {
        {std::cos(heading_global_frame), std::sin(heading_global_frame)},
        {-std::sin(heading_global_frame), std::cos(heading_global_frame)}};

    std::vector<std::vector<double>> temp_matrix(2,
                                                 std::vector<double>(2, 0.0));
    matrixMultiply(rotation_matrix, pdf_sigma, temp_matrix);

    // Rotation Matrix
    rotation_matrix = {
        {std::cos(heading_global_frame), -std::sin(heading_global_frame)},
        {std::sin(heading_global_frame), std::cos(heading_global_frame)}};

    pdf_sigma = {{0.0, 0.0}, {0.0, 0.0}};

    matrixMultiply(temp_matrix, rotation_matrix, pdf_sigma);

    mu_sigma.push_back(pdf_mu_x_global_frame);
    mu_sigma.push_back(pdf_mu_y_global_frame);
    mu_sigma.push_back(pdf_sigma[0][0]);
    mu_sigma.push_back(pdf_sigma[0][1]);
    mu_sigma.push_back(pdf_sigma[1][1]);
  }

  // Returns PDFs [t0, mu, sigma] for a given [scenario_clock_time,
  // scenario_clock_time + pdf_duration]
  std::array<std::vector<double>, 16>
  GetMuSigmas(const RoverState& tracked_rover_state,
              const double scenario_clock_time, const double pdf_duration) {
    std::array<std::vector<double>, 16> mu_sigma_for_all_deltas;

    // ROS_INFO("getManeuverPdf function called with t0 = %f for scenario = %d",
    // ego_rover_t0, this->scenario);

    double last_t, delta_t;        // secs
    double pdf_duration_remaining; // secs
    std::vector<double> mu_sigma;
    bool first_pdf;

    int max_pdf_delta_intervals = 16; // 10ms to 160ms

    for (int pdf_delta_idx = 0; pdf_delta_idx < max_pdf_delta_intervals;
         pdf_delta_idx++) {
      pdf_duration_remaining = pdf_duration;
      mu_sigma.clear();
      delta_t = (pdf_delta_idx + 1) * 0.01; // (File idx + 1) * 10ms
      first_pdf = true;
      double determinant_normalization_factor = 0.0;

      for (auto& pdf : pdf_for_all_maneuvers[pdf_delta_idx]) {
        last_t = pdf.t_0;

        // Append pdfs [scenario_clock_time, scenario_clock_time + pdf_duration]
        if ((pdf.t_0 >= scenario_clock_time) &&
            (pdf.t_1 <= (scenario_clock_time + pdf_duration))) {
          pushBackTransformedPdfToDoubleVector(
              mu_sigma, pdf, tracked_rover_state, first_pdf,
              determinant_normalization_factor);
          pdf_duration_remaining -= delta_t;
          first_pdf = false;
        }
        // Accept pdfs within delta_t interval of scenario_start_time
        else if (pdf.t_0 < scenario_clock_time &&
                 pdf.t_1 > scenario_clock_time) {
          pushBackTransformedPdfToDoubleVector(
              mu_sigma, pdf, tracked_rover_state, first_pdf,
              determinant_normalization_factor);
          pdf_duration_remaining -= delta_t;
          first_pdf = false;
        } else if (pdf.t_0 > (scenario_clock_time + pdf_duration)) {
          break;
        }
      }

      // While loop to return mu = (1000, 1000) sigma = [[0.1, 0], [0, 0.1]]
      // values if pdf_duration_remaining > 0
      PDF pdf;
      pdf.t_0 = last_t + delta_t;
      pdf.mu_0 = 1000;
      pdf.mu_1 = 1000;
      pdf.sigma_00 = 0.01;
      pdf.sigma_01 = 0.0;
      pdf.sigma_01 = 0.0;
      pdf.sigma_11 = 0.01;
      first_pdf = true;

      while ((pdf_duration_remaining - delta_t) > 0.0) {
        pushBackTransformedPdfToDoubleVector(mu_sigma, pdf, tracked_rover_state,
                                             first_pdf,
                                             determinant_normalization_factor);
        pdf.t_0 += delta_t;
        pdf_duration_remaining -= delta_t;
      }

      mu_sigma_for_all_deltas.at(pdf_delta_idx) = mu_sigma;
    }
    return mu_sigma_for_all_deltas;
  }
};
} // namespace roahm

#endif // ROAHM_RISK_PDF_H_