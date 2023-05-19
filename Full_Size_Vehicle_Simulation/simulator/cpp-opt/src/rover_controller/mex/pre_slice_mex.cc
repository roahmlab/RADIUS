#include "mex.h"
#include "pre_slice.hpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // [x0, y0, H1, H2, H4, t] = pre_slice_n_IA(num_zono, x0, y0, dx, dy,
  // u0v0r0_slice_beta, ...
  //                                          g_u0_x, g_u0_y, g_v0_x, g_v0_y,
  //                                          g_r0_x, g_r0_y,...
  //                                          block_inzono_list, rot_angle,
  //                                          mu_sigma,... g_p_x, grid_size)
  // output t gives the computation time (no memory allocation/copy/free
  // included) input = 0: number of zonotopes in this FRS serie
  //         1: [x0_zono1, x0_zono2, ...]
  //         2: [y0_zono1, y0_zono2, ...]
  //         3: [dx_zono1, dx_zono2, ...]
  //         4: [dy_zono1, dy_zono2, ...]
  //         5: coefficients of the sliceable generator w.r.t. u0v0r0
  //         6: [g_u0_x_zono1, g_u0_x_zono2, ...]
  //         7: [g_u0_y_zono1, g_u0_y_zono2, ...]
  //         8: [g_v0_x_zono1, g_v0_x_zono2, ...]
  //         9: [g_v0_y_zono1, g_v0_y_zono2, ...]
  //         10: [g_r0_x_zono1, g_r0_x_zono2, ...]
  //         11: [g_r0_y_zono1, g_r0_y_zono2, ...]
  //         12: [block_inzono_list_zono1, block_inzono_list_zono2, ...] NOTE,
  //         'outzono' block is set to be -1 13: [rot_angle_zono1,
  //         rot_angle_zono2, ...] NOTE, this vector stores the angle s.t.
  //         p-sliceable generator is rotated to [xxx; 0] along its xy
  //         dimensions 14: [mu_x_zono1, mu_y_zono1, sigma_1_zono1,
  //         sigma_2_zono1, sigma_4_zono1, mu_x_zono2, mu_y_zono2,
  //         sigma_1_zono2, sigma_2_zono2, sigma_4_zono2, ...] 15: [g_p_x_zono1,
  //         g_p_x_zono2, ...] 16: grid_size
  constexpr int kExpectedNumInputs = 17;
  constexpr int kExpectedNumOutputs = 6;
  if (nrhs != kExpectedNumInputs) {
    throw std::runtime_error("Wrong number of inputs, " + std::to_string(nrhs) +
                             " != " + std::to_string(kExpectedNumInputs));
  }

  if (nlhs != kExpectedNumOutputs) {
    throw std::runtime_error("Wrong number of outputs, " +
                             std::to_string(nlhs) +
                             " != " + std::to_string(kExpectedNumOutputs));
  }

  // get input data
  // int grid_size =
  // static_cast<int>(*static_cast<double*>(mxGetData(prhs[16]))); int num_zono
  // = static_cast<int>(*static_cast<double*>(mxGetData(prhs[0]))); double*
  // host_x0 = static_cast<double*>(mxGetData(prhs[1])); double* host_y0 =
  // static_cast<double*>(mxGetData(prhs[2])); double* host_dx =
  // static_cast<double*>(mxGetData(prhs[3])); double* host_dy =
  // static_cast<double*>(mxGetData(prhs[4])); double* host_u0v0r0_slice_beta =
  // static_cast<double*>(mxGetData(prhs[5])); double* host_g_u0_x =
  // static_cast<double*>(mxGetData(prhs[6])); double* host_g_u0_y =
  // static_cast<double*>(mxGetData(prhs[7])); double* host_g_v0_x =
  // static_cast<double*>(mxGetData(prhs[8])); double* host_g_v0_y =
  // static_cast<double*>(mxGetData(prhs[9])); double* host_g_r0_x =
  // static_cast<double*>(mxGetData(prhs[10])); double* host_g_r0_y =
  // static_cast<double*>(mxGetData(prhs[11])); double* host_block_inzono_list =
  // static_cast<double*>(mxGetData(prhs[12])); // this is stupid, MATLAB has to
  // load as double?! int* doesn't work??? double* host_rot_angle =
  // static_cast<double*>(mxGetData(prhs[13])); double* host_mu_sigma =
  // static_cast<double*>(mxGetData(prhs[14])); double* host_g_p_x =
  // static_cast<double*>(mxGetData(prhs[15]));
  int grid_size = (int)(*(double*)(mxGetData(prhs[16])));
  int num_zono = (int)(*(double*)(mxGetData(prhs[0])));
  double* host_x0 = (double*)(mxGetData(prhs[1]));
  double* host_y0 = (double*)(mxGetData(prhs[2]));
  double* host_dx = (double*)(mxGetData(prhs[3]));
  double* host_dy = (double*)(mxGetData(prhs[4]));
  double* host_u0v0r0_slice_beta = (double*)(mxGetData(prhs[5]));
  double* host_g_u0_x = (double*)(mxGetData(prhs[6]));
  double* host_g_u0_y = (double*)(mxGetData(prhs[7]));
  double* host_g_v0_x = (double*)(mxGetData(prhs[8]));
  double* host_g_v0_y = (double*)(mxGetData(prhs[9]));
  double* host_g_r0_x = (double*)(mxGetData(prhs[10]));
  double* host_g_r0_y = (double*)(mxGetData(prhs[11]));
  double* host_block_inzono_list =
      (double*)(mxGetData(prhs[12])); // this is stupid, MATLAB has to load as
                                      // double?! int* doesn't work???
  double* host_rot_angle = (double*)(mxGetData(prhs[13]));
  double* host_mu_sigma = (double*)(mxGetData(prhs[14]));
  double* host_g_p_x = (double*)(mxGetData(prhs[15]));

  auto touch_all = [](double* ptr, int num, std::string str) {
    std::cout << "Touching all " << num << " elements of " << str << "..."
              << std::endl;
    for (int i = 0; i < num; ++i) {
      if (num <= 5) {
        std::cout << str << "[" << i << "] = " << ptr[i] << std::endl;
      } else {
        ptr[i] = ptr[i];
      }
      // ptr[i] = ptr[i];
    }
    std::cout << "Touched all elements of " << str << "." << std::endl;
  };

  touch_all(host_x0, num_zono, "host_x0");
  touch_all(host_y0, num_zono, "host_y0");
  touch_all(host_dx, num_zono, "host_dx");
  touch_all(host_dy, num_zono, "host_dy");
  touch_all(host_u0v0r0_slice_beta, 3, "host_u0v0r0_slice_beta");
  touch_all(host_g_u0_x, num_zono, "host_g_u0_x");
  touch_all(host_g_u0_y, num_zono, "host_g_u0_y");
  touch_all(host_g_v0_x, num_zono, "host_g_v0_x");
  touch_all(host_g_v0_y, num_zono, "host_g_v0_y");
  touch_all(host_g_r0_x, num_zono, "host_g_r0_x");
  touch_all(host_g_r0_y, num_zono, "host_g_r0_y");
  touch_all(host_block_inzono_list, num_zono * grid_size * grid_size,
            "host_block_inzono_list");
  touch_all(host_rot_angle, num_zono, "host_rot_angle");
  touch_all(host_mu_sigma, 5 * num_zono, "host_mu_sigma");
  touch_all(host_g_p_x, num_zono, "host_g_p_x");

  const auto outputs = ::roahm::pre_slice::PreSliceImpl(
      grid_size, num_zono, host_x0, host_y0, host_dx, host_dy,
      host_u0v0r0_slice_beta, host_g_u0_x, host_g_u0_y, host_g_v0_x,
      host_g_v0_y, host_g_r0_x, host_g_r0_y, host_block_inzono_list,
      host_rot_angle, host_mu_sigma, host_g_p_x);

  // set output [x0, y0, H1, H2, H4, t]
  plhs[0] = mxCreateNumericMatrix(1, num_zono, mxDOUBLE_CLASS, mxREAL);
  double* x0 = (double*)mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(1, num_zono, mxDOUBLE_CLASS, mxREAL);
  double* y0 = (double*)mxGetData(plhs[1]);

  plhs[2] = mxCreateNumericMatrix(1, num_zono * grid_size * grid_size,
                                  mxDOUBLE_CLASS, mxREAL);
  double* H1 = (double*)mxGetData(plhs[2]);

  plhs[3] = mxCreateNumericMatrix(1, num_zono * grid_size * grid_size,
                                  mxDOUBLE_CLASS, mxREAL);
  double* H2 = (double*)mxGetData(plhs[3]);

  plhs[4] = mxCreateNumericMatrix(1, num_zono * grid_size * grid_size,
                                  mxDOUBLE_CLASS, mxREAL);
  double* H4 = (double*)mxGetData(plhs[4]);

  plhs[5] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  float* gpu_time = (float*)mxGetData(plhs[5]);

  for (int i = 0; i < num_zono; ++i) {
    x0[i] = (outputs.x0_.get())[i];
    y0[i] = (outputs.y0_.get())[i];
    std::cout << "x0[" << i << "] = " << x0[i] << std::endl;
    std::cout << "y0[" << i << "] = " << y0[i] << std::endl;
  }
  for (int i = 0; i < (num_zono * grid_size * grid_size); ++i) {
    H1[i] = (outputs.H1_.get())[i];
    H2[i] = (outputs.H2_.get())[i];
    H4[i] = (outputs.H4_.get())[i];
    if (i == 0) {
      std::cout << "H1[" << i << "] = " << H1[i] << std::endl;
      std::cout << "H2[" << i << "] = " << H2[i] << std::endl;
      std::cout << "H4[" << i << "] = " << H4[i] << std::endl;
    }
  }
  *gpu_time = *outputs.gpu_time_;
}
