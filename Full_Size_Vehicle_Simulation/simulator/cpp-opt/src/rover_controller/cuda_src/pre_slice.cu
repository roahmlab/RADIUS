#include <cstdio>
#include <cstdlib>

#include "cuda_interval_lib.hpp"
#include "pre_slice.hpp"

using CudaInterval = interval_gpu<double>;
#define MAX_ZONO_NUM 1000

#define BLOCK_SIZE1 512 // block size for slicing kernel
#define BLOCK_SIZE2 32  // block size for hessian kernel

__constant__ double dev_u0v0r0_slice_beta[3];
__constant__ int dev_num_zono[1];
__constant__ int dev_grid_size[1];
__constant__ double dev_dx[MAX_ZONO_NUM];
__constant__ double dev_dy[MAX_ZONO_NUM];
__constant__ double dev_rot_angle[MAX_ZONO_NUM];
__constant__ double dev_g_p_x[MAX_ZONO_NUM];

// __device__ inline void Hessian_fun(CudaInterval* H, Interval* x, double*
// mu_sigma, int zono_id) {
__device__ inline void Hessian_fun(CudaInterval* H, CudaInterval* x, double mu1,
                                   double mu2, double s1_1, double s2_1,
                                   double s2_2, int zono_id) {
  CudaInterval x1 = x[0];
  CudaInterval x2 = x[1];

  //     double mu1 = mu_sigma[0];
  //     double mu2 = mu_sigma[1];
  //     double s1_1 = mu_sigma[2];
  //     double s2_1 = mu_sigma[3];
  //     double s2_2 = mu_sigma[4];

  double R1_1 = cos(dev_rot_angle[zono_id]);
  double R1_2 = -sin(dev_rot_angle[zono_id]);
  double R2_1 = -R1_2;
  double R2_2 = R1_1;

  CudaInterval t2 = R1_1 * x1;
  CudaInterval t3 = R1_2 * x1;
  CudaInterval t4 = R2_1 * x2;
  CudaInterval t5 = R2_2 * x2;
  double t6 = s1_1 * s2_2;
  double t7 = s2_1 * s2_1;
  double t8 = 1.0 / 3.141592653589793;
  double t9 = -mu1;
  double t10 = -mu2;
  double t11 = mu1 / 2.0;
  double t12 = mu2 / 2.0;
  CudaInterval t13 = t2 / 2.0;
  CudaInterval t14 = t3 / 2.0;
  CudaInterval t15 = t4 / 2.0;
  CudaInterval t16 = t5 / 2.0;
  CudaInterval t17 = -t7;
  CudaInterval t18 = -t11;
  CudaInterval t19 = -t12;
  CudaInterval t21 = t2 + t4 + t9;
  CudaInterval t22 = t3 + t5 + t10;
  CudaInterval t20 = t6 + t17;
  CudaInterval t35 = t13 + t15 + t18;
  CudaInterval t36 = t14 + t16 + t19;
  CudaInterval t23 = 1.0 / t20;
  CudaInterval t24 = sqrt(t20);
  CudaInterval t25 = t24;
  CudaInterval t27 = (R1_2 * s1_1 * t23) / 2.0;
  CudaInterval t28 = (R1_1 * s2_1 * t23) / 2.0;
  CudaInterval t29 = (R1_1 * s2_2 * t23) / 2.0;
  CudaInterval t30 = (R1_2 * s2_1 * t23) / 2.0;
  CudaInterval t31 = (R2_2 * s1_1 * t23) / 2.0;
  CudaInterval t32 = (R2_1 * s2_1 * t23) / 2.0;
  CudaInterval t33 = (R2_1 * s2_2 * t23) / 2.0;
  CudaInterval t34 = (R2_2 * s2_1 * t23) / 2.0;
  CudaInterval t41 = s1_1 * t23 * t36;
  CudaInterval t42 = s2_1 * t23 * t35;
  CudaInterval t43 = s2_2 * t23 * t35;
  CudaInterval t44 = s2_1 * t23 * t36;
  CudaInterval t26 = 1.0 / t25;
  CudaInterval t37 = -t28;
  CudaInterval t38 = -t30;
  CudaInterval t39 = -t32;
  CudaInterval t40 = -t34;
  CudaInterval t45 = -t42;
  CudaInterval t46 = -t44;
  CudaInterval t47 = t27 + t37;
  CudaInterval t48 = t29 + t38;
  CudaInterval t49 = t31 + t39;
  CudaInterval t50 = t33 + t40;
  CudaInterval t59 = t41 + t45;
  CudaInterval t60 = t43 + t46;
  CudaInterval t51 = R2_2 * t47;
  CudaInterval t52 = R1_2 * t49;
  CudaInterval t53 = R2_1 * t48;
  CudaInterval t54 = R1_1 * t50;
  CudaInterval t55 = t22 * t47;
  CudaInterval t56 = t21 * t48;
  CudaInterval t57 = t22 * t49;
  CudaInterval t58 = t21 * t50;
  CudaInterval t61 = R1_2 * t59;
  CudaInterval t62 = R1_1 * t60;
  CudaInterval t63 = R2_2 * t59;
  CudaInterval t64 = R2_1 * t60;
  CudaInterval t65 = t22 * t59;
  CudaInterval t66 = t21 * t60;
  CudaInterval t67 = -t65;
  CudaInterval t68 = -t66;
  CudaInterval t69 = t51 + t52 + t53 + t54;
  CudaInterval t72 = t55 + t56 + t61 + t62;
  CudaInterval t73 = t57 + t58 + t63 + t64;
  CudaInterval t70 = t67 + t68;
  CudaInterval t71 = exp(t70);
  CudaInterval t74 = (t8 * t26 * t69 * t71) / 2.0;
  CudaInterval t76 = (t8 * t26 * t71 * t72 * t73) / 2.0;
  CudaInterval t75 = -t74;
  CudaInterval t77 = t75 + t76;
  H[0] = t8 * t26 * t71 * (R1_1 * t48 * 2.0 + R1_2 * t47 * 2.0) * (-1.0 / 2.0) +
         (t8 * t26 * t71 * (t72 * t72)) / 2.0;
  H[1] = t77;
  H[2] = t8 * t26 * t71 * (R2_1 * t50 * 2.0 + R2_2 * t49 * 2.0) * (-1.0 / 2.0) +
         (t8 * t26 * t71 * (t73 * t73)) / 2.0;
}

__global__ void __launch_bounds__(BLOCK_SIZE2, 1)
    hessian_kernel(double* x0, double* y0, double* mu_sigma, double* H1,
                   double* H2, double* H4, double* block_inzono_list) {
  int zono_id = blockIdx.y;
  int patch_id = (zono_id * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;

  // get mu, sigma, and xy0 for zonotope #zono_id
  __shared__ volatile double shared_mu_sigma[5];
  __shared__ volatile double shared_xy0[2];
  if (threadIdx.x < 5) {
    shared_mu_sigma[threadIdx.x] = mu_sigma[zono_id * 5 + threadIdx.x];
    if (threadIdx.x == 0) {
      shared_xy0[0] = x0[zono_id];
      shared_xy0[1] = y0[zono_id];
    }
  }
  __syncthreads();

  // actual computation
  CudaInterval x[2];
  CudaInterval H[3];
  int pid_inzono = (int)block_inzono_list[patch_id];
  if (pid_inzono < 0)
    return;
  int xid = pid_inzono % dev_grid_size[0];
  int yid = pid_inzono / dev_grid_size[0];
  x[1] = CudaInterval(shared_xy0[1] + yid * dev_dy[zono_id],
                      shared_xy0[1] + (yid + 1) * dev_dy[zono_id]);

  double xlo = shared_xy0[0] + xid * dev_dx[zono_id] - dev_g_p_x[zono_id];
  double dxx = (2 * dev_g_p_x[zono_id] + dev_dx[zono_id]) / 12;
  double H_ub[3] = {nanf(""), nanf(""), nanf("")};
  for (int i = 0; i < 12; i++) {
    x[0] = CudaInterval(xlo, xlo + dxx);
    xlo += dxx;
    Hessian_fun(H, x, shared_mu_sigma[0], shared_mu_sigma[1],
                shared_mu_sigma[2], shared_mu_sigma[3], shared_mu_sigma[4],
                zono_id);
    H_ub[0] = fmax(H_ub[0], H[0].upper());
    H_ub[1] = fmax(H_ub[1], H[1].upper());
    H_ub[2] = fmax(H_ub[2], H[2].upper());
  }

  H1[patch_id] = H_ub[0];
  H2[patch_id] = H_ub[1];
  H4[patch_id] = H_ub[2];
}

__global__ void __launch_bounds__(BLOCK_SIZE1, 1)
    slice_u0v0r0_kernel(double* x0, double* y0, double* g_u0_x, double* g_u0_y,
                        double* g_v0_x, double* g_v0_y, double* g_r0_x,
                        double* g_r0_y) {
  int zono_id = threadIdx.x + blockDim.x * blockIdx.x;
  if (zono_id >= dev_num_zono[0])
    return;

  x0[zono_id] = x0[zono_id] + dev_u0v0r0_slice_beta[0] * g_u0_x[zono_id] +
                dev_u0v0r0_slice_beta[1] * g_v0_x[zono_id] +
                dev_u0v0r0_slice_beta[2] * g_r0_x[zono_id];
  y0[zono_id] = y0[zono_id] + dev_u0v0r0_slice_beta[0] * g_u0_y[zono_id] +
                dev_u0v0r0_slice_beta[1] * g_v0_y[zono_id] +
                dev_u0v0r0_slice_beta[2] * g_r0_y[zono_id];
}

inline void initializeDeviceArray(double** dev_arr, int size) {
  cudaMalloc((void**)dev_arr, size * sizeof(double));
  cudaMemset(*dev_arr, 0, size * sizeof(double));
}

inline void initializeDeviceArray(double** dev_arr, const double* host_arr,
                                  int size) {
  cudaMalloc((void**)dev_arr, size * sizeof(double));
  cudaMemcpy(*dev_arr, host_arr, size * sizeof(double), cudaMemcpyHostToDevice);
}

inline int divideup(int a, int b) {
  int c = a / b;
  if (c * b == a)
    return c;
  return c + 1;
}

namespace roahm::pre_slice {
PreSliceOutputs
PreSliceImpl(const int grid_size, const int num_zono, const double* host_x0,
             const double* host_y0, const double* host_dx,
             const double* host_dy, const double* host_u0v0r0_slice_beta,
             const double* host_g_u0_x, const double* host_g_u0_y,
             const double* host_g_v0_x, const double* host_g_v0_y,
             const double* host_g_r0_x, const double* host_g_r0_y,
             const double* host_block_inzono_list, const double* host_rot_angle,
             const double* host_mu_sigma, const double* host_g_p_x) {
  cudaMemcpyToSymbol(dev_grid_size, &grid_size, sizeof(int));
  cudaMemcpyToSymbol(dev_num_zono, &num_zono, sizeof(int));

  double* dev_x0 = nullptr;
  initializeDeviceArray(&dev_x0, host_x0, num_zono);

  double* dev_y0 = nullptr;
  initializeDeviceArray(&dev_y0, host_y0, num_zono);

  cudaMemcpyToSymbol(dev_dx, host_dx, num_zono * sizeof(double));
  cudaMemcpyToSymbol(dev_dy, host_dy, num_zono * sizeof(double));
  cudaMemcpyToSymbol(dev_u0v0r0_slice_beta, host_u0v0r0_slice_beta,
                     3 * sizeof(double));

  double* dev_g_u0_x = nullptr;
  initializeDeviceArray(&dev_g_u0_x, host_g_u0_x, num_zono);

  double* dev_g_u0_y = nullptr;
  initializeDeviceArray(&dev_g_u0_y, host_g_u0_y, num_zono);

  double* dev_g_v0_x = nullptr;
  initializeDeviceArray(&dev_g_v0_x, host_g_v0_x, num_zono);

  double* dev_g_v0_y = nullptr;
  initializeDeviceArray(&dev_g_v0_y, host_g_v0_y, num_zono);

  double* dev_g_r0_x = nullptr;
  initializeDeviceArray(&dev_g_r0_x, host_g_r0_x, num_zono);

  double* dev_g_r0_y = nullptr;
  initializeDeviceArray(&dev_g_r0_y, host_g_r0_y, num_zono);

  double* dev_block_inzono_list = nullptr;
  initializeDeviceArray(&dev_block_inzono_list, host_block_inzono_list,
                        num_zono * grid_size * grid_size);

  cudaMemcpyToSymbol(dev_rot_angle, host_rot_angle, num_zono * sizeof(double));

  double* dev_mu_sigma = nullptr;
  initializeDeviceArray(&dev_mu_sigma, host_mu_sigma, 5 * num_zono);

  cudaMemcpyToSymbol(dev_g_p_x, host_g_p_x, num_zono * sizeof(double));

  // TODO
  // set output [x0, y0, H1, H2, H4, t]
  PreSliceOutputs outputs{
      std::shared_ptr<double>{new double[num_zono]},
      std::shared_ptr<double>{new double[num_zono]},
      std::shared_ptr<double>{new double[num_zono * grid_size * grid_size]},
      std::shared_ptr<double>{new double[num_zono * grid_size * grid_size]},
      std::shared_ptr<double>{new double[num_zono * grid_size * grid_size]},
      std::shared_ptr<float>{new float[1]},
  };
  double* x0 = outputs.x0_.get();
  double* y0 = outputs.y0_.get();
  double* H1 = outputs.H1_.get();
  double* H2 = outputs.H2_.get();
  double* H4 = outputs.H4_.get();
  float* gpu_time = outputs.gpu_time_.get();

  cudaEvent_t start, stop;
  // define output in device
  double* dev_H1 = nullptr;
  initializeDeviceArray(&dev_H1, num_zono * grid_size * grid_size);
  double* dev_H2 = nullptr;
  initializeDeviceArray(&dev_H2, num_zono * grid_size * grid_size);
  double* dev_H4 = nullptr;
  initializeDeviceArray(&dev_H4, num_zono * grid_size * grid_size);

  //////  MAIN STARTS HERE //////
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  slice_u0v0r0_kernel<<<divideup(num_zono, BLOCK_SIZE1), BLOCK_SIZE1>>>(
      dev_x0, dev_y0, dev_g_u0_x, dev_g_u0_y, dev_g_v0_x, dev_g_v0_y,
      dev_g_r0_x, dev_g_r0_y);

  dim3 dimGrid(divideup(grid_size * grid_size, BLOCK_SIZE2), num_zono);
  hessian_kernel<<<dimGrid, BLOCK_SIZE2>>>(dev_x0, dev_y0, dev_mu_sigma, dev_H1,
                                           dev_H2, dev_H4,
                                           dev_block_inzono_list);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(gpu_time, start, stop);
  // printf("Time spent: %.5f ms\n", *gpu_time);
  //////  MAIN ENDS HERE //////

  // copy output
  cudaMemcpy(x0, dev_x0, num_zono * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(y0, dev_y0, num_zono * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(H1, dev_H1, num_zono * grid_size * grid_size * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(H2, dev_H2, num_zono * grid_size * grid_size * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(H4, dev_H4, num_zono * grid_size * grid_size * sizeof(double),
             cudaMemcpyDeviceToHost);

  // free memory
  cudaFree(dev_x0);
  cudaFree(dev_y0);
  cudaFree(dev_g_u0_x);
  cudaFree(dev_g_u0_y);
  cudaFree(dev_g_v0_x);
  cudaFree(dev_g_v0_y);
  cudaFree(dev_g_r0_x);
  cudaFree(dev_g_r0_y);
  cudaFree(dev_block_inzono_list);
  cudaFree(dev_mu_sigma);
  cudaFree(dev_H1);
  cudaFree(dev_H2);
  cudaFree(dev_H4);
  return outputs;
}

PreSliceOutputs PreSlice(const ::roahm::CudaInfo& cuda_info,
                         const std::array<double, 3>& u0v0r0_slice_beta_in,
                         const double* const mu_sigma) {
  return PreSliceImpl(
      cuda_info.grid_size_, cuda_info.num_zono_, cuda_info.grid_x0_.data(),
      cuda_info.grid_y0_.data(), cuda_info.grid_dx_.data(),
      cuda_info.grid_dy_.data(), u0v0r0_slice_beta_in.data(),
      cuda_info.g_u0_x_.data(), cuda_info.g_u0_y_.data(),
      cuda_info.g_v0_x_.data(), cuda_info.g_v0_y_.data(),
      cuda_info.g_r0_x_.data(), cuda_info.g_r0_y_.data(),
      cuda_info.block_inzono_list_.data(), cuda_info.rot_angle_.data(),
      mu_sigma, cuda_info.g_p_x_.data());
}

} // namespace roahm::pre_slice

#undef MAX_ZONO_NUM
#undef BLOCK_SIZE1
#undef BLOCK_SIZE2
