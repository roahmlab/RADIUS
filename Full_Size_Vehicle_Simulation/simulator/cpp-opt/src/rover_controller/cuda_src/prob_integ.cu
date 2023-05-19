//
// Header PROB_FUNCTION_INFO_FRS_FLZONO_H Begin
//

#ifndef PROB_FUNCTION_INFO_FRS_FLZONO_H
#define PROB_FUNCTION_INFO_FRS_FLZONO_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define MAX_ZONO_NUM 1200

__constant__ int dev_num_zono[1];
__constant__ int dev_grid_size[1];
__constant__ double dev_dx[MAX_ZONO_NUM];
__constant__ double dev_dy[MAX_ZONO_NUM];
__constant__ double dev_rot_angle[MAX_ZONO_NUM];
__constant__ double dev_g_p_x[MAX_ZONO_NUM];
__constant__ double dev_p_slice_beta[1];
__constant__ double dev_cg_p[2];

// probablity density function
__device__ inline void f_fun(double& f, double* x, double mu1, double mu2,
                             double s1_1, double s2_1, double s2_2,
                             int zono_idx) {
  double x1 = x[0];
  double x2 = x[1];

  //     double mu1 = mu_sigma[0];
  //     double mu2 = mu_sigma[1];
  //     double s1_1 = mu_sigma[2];
  //     double s2_1 = mu_sigma[3];
  //     double s2_2 = mu_sigma[4];

  double R1_1 = cos(dev_rot_angle[zono_idx]);
  double R1_2 = -sin(dev_rot_angle[zono_idx]);
  double R2_1 = -R1_2;
  double R2_2 = R1_1;

  double t2 = s1_1 * s2_2;
  double t3 = s2_1 * s2_1;
  double t4 = mu1 / 2.0;
  double t5 = mu2 / 2.0;
  double t6 = (R1_1 * x1) / 2.0;
  double t7 = (R1_2 * x1) / 2.0;
  double t8 = (R2_1 * x2) / 2.0;
  double t9 = (R2_2 * x2) / 2.0;
  double t10 = -t3;
  double t11 = -t4;
  double t12 = -t5;
  double t13 = t2 + t10;
  double t15 = t6 + t8 + t11;
  double t16 = t7 + t9 + t12;
  double t14 = 1.0 / t13;
  f = (1.0 / sqrt(t13) *
       exp(-(s1_1 * t14 * t16 - s2_1 * t14 * t15) *
               (-mu2 + R1_2 * x1 + R2_2 * x2) +
           (s2_1 * t14 * t16 - s2_2 * t14 * t15) *
               (-mu1 + R1_1 * x1 + R2_1 * x2))) /
      (3.141592653589793 * 2.0);
}

// first order derivative of probablity density function
__device__ inline void Jacobian_fun(double* J, double* x, double mu1,
                                    double mu2, double s1_1, double s2_1,
                                    double s2_2, int zono_idx) {
  double x1 = x[0];
  double x2 = x[1];

  //     double mu1 = mu_sigma[0];
  //     double mu2 = mu_sigma[1];
  //     double s1_1 = mu_sigma[2];
  //     double s2_1 = mu_sigma[3];
  //     double s2_2 = mu_sigma[4];

  double R1_1 = cos(dev_rot_angle[zono_idx]);
  double R1_2 = -sin(dev_rot_angle[zono_idx]);
  double R2_1 = -R1_2;
  double R2_2 = R1_1;

  double t2 = R1_1 * x1;
  double t3 = R1_2 * x1;
  double t4 = R2_1 * x2;
  double t5 = R2_2 * x2;
  double t6 = s1_1 * s2_2;
  double t7 = s2_1 * s2_1;
  double t8 = 1.0 / 3.141592653589793;
  double t9 = -mu1;
  double t10 = -mu2;
  double t11 = mu1 / 2.0;
  double t12 = mu2 / 2.0;
  double t13 = t2 / 2.0;
  double t14 = t3 / 2.0;
  double t15 = t4 / 2.0;
  double t16 = t5 / 2.0;
  double t17 = -t7;
  double t18 = -t11;
  double t19 = -t12;
  double t21 = t2 + t4 + t9;
  double t22 = t3 + t5 + t10;
  double t20 = t6 + t17;
  double t25 = t13 + t15 + t18;
  double t26 = t14 + t16 + t19;
  double t23 = 1.0 / t20;
  double t24 = 1.0 / sqrt(t20);
  double t27 = s1_1 * t23 * t26;
  double t28 = s2_1 * t23 * t25;
  double t29 = s2_2 * t23 * t25;
  double t30 = s2_1 * t23 * t26;
  double t31 = -t28;
  double t32 = -t30;
  double t33 = t27 + t31;
  double t34 = t29 + t32;
  double t35 = t22 * t33;
  double t36 = t21 * t34;
  double t37 = -t35;
  double t38 = -t36;
  double t39 = t37 + t38;
  double t40 = exp(t39);
  J[0] = t8 * t24 * t40 *
         (R1_1 * t34 + R1_2 * t33 +
          t22 * ((R1_2 * s1_1 * t23) / 2.0 - (R1_1 * s2_1 * t23) / 2.0) +
          t21 * ((R1_1 * s2_2 * t23) / 2.0 - (R1_2 * s2_1 * t23) / 2.0)) *
         (-1.0 / 2.0);
  J[1] = t8 * t24 * t40 *
         (R2_1 * t34 + R2_2 * t33 +
          t22 * ((R2_2 * s1_1 * t23) / 2.0 - (R2_1 * s2_1 * t23) / 2.0) +
          t21 * ((R2_1 * s2_2 * t23) / 2.0 - (R2_2 * s2_1 * t23) / 2.0)) *
         (-1.0 / 2.0);
}

// second order derivative of probablity density function
__device__ inline void Hessian_fun(double* H, double* x, double mu1, double mu2,
                                   double s1, double s2, double s3,
                                   int zono_idx) {
  // this function only returns ddf/dxdx and ddf/dxdy.
  double x1 = x[0];
  double x2 = x[1];

  //     double mu1 = mu_sigma[0];
  //     double mu2 = mu_sigma[1];
  //     double s1 = mu_sigma[2];
  //     double s2 = mu_sigma[3];
  //     double s3 = mu_sigma[4];

  double r1 = cos(dev_rot_angle[zono_idx]);
  double r2 = sin(dev_rot_angle[zono_idx]);
  double r3 = -r2;
  double r4 = r1;

  double t2 = s1 * s3;
  double t3 = r1 * x1;
  double t4 = r2 * x2;
  double t5 = r3 * x1;
  double t6 = r4 * x2;
  double t7 = s2 * s2;
  double t8 = 1.0 / 3.141592653589793;
  double t9 = -mu1;
  double t10 = -mu2;
  double t11 = mu1 / 2.0;
  double t12 = mu2 / 2.0;
  double t13 = -t7;
  double t14 = t3 / 2.0;
  double t15 = t4 / 2.0;
  double t16 = t5 / 2.0;
  double t17 = t6 / 2.0;
  double t18 = -t11;
  double t19 = -t12;
  double t21 = t3 + t4 + t9;
  double t22 = t5 + t6 + t10;
  double t20 = t2 + t13;
  double t35 = t14 + t15 + t18;
  double t36 = t16 + t17 + t19;
  double t23 = 1.0 / t20;
  double t24 = sqrt(t20);
  double t25 = t24;
  double t27 = (r1 * s2 * t23) / 2.0;
  double t28 = (r1 * s3 * t23) / 2.0;
  double t29 = (r2 * s2 * t23) / 2.0;
  double t30 = (r3 * s1 * t23) / 2.0;
  double t31 = (r2 * s3 * t23) / 2.0;
  double t32 = (r3 * s2 * t23) / 2.0;
  double t33 = (r4 * s1 * t23) / 2.0;
  double t34 = (r4 * s2 * t23) / 2.0;
  double t41 = s2 * t23 * t35;
  double t42 = s3 * t23 * t35;
  double t43 = s1 * t23 * t36;
  double t44 = s2 * t23 * t36;
  double t26 = 1.0 / t25;
  double t37 = -t30;
  double t38 = -t32;
  double t39 = -t33;
  double t40 = -t34;
  double t45 = -t43;
  double t46 = -t44;
  double t47 = t27 + t37;
  double t48 = t28 + t38;
  double t49 = t29 + t39;
  double t50 = t31 + t40;
  double t54 = t41 + t45;
  double t55 = t42 + t46;
  double t51 = t21 * t48;
  double t52 = t22 * t47;
  double t56 = r1 * t55;
  double t57 = r3 * t54;
  double t59 = t21 * t55;
  double t60 = t22 * t54;
  double t53 = -t52;
  double t58 = -t57;
  double t61 = -t59;
  double t62 = t60 + t61;
  double t64 = t51 + t53 + t56 + t58;
  double t63 = exp(t62);
  H[0] = (t8 * t26 * t63 * (t64 * t64)) / 2.0 -
         (t8 * t26 * t63 * (r1 * t48 * 2.0 - r3 * t47 * 2.0)) / 2.0;
  H[1] =
      t8 * t26 * t63 * (r2 * t48 + r1 * t50 - r4 * t47 - r3 * t49) *
          (-1.0 / 2.0) +
      (t8 * t26 * t63 * t64 * (r2 * t55 - r4 * t54 + t21 * t50 - t22 * t49)) /
          2.0;
}

__device__ inline void d3f_fun(double* d3f, double* x, double mu1, double mu2,
                               double s1, double s2, double s3, int zono_idx) {
  // this function only returns dddf/dxdxdx and dddf/dxdydx.
  double x1 = x[0];
  double x2 = x[1];

  //     double mu1 = mu_sigma[0];
  //     double mu2 = mu_sigma[1];
  //     double s1 = mu_sigma[2];
  //     double s2 = mu_sigma[3];
  //     double s3 = mu_sigma[4];

  double r1 = cos(dev_rot_angle[zono_idx]);
  double r2 = sin(dev_rot_angle[zono_idx]);
  double r3 = -r2;
  double r4 = r1;

  double t2 = s1 * s3;
  double t3 = r1 * x1;
  double t4 = r2 * x2;
  double t5 = r3 * x1;
  double t6 = r4 * x2;
  double t7 = s2 * s2;
  double t8 = 1.0 / 3.141592653589793;
  double t9 = -mu1;
  double t10 = -mu2;
  double t11 = mu1 / 2.0;
  double t12 = mu2 / 2.0;
  double t13 = -t7;
  double t14 = t3 / 2.0;
  double t15 = t4 / 2.0;
  double t16 = t5 / 2.0;
  double t17 = t6 / 2.0;
  double t18 = -t11;
  double t19 = -t12;
  double t21 = t3 + t4 + t9;
  double t22 = t5 + t6 + t10;
  double t20 = t2 + t13;
  double t35 = t14 + t15 + t18;
  double t36 = t16 + t17 + t19;
  double t23 = 1.0 / t20;
  double t24 = sqrt(t20);
  double t25 = t24;
  double t27 = (r1 * s2 * t23) / 2.0;
  double t28 = (r1 * s3 * t23) / 2.0;
  double t29 = (r2 * s2 * t23) / 2.0;
  double t30 = (r3 * s1 * t23) / 2.0;
  double t31 = (r2 * s3 * t23) / 2.0;
  double t32 = (r3 * s2 * t23) / 2.0;
  double t33 = (r4 * s1 * t23) / 2.0;
  double t34 = (r4 * s2 * t23) / 2.0;
  double t41 = s2 * t23 * t35;
  double t42 = s3 * t23 * t35;
  double t43 = s1 * t23 * t36;
  double t44 = s2 * t23 * t36;
  double t26 = 1.0 / t25;
  double t37 = -t30;
  double t38 = -t32;
  double t39 = -t33;
  double t40 = -t34;
  double t45 = -t43;
  double t46 = -t44;
  double t47 = t27 + t37;
  double t48 = t28 + t38;
  double t49 = t29 + t39;
  double t50 = t31 + t40;
  double t60 = t41 + t45;
  double t61 = t42 + t46;
  double t51 = r1 * t48 * 2.0;
  double t52 = r3 * t47 * 2.0;
  double t54 = t21 * t48;
  double t55 = t21 * t50;
  double t56 = t22 * t47;
  double t57 = t22 * t49;
  double t62 = r1 * t61;
  double t63 = r3 * t60;
  double t64 = r2 * t61;
  double t65 = r4 * t60;
  double t68 = t21 * t61;
  double t69 = t22 * t60;
  double t53 = -t52;
  double t58 = -t56;
  double t59 = -t57;
  double t66 = -t63;
  double t67 = -t65;
  double t70 = -t68;
  double t71 = t51 + t53;
  double t72 = t69 + t70;
  double t74 = t54 + t58 + t62 + t66;
  double t75 = t55 + t59 + t64 + t67;
  double t73 = exp(t72);

  d3f[0] = t8 * t26 * t73 * (t74 * t74 * t74) * (-1.0 / 2.0) +
           t8 * t26 * t71 * t73 * t74 * (3.0 / 2.0);
  d3f[1] = t8 * t26 * t73 * (t74 * t74) * t75 * (-1.0 / 2.0) +
           (t8 * t26 * t71 * t73 * t75) / 2.0 +
           t8 * t26 * t73 * t74 * (r2 * t48 + r1 * t50 - r4 * t47 - r3 * t49);
}

// get the interval bound of the intergral of the probability density function
// over a 2d triangle
__device__ inline void integrator(double* constr, double* xy0, double* H,
                                  double mu1, double mu2, double s1, double s2,
                                  double s3, int zono_idx, int uplo) {
  double f0 = 0;
  f_fun(f0, xy0, mu1, mu2, s1, s2, s3, zono_idx);

  double df0[2] = {0};
  Jacobian_fun(df0, xy0, mu1, mu2, s1, s2, s3, zono_idx);

  double d2f0[2] = {0};
  Hessian_fun(d2f0, xy0, mu1, mu2, s1, s2, s3, zono_idx);

  double d3f0[2] = {0};
  d3f_fun(d3f0, xy0, mu1, mu2, s1, s2, s3, zono_idx);

  double df1 = df0[0];
  double df2 = df0[1];

  double d2f1 = d2f0[0];
  double d2f2 = d2f0[1];

  double d3f1 = d3f0[0];
  double d3f2 = d3f0[1];

  double H1_1 = H[0];
  double H2_1 = H[1];
  double H2_2 = H[2];

  double dx = dev_dx[zono_idx];
  double dy = dev_dy[zono_idx];
  double dxdp = dev_g_p_x[zono_idx] / dev_cg_p[1];

  double t2 = dxdp * dxdp;

  constr[0] = (dx * dy *
               (f0 + (df1 * dx * uplo) / 3.0 + (df2 * dy * uplo) / 3.0 +
                (H1_1 * (dx * dx)) / 1.2E+1 + (H2_2 * (dy * dy)) / 1.2E+1 +
                (H2_1 * dx * dy) / 1.2E+1)) /
              2.0;
  constr[1] = (dx * dy *
               (df1 * dxdp +
                dxdp * ((d2f1 * dx * uplo) / 3.0 + (d2f2 * dy * uplo) / 3.0))) /
              2.0;
  constr[2] = (dx * dy *
               (d2f1 * t2 +
                t2 * ((d3f1 * dx * uplo) / 3.0 + (d3f2 * dy * uplo) / 3.0))) /
              2.0;

  //     constr[0] = df0[1];
  //     constr[1] = df0[2];
}

#endif

//
// Header PROB_FUNCTION_INFO_FRS_FLZONO_H End
//

#include <unistd.h>

#include <cstdio>
#include <cstdlib>

#include "prob_integ.hpp"

// #define MAX_PATCHES_NUM 1024    // This is maximum patch number per zono!
// Otherwise need to change the sum_over kernel
#define BLOCK_SIZE1 32   // block size for integration kernel
#define BLOCK_SIZE2 1024 // block size for sum_over kernel

__global__ void __launch_bounds__(BLOCK_SIZE1, 1)
    eval_patch_kernel(double* x0, double* y0, double* block_inzono_list,
                      double* H1, double* H2, double* H4, double* mu_sigma,
                      double* constr, double* dconstr, double* d2constr) {
  int patch_id =
      (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;

  double H[3] = {H1[patch_id], H2[patch_id], H4[patch_id]};

  // throw xy0, mu and sigma into shared memory
  __shared__ volatile double shared_mu_sigma[5];
  __shared__ volatile double shared_xy0[2];

  if (threadIdx.x < 5) {
    shared_mu_sigma[threadIdx.x] = mu_sigma[blockIdx.y * 5 + threadIdx.x];
    if (threadIdx.x == 0) {
      shared_xy0[0] =
          x0[blockIdx.y] + dev_g_p_x[blockIdx.y] * dev_p_slice_beta[0];
      shared_xy0[1] = y0[blockIdx.y];
    }
  }
  __syncthreads();

  // actual computation
  double x[2]; // vertex of the right triangle in a simplex
  int pid_inzono = (int)block_inzono_list[patch_id];
  if (pid_inzono < 0)
    return;
  int xid = pid_inzono % dev_grid_size[0];
  int yid = pid_inzono / dev_grid_size[0];

  // lower right triangle
  double constr1[3];
  x[0] = shared_xy0[0] + xid * dev_dx[blockIdx.y];
  x[1] = shared_xy0[1] + yid * dev_dy[blockIdx.y];
  integrator(constr1, x, H, shared_mu_sigma[0], shared_mu_sigma[1],
             shared_mu_sigma[2], shared_mu_sigma[3], shared_mu_sigma[4],
             blockIdx.y, 1);

  //     double bla[2];
  //     Hessian_fun(bla, x, shared_mu_sigma[0], shared_mu_sigma[1],
  //     shared_mu_sigma[2],
  //                     shared_mu_sigma[3], shared_mu_sigma[4], zono_id);

  // upper right triangle
  double constr2[3];
  x[0] += dev_dx[blockIdx.y];
  x[1] += dev_dy[blockIdx.y];
  integrator(constr2, x, H, shared_mu_sigma[0], shared_mu_sigma[1],
             shared_mu_sigma[2], shared_mu_sigma[3], shared_mu_sigma[4],
             blockIdx.y, -1);

  constr[patch_id] = constr1[0] + constr2[0];
  dconstr[patch_id] = constr1[1] + constr2[1];
  d2constr[patch_id] = constr1[2] + constr2[2];

  //     constr[patch_id] = H[0];
  //     dconstr[patch_id] = H[1];
  //     d2constr[patch_id] = H[2];
}

__global__ void sum_over(double* input, double* output, int len) {
  int tid = threadIdx.x;
  int bid = blockIdx.x;
  int idx = bid * BLOCK_SIZE2 + tid;

  __shared__ volatile double shared_sum[BLOCK_SIZE2];

  if (idx < len)
    shared_sum[tid] = input[idx];
  else
    shared_sum[tid] = 0;
  __syncthreads();

  if (BLOCK_SIZE2 >= 1024) {
    if (tid < 512) {
      shared_sum[tid] += shared_sum[tid + 512];
    }
    __syncthreads();
  }
  if (BLOCK_SIZE2 >= 512) {
    if (tid < 256) {
      shared_sum[tid] += shared_sum[tid + 256];
    }
    __syncthreads();
  }
  if (BLOCK_SIZE2 >= 256) {
    if (tid < 128) {
      shared_sum[tid] += shared_sum[tid + 128];
    }
    __syncthreads();
  }
  if (BLOCK_SIZE2 >= 128) {
    if (tid < 64) {
      shared_sum[tid] += shared_sum[tid + 64];
    }
    __syncthreads();
  }
  if (tid < 32) {
    if (BLOCK_SIZE2 >= 64)
      shared_sum[tid] += shared_sum[tid + 32];
    if (BLOCK_SIZE2 >= 32)
      shared_sum[tid] += shared_sum[tid + 16];
    if (BLOCK_SIZE2 >= 16)
      shared_sum[tid] += shared_sum[tid + 8];
    if (BLOCK_SIZE2 >= 8)
      shared_sum[tid] += shared_sum[tid + 4];
    if (BLOCK_SIZE2 >= 4)
      shared_sum[tid] += shared_sum[tid + 2];
    if (BLOCK_SIZE2 >= 2)
      shared_sum[tid] += shared_sum[tid + 1];
  }

  if (tid == 0) {
    output[bid] = shared_sum[0];
  }
}

// helper functions
inline void initializeDeviceArray(double** dev_arr, int size) {
  cudaMalloc((void**)dev_arr, size * sizeof(double));
  cudaMemset(*dev_arr, 0, size * sizeof(double));
}

// inline void initializeDeviceArray(double** dev_arr, double* host_arr, int
// size) {
//     cudaMalloc((void**)dev_arr, size * sizeof(double));
//     cudaMemcpy(*dev_arr, host_arr, size * sizeof(double),
//     cudaMemcpyHostToDevice);
// }
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

/// Performs the integration of the probability function over the zonotopes.
///
/// @param[in] inputs The inputs to the integration, see the member variable
/// documentation for more information.
///
/// @param[out] outputs The outputs of the integration, see the member variable
/// documentation for more information.
namespace roahm::prob_integ {
::roahm::ProbIntegrationOutputs
ProbIntegration(const ::roahm::ProbIntegrationInputs& inputs) {
  const int num_zono = inputs.num_zonos_;
  const double* host_x0 = inputs.x0_;
  const double* host_y0 = inputs.y0_;
  const double* host_dx = inputs.dx_;
  const double* host_dy = inputs.dy_;
  const double* host_block_inzono_list = inputs.block_inzono_list_;
  const double* host_H1 = inputs.H1_;
  const double* host_H2 = inputs.H2_;
  const double* host_H4 = inputs.H4_;
  const double* host_rot_angle = inputs.rot_angle_;
  const double* host_mu_sigma = inputs.mu_sigma_;
  const double* cg_p = inputs.cg_p_;
  const double* host_g_p_x = inputs.g_p_x_;
  const int grid_size = inputs.grid_size_;
  const double p = inputs.p_;
  const int total_num_patch = num_zono * grid_size * grid_size;

  cudaMemcpyToSymbol(dev_grid_size, &grid_size, sizeof(int));
  cudaMemcpyToSymbol(dev_num_zono, &num_zono, sizeof(int));

  double* dev_x0 = nullptr;
  initializeDeviceArray(&dev_x0, host_x0, num_zono);

  double* dev_y0 = nullptr;
  initializeDeviceArray(&dev_y0, host_y0, num_zono);

  cudaMemcpyToSymbol(dev_dx, host_dx, num_zono * sizeof(double));

  cudaMemcpyToSymbol(dev_dy, host_dy, num_zono * sizeof(double));

  double* dev_block_inzono_list = nullptr;
  initializeDeviceArray(&dev_block_inzono_list, host_block_inzono_list,
                        num_zono * grid_size * grid_size);

  double* dev_H1 = nullptr;
  initializeDeviceArray(&dev_H1, host_H1, total_num_patch);

  double* dev_H2 = nullptr;
  initializeDeviceArray(&dev_H2, host_H2, total_num_patch);

  double* dev_H4 = nullptr;
  initializeDeviceArray(&dev_H4, host_H4, total_num_patch);

  cudaMemcpyToSymbol(dev_rot_angle, host_rot_angle, num_zono * sizeof(double));

  double* dev_mu_sigma = nullptr;
  initializeDeviceArray(&dev_mu_sigma, host_mu_sigma, 5 * num_zono);

  cudaMemcpyToSymbol(dev_cg_p, cg_p, 2 * sizeof(double));

  cudaMemcpyToSymbol(dev_g_p_x, host_g_p_x, num_zono * sizeof(double));

  cudaMemcpyToSymbol(dev_p_slice_beta, &p, sizeof(double));

  // Setup output structure
  ::roahm::ProbIntegrationOutputs outputs{};

  // Point variables to the output structur member variables, not cleanest
  // but prevents needing to change anything below.
  double* constr = &outputs.constraint_val_;
  double* dconstr = &outputs.d_constraint_val_;
  double* d2constr = &outputs.d2_constraint_val_;
  float* gpu_time = &outputs.computation_time_;

  cudaEvent_t start, stop;

  // define output in device
  double* dev_constr_1 = nullptr;
  initializeDeviceArray(&dev_constr_1, total_num_patch);
  double* dev_dconstr_1 = nullptr;
  initializeDeviceArray(&dev_dconstr_1, total_num_patch);
  double* dev_d2constr_1 = nullptr;
  initializeDeviceArray(&dev_d2constr_1, total_num_patch);

  int sum_length = total_num_patch;
  int num_sum_block = divideup(sum_length, BLOCK_SIZE2);
  double* dev_constr_2 = nullptr;
  initializeDeviceArray(&dev_constr_2, num_sum_block);
  double* dev_dconstr_2 = nullptr;
  initializeDeviceArray(&dev_dconstr_2, num_sum_block);
  double* dev_d2constr_2 = nullptr;
  initializeDeviceArray(&dev_d2constr_2, num_sum_block);

  //////  MAIN STARTS HERE //////
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  dim3 dimGrid(divideup(grid_size * grid_size, BLOCK_SIZE1), num_zono);
  eval_patch_kernel<<<dimGrid, BLOCK_SIZE1>>>(
      dev_x0, dev_y0, dev_block_inzono_list, dev_H1, dev_H2, dev_H4,
      dev_mu_sigma, dev_constr_1, dev_dconstr_1, dev_d2constr_1);
  cudaThreadSynchronize();
  cudaDeviceSynchronize();

  // sum up all integrations
  bool final_summation_in_1 = true;
  while (sum_length > 1) {
    sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_constr_1, dev_constr_2,
                                             sum_length);
    sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_dconstr_1, dev_dconstr_2,
                                             sum_length);
    sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_d2constr_1, dev_d2constr_2,
                                             sum_length);
    final_summation_in_1 = false;

    sum_length = num_sum_block;
    num_sum_block = divideup(sum_length, BLOCK_SIZE2);
    if (sum_length > 1) {
      sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_constr_2, dev_constr_1,
                                               sum_length);
      sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_dconstr_2, dev_dconstr_1,
                                               sum_length);
      sum_over<<<num_sum_block, BLOCK_SIZE2>>>(dev_d2constr_2, dev_d2constr_1,
                                               sum_length);
      sum_length = num_sum_block;
      final_summation_in_1 = true;
    }
  }

  // copy output
  if (final_summation_in_1) {
    cudaMemcpy(constr, dev_constr_1, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dconstr, dev_dconstr_1, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(d2constr, dev_d2constr_1, sizeof(double),
               cudaMemcpyDeviceToHost);
  } else {
    cudaMemcpy(constr, dev_constr_2, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(dconstr, dev_dconstr_2, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(d2constr, dev_d2constr_2, sizeof(double),
               cudaMemcpyDeviceToHost);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(gpu_time, start, stop);
  // printf("Time spent: %.5f ms\n", *gpu_time);
  //////  MAIN ENDS HERE //////

  // free memory
  cudaFree(dev_x0);
  cudaFree(dev_y0);
  cudaFree(dev_block_inzono_list);
  cudaFree(dev_mu_sigma);
  cudaFree(dev_H1);
  cudaFree(dev_H2);
  cudaFree(dev_H4);
  cudaFree(dev_constr_1);
  cudaFree(dev_dconstr_1);
  cudaFree(dev_d2constr_1);
  cudaFree(dev_constr_2);
  cudaFree(dev_dconstr_2);
  cudaFree(dev_d2constr_2);
  return outputs;
}
} // namespace roahm::prob_integ

#undef BLOCK_SIZE1
#undef BLOCK_SIZE2