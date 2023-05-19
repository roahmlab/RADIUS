#include "mex.h"
#include "prob_integ.hpp"
#include "prob_integ_io.hpp"
/// Compute the constraint value, its gradient, and its Hessian given grid
/// dimensions, block locations that intersect the zonotope, and obstacle PDFs
///
/// Call from MATLAB as:
/// [constr, dconstr, d2constr, t] = prob_integration_frs_flzono(num_zono, ...
///   x0, y0, dx, dy, block_inzono_list, H1, H2, H4, rot_angle, mu_sigma, ...
///   cg_p, g_p_x, grid_size, p)
///
/// @param[in]
/// prhs[0]: (1 x 1 double)
///   the number of zonotopes in this FRS.
/// prhs[1]: (num_zonos x 1 double)
///   [x0_zono1, x0_zono2, ...] the x origin of each zonotope's grid.
///   Each x0_zono# is a 1x1 double.
/// prhs[2]: (num_zonos x 1 double)
///   [y0_zono1, y0_zono2, ...] the y offset/origin of each zonotope's grid.
///   Each y0_zono# is a 1x1 double.
/// prhs[3]: (num_zonos x 1 double)
///   [dx_zono1, dx_zono2, ...] the x grid spacing (column width) of the grid
///   on each zonotope. Each dx_zono# is a 1x1 double.
/// prhs[4]: (num_zonos x 1 double)
///   [dy_zono1, dy_zono2, ...] the y grid spacing (row height) of the grid
///   on each zonotope. Each dy_zono# is a 1x1 double.
/// prhs[5]: (num_zonos x grid_size x grid_size double)
///   [block_inzonolist_zono1, block_inzonolist_zono2, ...] the list of blocks
///   contained in each zonotope, each block outside the zonotope is -1.
///   Each block_inzonolist_zono# is a (grid_size x grid_size) array of vectors,
///   with ordering TODO:
/// prhs[6]: (num_zono * grid_size * grid_size double)
///   [H1_zono1, H1_zono2, ...] TODO: more explanation
///   each H1_zono# contains values of H1 for all patches (grid_size x grid_size
///   double)
/// prhs[7]: (num_zono * grid_size * grid_size double)
///   [H2_zono1, H2_zono2, ...] TODO: more explanation
/// prhs[8]: (num_zono * grid_size * grid_size double)
///   [H4_zono1, H4_zono2, ...] TODO: more explanation
/// prhs[9]: (num_zonos x 1 double)
///   [rot_angle_zono1, rot_angle_zono2, ...] the angles to rotate such that the
///   p-sliceable (trajectory parameter sliceable) generator of zono# is rotated
///   to [x; 0] along its xy dimension. Each rot_angle_zono# is a 1x1 double.
/// prhs[10]: (5 x num_zonos double)
///   [mu_sigma_zono1, mu_sigma_zono2, ...] each mu_sigma_zono# is a 1x5 array,
///   [mu_x, mu_y, sigma_1, sigma_2, sigma_4]
/// prhs[11]: (2 x 1 double)
///   [cp, gp] cp is the center of the parameter dimension, gp is the
///   generator of the parameter dimension, used for slicing
/// prhs[12]: (num_zonos x 1 double)
///   [g_p_x_zono1, g_p_x_zono2, ...] TODO: more explanation
/// prhs[13]: (1 x 1 double)
///   grid_size, the number of grid points in each dimension
/// prhs[14]: (1 x 1 double)
///   p, the trajectory parameter to slice at
///
/// @param[out]
/// plhs[0]: (1 x 1 double)
///   the chance constraint value
/// plhs[1]: (1 x 1 double)
///   the derivative of the chance constraint value w.r.t. the slice parameter
/// plhs[2]: (1 x 1 double)
///   the second derivative of the chance constraint value w.r.t. the slice
///   parameter
/// plhs[3]: (1 x 1 double)
///   the computation time, not including memory allocations, copying to/from
///   device, or freeing memory
///   TODO: this does actually include some copy time, at least.
void mexFunction(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]) { // main function
  //
  // Setup input data
  //
  ::roahm::ProbIntegrationInputs inputs;
  inputs.num_zonos_ = (int)(*(double*)mxGetData(prhs[0]));
  inputs.x0_ = (double*)mxGetData(prhs[1]);
  inputs.y0_ = (double*)mxGetData(prhs[2]);
  inputs.dx_ = (double*)mxGetData(prhs[3]);
  inputs.dy_ = (double*)mxGetData(prhs[4]);
  inputs.block_inzono_list_ =
      (double*)mxGetData(prhs[5]); // this is stupid, MATLAB has to load as
                                   // double?! int* doesn't work???
  inputs.H1_ = (double*)mxGetData(prhs[6]);
  inputs.H2_ = (double*)mxGetData(prhs[7]);
  inputs.H4_ = (double*)mxGetData(prhs[8]);
  inputs.rot_angle_ = (double*)mxGetData(prhs[9]);
  inputs.mu_sigma_ = (double*)mxGetData(prhs[10]);
  inputs.cg_p_ = (double*)mxGetData(prhs[11]);
  inputs.g_p_x_ = (double*)mxGetData(prhs[12]);
  double p = (double)(*(double*)mxGetData(prhs[14]));
  inputs.p_ = (p - inputs.cg_p_[0]) / inputs.cg_p_[1];
  inputs.grid_size_ = (int)(*(double*)mxGetData(prhs[13]));

  // Run the computation
  auto out_vals = ::roahm::prob_integ::ProbIntegration(inputs);

  //
  // Get data pointers to MATLAB output
  //
  plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  double* constr = (double*)mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  double* dconstr = (double*)mxGetData(plhs[1]);

  plhs[2] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  double* d2constr = (double*)mxGetData(plhs[2]);

  plhs[3] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  float* gpu_time = (float*)mxGetData(plhs[3]);

  //
  // Write data into MATLAB output
  //
  *constr = out_vals.constraint_val_;
  *dconstr = out_vals.d_constraint_val_;
  *d2constr = out_vals.d2_constraint_val_;
  *gpu_time = out_vals.computation_time_;
}