// /**
//  * @file dynamics_params.cc
//  * @author your name (you@domain.com)
//  * @brief
//  * @version 0.1
//  * @date 2024-12-28
//  *
//  * @copyright Copyright (c) 2024
//  *
//  */

// #include "lupnt/dynamics/dynamics_params.h"

// namespace lupnt {

//   const std::map<ParamsEstOption, std::string> params_est_option2string
//       = {{ParamsEstOption::TrueFixed, "TrueFixed"},
//          {ParamsEstOption::Estimated, "Estimated"},
//          {ParamsEstOption::Consider, "Consider"}};

//   const std::map<std::string, ParamsEstOption> string2params_est_option
//       = {{"TrueFixed", ParamsEstOption::TrueFixed},
//          {"Estimated", ParamsEstOption::Estimated},
//          {"Consider", ParamsEstOption::Consider}};

//   /**
//    * @brief
//    *
//    * @param x0  state + parameters (est or consider)
//    * @param t0  initial time
//    * @param tf  final time
//    * @param stm  state transition matrix
//    * @return VecX
//    */
//   VecX DynamicsWithParams::Propagate(const VecX& x0, Real t0, Real tf, MatXd* stm) {
//     // state size
//     int n = x0.size();
//     int n_param_est = GetEstParamSize();
//     int n_param_consider = GetConsiderParamSize();
//     int n_param = n_param_est + n_param_consider;
//     int state_size = n - n_param;

//     // To compute the jacobians, x_with_params need to capture all estimated and considered
//     // parameters
//     auto func = [=, this](const VecX& x_with_est_params) {
//       VecX x = x_with_est_params.head(state_size);

//       if (use_params_) {
//         // extract the parameters from x_with_est_params
//         dynamics_->SetParams(params_.GetDynamicsParam());
//       }
//       VecX xf = dynamics_->Propagate(x, t0, tf, nullptr);  // Do not compute stm here

//       if (use_params_) {
//         DynamicsParam new_params
//             = dynamics_->GetParams();   // Get the updated parameters (often its constant)
//         UpdateParamValues(new_params);  // Update the parameters in the class
//         VecX new_params_vec
//             = GetEstConsiderParamVector(new_params);  // Get the updated parameters in vector
//             form
//         xf.conservativeResize(
//             n + new_params_vec.size());  // Resize the state vector to include the parameters
//         xf.tail(new_params_vec.size()) = new_params_vec;
//       }
//       return xf;
//     };

//     VecX x0_tmp = x0.cast<double>();
//     VecX xf;
//     if (stm != nullptr) {
//       *stm = jacobian(func, wrt(x0_tmp), at(x0_tmp), xf);
//     } else {
//       xf = func(x0_tmp);
//     }

//     return xf;
//   };

// }  // namespace lupnt
