#pragma once

#include <map>
#include <set>

#include "lupnt/core/constants.h"

namespace lupnt {
namespace task_scheduling {

enum TaskType {
  User,
  SunPointing,
  Downlink,
};

struct Opportunity {
  int task_id;

  real start_time;
  real end_time;
  real duration;

  real power_consumption;
  real data_generation;

  real reward;
  TaskType type;
};

struct State {
  real time;
  real power;
  real data;
  std::set<int> completed_tasks;
};

struct Action {
  int task_id;
};

// class SatelliteSchedulingModel {
//  private:
//   std::map<int, Opportunity> opportunities_;
//   real gamma_;
//   real horizon_;

//  public:
//   std::vector<Action> GetActions(const State& state) {
//     std::vector<Action> actions;
//     for (const auto& [task_id, opportunity] : opportunities_) {
//       if (state.completed_tasks.count(task_id) == 0) {
//         actions.push_back(Action{task_id});
//       }
//     }
//     return actions;
//   }
// };

// std::tuple<Action, real> SelectActions(const SatelliteSchedulingModel&
// problem,
//                                        State& state, int depth, double gamma)
//                                        {
//   if (depth == 0) {
//     return std::make_tuple(Action{-1}, 0);
//   }
//   Action best_action(-1);
//   real best_value = -std::numeric_limits<double>::infinity();
//   for (const auto& action : GetActions(problem, state)) {
//     real value = Rewar:
//   }
//   return std::make_tuple(best_action, best_value);
// }

}  // namespace task_scheduling

}  // namespace lupnt
