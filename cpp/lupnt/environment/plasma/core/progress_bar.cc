/**
 * @file progress_bar.cpp
 * @author Keidai Iiyama
 * @brief Implementation of the progress bar class
 * @version 0.1
 * @date 2025-02-17
 *
 * This file was retrieved from the stack overflow:
 * https://stackoverflow.com/questions/28050669/can-i-report-progress-for-openmp-tasks
 */

#include "lupnt/environment/plasma/core/progress_bar.h"

namespace pecsim {
  /// Start the timer. Throws an exception if timer was already running.
  void Timer::start() {
    if (running) throw std::runtime_error("Timer was already started!");
    running = true;
    start_time = clock::now();
  }

  /// Stop the timer. Throws an exception if timer was already stopped.
  /// Calling this adds to the timer's accumulated time.
  ///@return The accumulated time in seconds.
  double Timer::stop() {
    if (!running) throw std::runtime_error("Timer was already stopped!");

    accumulated_time += lap();
    running = false;

    return accumulated_time;
  }

  /// Returns the timer's accumulated time. Throws an exception if the timer is
  /// running.
  double Timer::accumulated() {
    if (running) throw std::runtime_error("Timer is still running!");
    return accumulated_time;
  }

  /// Returns the time between when the timer was started and the current
  /// moment. Throws an exception if the timer is not running.
  double Timer::lap() {
    if (!running) throw std::runtime_error("Timer was not started!");
    return std::chrono::duration_cast<second>(clock::now() - start_time).count();
  }

  void ProgressBar::start(uint32_t total_work) {
    timer = Timer();
    timer.start();
    this->total_work = total_work;
    next_update = 0;
    call_diff = total_work / 200;
    old_percent = 0;
    work_done = 0;
    clearConsoleLine();
  }

  ///@brief Update the visible progress bar, but only if enough work has been
  /// done.
  ///
  /// Define the global `NOPROGRESS` flag to prevent this from having an
  /// effect. Doing so may speed up the program's execution.
  void ProgressBar::update(uint32_t work_done0) {
// Provide simple way of optimizing out progress updates
#ifdef NOPROGRESS
    return;
#endif

    // Quick return if this isn't the main thread
    if (omp_get_thread_num() != 0) return;

    // Update the amount of work done
    work_done = work_done0;

    // Quick return if insufficient progress has occurred
    if (work_done < next_update) return;

    // Update the next time at which we'll do the expensive update stuff
    next_update += call_diff;

    // Use a uint16_t because using a uint8_t will cause the result to print as
    // a character instead of a number
    uint16_t percent = (uint8_t)(work_done * omp_get_num_threads() * 100 / total_work);

    // Handle overflows
    if (percent > 100) percent = 100;

    // In the case that there has been no update (which should never be the
    // case, actually), skip the expensive screen print
    if (percent == old_percent) return;

    // Update old_percent accordingly
    old_percent = percent;

    // Print an update string which looks like this:
    //   [================================================  ] (96% - 1.0s - 4
    //   threads)
    std::cerr << "\r\033[2K[" << std::string(percent / 2, '=') << std::string(50 - percent / 2, ' ')
              << "] (" << percent << "% - " << std::fixed << std::setprecision(1)
              << timer.lap() / percent * (100 - percent) << "s - " << omp_get_num_threads()
              << " threads)" << std::flush;
  }

}  // namespace pecsim
