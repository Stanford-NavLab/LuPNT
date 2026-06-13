#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <thread>

#if defined(_OPENMP) && defined(__has_include)
#  if __has_include(<omp.h>)
#    include <omp.h>
#    define PECSIM_HAS_OPENMP 1
#  endif
#elif defined(_OPENMP)
#  include <omp.h>
#  define PECSIM_HAS_OPENMP 1
#endif

#ifndef PECSIM_HAS_OPENMP
// Fallback shims used when OpenMP headers are unavailable.
#  define omp_get_thread_num() 0
#  define omp_get_num_threads() 1
#endif

namespace pecsim {

  ///@brief Used to time how intervals in code.
  ///
  /// Such as how long it takes a given function to run, or how long I/O has
  /// taken.
  class Timer {
  private:
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, std::ratio<1> > second;

    std::chrono::time_point<clock> start_time;  ///< Last time the timer was started
    double accumulated_time;                    ///< Accumulated running time since creation
    bool running;                               ///< True when the timer is running

  public:
    Timer() {
      accumulated_time = 0;
      running = false;
    }

    /// Start the timer. Throws an exception if timer was already running.
    void start();

    /// Stop the timer. Throws an exception if timer was already stopped.
    /// Calling this adds to the timer's accumulated time.
    ///@return The accumulated time in seconds.
    double stop();

    /// Returns the timer's accumulated time. Throws an exception if the timer is
    /// running.
    double accumulated();

    /// Returns the time between when the timer was started and the current
    /// moment. Throws an exception if the timer is not running.
    double lap();

    /// Stops the timer and resets its accumulated time. No exceptions are thrown
    /// ever.
    void reset() {
      accumulated_time = 0;
      running = false;
    }
  };

  //@brief Manages a console-based progress bar to keep the user entertained.
  ///
  /// Defining the global `NOPROGRESS` will
  /// disable all progress operations, potentially speeding up a program. The look
  /// of the progress bar is shown in ProgressBar.hpp.
  class ProgressBar {
  private:
    uint32_t total_work;   ///< Total work to be accomplished
    uint32_t next_update;  ///< Next point to update the visible progress bar
    uint32_t call_diff;    ///< Interval between updates in work units
    uint32_t work_done;
    uint16_t old_percent;  ///< Old percentage value (aka: should we update the
                           ///< progress bar) TODO: Maybe that we do not need this
    Timer timer;           ///< Used for generating ETA

    /// Clear current line on console so a new progress bar can be written
    void clearConsoleLine() const { std::cerr << "\r\033[2K" << std::flush; }

  public:
    ///@brief Start/reset the progress bar.
    ///@param total_work  The amount of work to be completed, usually specified in
    /// cells.
    void start(uint32_t total_work);

    ///@brief Update the visible progress bar, but only if enough work has been
    /// done.
    ///
    /// Define the global `NOPROGRESS` flag to prevent this from having an
    /// effect. Doing so may speed up the program's execution.
    void update(uint32_t work_done0);

    /// Increment by one the work done and update the progress bar
    ProgressBar& operator++() {
      // Quick return if this isn't the main thread
      if (omp_get_thread_num() != 0) return *this;

      work_done++;
      update(work_done);
      return *this;
    }

    /// Stop the progress bar. Throws an exception if it wasn't started.
    ///@return The number of seconds the progress bar was running.
    double stop() {
      clearConsoleLine();

      timer.stop();
      return timer.accumulated();
    }

    ///@return Return the time the progress bar ran for.
    double time_it_took() { return timer.accumulated(); }

    uint32_t cellsProcessed() const { return work_done; }
  };

}  // namespace pecsim
