/**
 * File:    simple_clock.h
 *
 * Date:    22.09.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_UTILS_SIMPLE_CLOCK_H_
#define TRIVIS_UTILS_SIMPLE_CLOCK_H_

#include <ctime>
#include <chrono>

namespace trivis::utils {

namespace time = std::chrono;
using Clock = time::high_resolution_clock;
using TimePoint = time::high_resolution_clock::time_point;

class SimpleClock {
 public:

  SimpleClock();

  void Pause();

  void Unpause();

  void Restart();

  [[nodiscard]] double TimeInSeconds() const;

  static TimePoint Now() noexcept;

 private:

  bool is_running_ = true;
  double time_in_seconds_ = 0.0;
  TimePoint start_;

};

inline SimpleClock::SimpleClock() : start_(Clock::now()) {

}

inline void SimpleClock::Pause() {
  auto end = Clock::now();
  time_in_seconds_ += static_cast<double>(time::duration_cast<time::nanoseconds>(end - start_).count()) * 1e-9;
  is_running_ = false;
}

inline void SimpleClock::Unpause() {
  start_ = Clock::now();
  is_running_ = true;
}

inline void SimpleClock::Restart() {
  time_in_seconds_ = 0.0;
  start_ = Clock::now();
  is_running_ = true;
}

[[nodiscard]] inline double SimpleClock::TimeInSeconds() const {
  auto end = Clock::now();
  if (is_running_) {
    return time_in_seconds_ +
        static_cast<double>(time::duration_cast<time::nanoseconds>(end - start_).count()) * 1e-9;
  } else {
    return time_in_seconds_;
  }
}

inline TimePoint SimpleClock::Now() noexcept {
  return time::high_resolution_clock::now();
}

}

#endif //TRIVIS_UTILS_SIMPLE_CLOCK_H_
