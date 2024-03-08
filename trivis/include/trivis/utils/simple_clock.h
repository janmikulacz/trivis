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

using Clock = std::chrono::steady_clock;
using TimePoint = Clock::time_point;

class SimpleClock {
public:

    SimpleClock();

    void Pause();

    void Unpause();

    void Restart();

    [[nodiscard]] double TimeInSeconds() const;

    static TimePoint Now() noexcept;

private:

    bool _is_running = true;
    double _time_in_seconds = 0.0;
    TimePoint _start;

    [[nodiscard]] long ElapsedInNanoSeconds() const;
};

inline SimpleClock::SimpleClock() : _start(Clock::now()) {

}

inline void SimpleClock::Pause() {
    _time_in_seconds += static_cast<double>(ElapsedInNanoSeconds()) * 1e-9;
    _is_running = false;
}

inline void SimpleClock::Unpause() {
    _start = Clock::now();
    _is_running = true;
}

inline void SimpleClock::Restart() {
    _time_in_seconds = 0.0;
    Unpause();
}

[[nodiscard]] inline double SimpleClock::TimeInSeconds() const {
    auto end = Clock::now();
    if (_is_running) {
        return _time_in_seconds + static_cast<double>(ElapsedInNanoSeconds()) * 1e-9;
    } else {
        return _time_in_seconds;
    }
}

inline TimePoint SimpleClock::Now() noexcept {
    return Clock::now();
}

inline long SimpleClock::ElapsedInNanoSeconds() const {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - _start).count();
    }

}

#endif //TRIVIS_UTILS_SIMPLE_CLOCK_H_
