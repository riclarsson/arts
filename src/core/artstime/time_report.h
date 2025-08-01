#pragma once

#include <chrono>
#include <source_location>
#include <unordered_map>
#include <vector>

namespace arts {
using time_t          = std::chrono::time_point<std::chrono::system_clock>;
using StartEnd        = std::pair<time_t, time_t>;
using SingleCoreTimer = std::unordered_map<std::string, std::vector<StartEnd>>;
using TimeReport      = std::unordered_map<int, SingleCoreTimer>;

struct profiler {
  std::string name;
  time_t start;

  profiler(std::string&& key);
  profiler(std::source_location loc = std::source_location::current());

  profiler(const profiler&)            = delete;
  profiler(profiler&&)                 = delete;
  profiler& operator=(const profiler&) = delete;
  profiler& operator=(profiler&&)      = delete;

  ~profiler();
};

TimeReport get_report(bool clear = true);
}  // namespace arts

#if ARTS_PROFILING
#define ARTS_TIME_REPORT arts::profiler _arts_prof_var_name_{};
#define ARTS_NAMED_TIME_REPORT(_name_var_) \
  arts::profiler _arts_named__prof_var_name_{_name_var_};
#else
#define ARTS_TIME_REPORT
#define ARTS_NAMED_TIME_REPORT(_name_var_)
#endif
