#ifndef PAR_HELPERS_H
#define PAR_HELPERS_H

#include <utility> // for std::as_const
#include <csignal>

#include <Eigen/Core>
#include "tatami/utils/parallelize.hpp"
#include "tatami/base/Matrix.hpp"

namespace {
volatile std::sig_atomic_t GOT_SIGTERM = 0;
inline void sigterm_handle(int signum) {
  if (signum == SIGTERM) {
    GOT_SIGTERM = 1;
  }
}
} // namespace

inline bool check_interrupt() { return GOT_SIGTERM == 1; }

template <class Fn> inline void run_par(Fn &f, tatami::NumericMatrix::index_type tasks, int workers) {
  const auto nt = Eigen::nbThreads();
  Eigen::setNbThreads(1);
  const auto old_hdl = std::signal(SIGINT, sigterm_handle); // install our custom threadsafe signal handler
  tatami::parallelize([&f = std::as_const(f)](const auto thread, const auto start, const auto length) -> void { f(start, length); }, tasks, workers);
  std::signal(SIGINT, old_hdl); // replace w/ existing signal handler
  Eigen::setNbThreads(nt);
}

#endif