/*
 * spinlock.hpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_SLOCK_HPP
#define DISTRIBUTEDNE_SLOCK_HPP

#include <atomic>

class SLock {
  std::atomic_flag locked = ATOMIC_FLAG_INIT;
 public:
  void lock() {
    while (locked.test_and_set(std::memory_order_acquire)) { }
  }
  void unlock() {
    locked.clear(std::memory_order_release);
  }
};

using SpinLocks = std::vector<SLock>;

#endif //DISTRIBUTEDNE_SLOCK_HPP
