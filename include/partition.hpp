/*
 * partition.hpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_PARTITION_HPP
#define DISTRIBUTEDNE_PARTITION_HPP

#include <assert.h>
#include <string.h>

#include "type.hpp"

#if !defined(MAX_NUM_PARTS)
#define MAX_NUM_PARTS 1024
#endif

const int P_BITS_SIZE = MAX_NUM_PARTS/8/sizeof(uint64_t);

struct PartitionSet {
  uint64_t bitSet_[P_BITS_SIZE];

  PartitionSet() {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      bitSet_[i] = uint64_t(0);
    }
  };

  void operator=(const PartitionSet& r) {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      bitSet_[i] = r.bitSet_[i];
    }
  };

  PartitionSet(const PartitionSet& set) {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      bitSet_[i] = set.bitSet_[i];
    }
  };

  int countBit() {
    int ret = 0;
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      uint64_t set = bitSet_[i];
      ret += __builtin_popcountll(set);
    }
    return ret;
  };

  int next() {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      if (bitSet_[i] == 0) continue;
      int ret = ffsll(bitSet_[i]) - 1;
      bitSet_[i] &= ~(uint64_t(1) << ret); /* clear bit */
      return ret + i*sizeof(uint64_t)*8;
    }
    return -1;
  };

  bool has(int i) {
    return (bitSet_[i/(sizeof(uint64_t)*8)] & (uint64_t(1) << uint64_t(i%(sizeof(uint64_t)*8))));
  };

  void set(int i) {
    assert(P_BITS_SIZE > i/(sizeof(uint64_t )*8));
    bitSet_[i/(sizeof(uint64_t)*8)] |= (uint64_t(1) << uint64_t(i % (sizeof(uint64_t)*8)));
  };

  static PartitionSet OR(const PartitionSet& x, const PartitionSet& y) {
    PartitionSet ret;
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      ret.bitSet_[i] = x.bitSet_[i] | y.bitSet_[i];
    }
    return ret;
  };

  static PartitionSet XOR(const PartitionSet& x, const PartitionSet& y) {
    PartitionSet ret;
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      ret.bitSet_[i] = x.bitSet_[i] ^ y.bitSet_[i];
    }
    return ret;
  };

  static PartitionSet InterSection(const PartitionSet& x, const PartitionSet& y) {
    PartitionSet ret;
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      ret.bitSet_[i] = x.bitSet_[i] & y.bitSet_[i];
    }
    return ret;
  };

  static bool HasInterSection(const PartitionSet& x, const PartitionSet& y) {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      if (x.bitSet_[i] & y.bitSet_[i]) return true;
    }
    return false;
  };

  static int GetOneInterSection(const PartitionSet& x, const PartitionSet& y) {
    for (int i = 0; i < P_BITS_SIZE; ++i) {
      uint64_t is = x.bitSet_[i] & y.bitSet_[i];
      if (is) return (ffsll(is) - 1) + i*sizeof(uint64_t)*8;
    }
    return -1;
  }
};

#endif //DISTRIBUTEDNE_PARTITION_HPP