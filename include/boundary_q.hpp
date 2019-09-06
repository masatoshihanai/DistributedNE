/*
 * boundary_q.hpp 
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_BOUNDARY_Q_HPP
#define DISTRIBUTEDNE_BOUNDARY_Q_HPP

#include <unordered_map>
#include <map>
#include <parallel/algorithm>
#include "type.hpp"
#include "partition.hpp"

class BoundaryQueue {
  VertScores queue_;
  uint64_t idx_ = 0;

 public:
  BoundaryQueue() {};

  VertT smallestValuedID() {
    return queue_[idx_].id_;
  };

  VertT smallestValue() {
    return queue_[idx_].score_;
  };

  void popSmallest() {
    ++idx_;
  };

  VertT size() {
    return queue_.size() - idx_;
  };

  void insert(VertScores& scores) {
    auto cmp = [](const VertScore& r, const VertScore& l){ return r.score_ < l.score_;};
    std::sort(scores.begin(), scores.end(), cmp);

    VertScores newQueue_;
    newQueue_.reserve(queue_.size() + scores.size());
    uint64_t i = idx_; uint64_t j = 0;
    std::unordered_set<VertT> deplicated;
    while (i < queue_.size() || j < scores.size()) {
      if (i == queue_.size()) {
        if (deplicated.count(scores[j].id_) == 0) {
          if (scores[j].score_ != 0) {
            newQueue_.push_back(scores[j]);
          }
          deplicated.insert(scores[j].id_);
        }
        ++j;
      } else if (j == scores.size()) {
        if (deplicated.count(queue_[i].id_) == 0) {
          if (queue_[i].score_ != 0) {
            newQueue_.push_back(queue_[i]);
          }
          deplicated.insert(queue_[i].id_);
        }
        ++i;
      } else if (queue_[i].score_ > scores[j].score_) {
        if (deplicated.count(scores[j].id_) == 0) {
          if (scores[j].score_ != 0) {
            newQueue_.push_back(scores[j]);
          }
          deplicated.insert(scores[j].id_);
        }
        ++j;
      } else {
        if (deplicated.count(queue_[i].id_) == 0) {
          if (queue_[i].score_ != 0) {
            newQueue_.push_back(queue_[i]);
          }
          deplicated.insert(queue_[i].id_);
        }
        ++i;
      }
    }
    idx_ = 0;
    newQueue_.swap(queue_);
  };
};

#endif //DISTRIBUTEDNE_BOUNDARY_Q_HPP
