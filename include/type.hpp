/*
 * type.hpp 
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_TYPE_HPP
#define DISTRIBUTEDNE_TYPE_HPP

#include <iostream>
#include <cmath>
#include <deque>
#include <math.h>
#include <limits>
#include <vector>

bool VERBOSE = false;
bool VERTEX_BALANCE = false;
double BALANCE_FACTOR = 1.01;

#ifdef UINT32_VERTT
using VertT = uint32_t;
#define MPI_VERT_T MPI_UINT32_T
#else
using VertT = uint64_t;
#define MPI_VERT_T MPI_UINT64_T
#endif

const static VertT VERT_MAX = std::numeric_limits<VertT>::max();
using Score = VertT;
using Partition = uint16_t;
const static Partition PART_MAX = std::numeric_limits<Partition>::max();


using VertVec = std::vector<VertT>;
using VertItr = typename std::deque<VertT>::iterator;
using VertBuf = std::vector<std::vector<VertT>>;

struct VPair {
  VertT first_;
  VertT second_;

  VPair() {};
  VPair(VertT first, VertT second): first_(first), second_(second) {};
  VPair(const VPair& pair) {
    first_ = pair.first_;
    second_ = pair.second_;
  }

  bool operator<(const VPair& r) const {
    if (first_ == r.first_) {
      return second_ < r.second_;
    }
    return first_ < r.first_;
  }

  bool operator>(const VPair& r) const {
    if (first_ == r.first_) {
      return second_ > r.second_;
    }
    return first_ > r.first_;
  }

  void operator=(const VPair& r) {
    first_ = r.first_;
    second_ = r.second_;
  }
};
using VPairs = std::vector<VPair>;
using VPairsBuf = std::vector<VPairs>;

struct VertParti {
  VertT id_;
  Partition parti_;

  VertParti() {};
  VertParti(VertT first, Partition second): id_(first), parti_(second) {};
  VertParti(const VertParti& vp) {
    id_ = vp.id_;
    parti_ = vp.parti_;
  }

  bool operator<(const VertParti& r) const {
    return id_ < r.id_;
  }

  bool operator>(const VertParti& r) const {
    return id_ > r.id_;
  }

  bool operator==(const VertParti& r) {
    return id_ == r.id_;
  }

  void operator=(const VertParti& r) {
    id_ = r.id_;
    parti_ = r.parti_;
  }
};
using VPs = std::vector<VertParti>;
using VPBuf = std::vector<VPs>;

struct VertScore {
  VertT id_;
  Score score_;

  VertScore() {};
  VertScore(VertT first, Score second)
      : id_(first), score_(second) {};
  VertScore(const VertScore& r)
      : id_(r.id_), score_(r.score_){};

  bool operator<(const VertScore& r) const {
    return id_ < r.id_;
  }

  bool operator==(const VertScore& r) const {
    return id_ == r.id_;
  }

  void operator=(const VertScore& r) {
    id_ = r.id_;
    score_ = r.score_;
  }
};
using VertScores = std::vector<VertScore>;
using VertScoresBuf = std::vector<std::vector<VertScore>>;

struct Edge {
  VertT from_;
  VertT to_;
  bool reverse_ = false;

  Edge() {};
  Edge(VertT from, VertT to): from_(from), to_(to), reverse_(false) {};
  Edge(VertT from, VertT to, bool isReverse): from_(from), to_(to), reverse_(isReverse) {};

  bool operator<(const Edge& r) const {
    if (from_ == r.from_ && to_ == r.to_) {
      return reverse_;
    } else if (from_ == r.from_) {
      return to_ < r.to_;
    } else {
      return from_ < r.from_;
    }
  }

  bool operator>(const Edge& r) const {
    if (from_ == r.from_ && to_ == r.to_) {
      return r.reverse_;
    } else if (from_ == r.from_) {
      return to_ > r.to_;
    } else {
      return from_ > r.from_;
    }
  }

  bool operator==(const Edge& r) const {
    return (from_ == r.from_) && (to_ == r.to_);
  }

  bool operator!=(const Edge& r) const {
    return (from_ != r.from_) || (to_ != r.to_);
  }

  void operator=(const Edge& r) {
    to_ = r.to_;
    from_ = r.from_;
    reverse_ = r.reverse_;
    return;
  }
} __attribute__((packed));
using EdgeVec = typename std::vector<Edge>;
using EdgeItr = typename std::vector<Edge>::iterator;
using EdgeBuf = typename std::vector<std::vector<Edge>>;
using EdgeNum = uint64_t;
const static EdgeNum ENUM_MAX = std::numeric_limits<EdgeNum>::max();
#define MPI_EDGE_NUM MPI_UINT64_T

struct MetaD {
  static int rankSize_;
  static int thisRank_;
  static VertT smallestCommon_;
  static std::vector<VertT> pGtoL_;
  static Partition cols;
  static Partition rows;
  static Partition lcm;

  static void init(int rankSize, int rank) {
    rankSize_ = rankSize;
    thisRank_ = rank;
    /* rows is always bigger than cols */
    rows = std::ceil(std::sqrt(rankSize_));
    while (rankSize_ % rows != 0) {
      ++rows;
    }
    cols = rankSize_ / rows;

    /* get lowest common multiple */
    lcm = cols;
    for (; lcm < rankSize_; ++lcm) {
      if (lcm % cols == 0 && lcm % rows == 0) {
        break;
      }
    }

    /* get smallest common global id between row's and col's direction */
    for (smallestCommon_ = 0; smallestCommon_ < rankSize_; ++smallestCommon_) {
      if (smallestCommon_ % cols == rank / rows) {
        if (smallestCommon_ % rows == rank % rows) {
          break;
        }
      }
      if (smallestCommon_ % rows == rank % rows) {
        if (smallestCommon_ % cols == rank / rows) {
          break;
        }
      }
    }

    /* get particular solution of GtoL */
    VertT lid = VERT_MAX;
    VertT prev = VERT_MAX;
    pGtoL_.reserve(rankSize_);
    for (int vid = 0; vid < rankSize_; ++vid) {
      lid = GtoL(vid);
      if (prev != lid) {
        pGtoL_.push_back(vid);
      }
      prev = lid;
    }
  };

  static inline Partition getMaster(const Edge& edge) {
    /* 2D-Allocation */
    Partition col = edge.from_ % cols;
    Partition row = edge.to_ % rows;
    return col * rows + row;
  };

  static inline Partition getMaster(const VertT& from, const VertT& to) {
    /* 2D-Allocation */
    Partition col = from % cols;
    Partition row = to % rows;
    return col * rows + row;
  };

  static inline VertT GtoL(VertT gid) {
    return (gid - thisRank_ % rows + rows) / rows
           + (gid - thisRank_ / rows + cols) / cols
           - 1
           - (smallestCommon_ == rankSize_ ? 0 : (gid - smallestCommon_ + lcm) / lcm);
  };

  static inline VertT LtoG(VertT lid) {
    return (lid / pGtoL_.size()) * rankSize_
           + pGtoL_[lid % pGtoL_.size()];
  };

  class PItr {
    bool isColDirection;
    Partition counter;
    VertT vid;

   public:
    Partition operator*() {
      if (isColDirection) {
        return (vid % rows) + counter*rows;
      } else {
        return (vid % cols)*rows + counter;
      }
    };

    void operator++() {
      ++counter;
      if (isColDirection && counter == cols) {
        isColDirection = false;
        counter = 0;
        return;
      }
    };

    bool operator!=(const PItr& l) {
      return isColDirection != l.isColDirection || counter != l.counter;
    };

    PItr() {
      vid = 0;
      counter = 0;
      isColDirection = true;
    };

    PItr(VertT id, Partition c, bool b)
        : vid(id), counter(c), isColDirection(b) {};
  };

  static PItr beginPartitions(VertT vid) {
    return PItr(vid, 0, true);
  };

  static PItr endPartitions(VertT vid) {
    return PItr(vid, rows, false);
  };
};

int MetaD::rankSize_;
int MetaD::thisRank_;
VertT MetaD::smallestCommon_;
std::vector<VertT> MetaD::pGtoL_;
Partition MetaD::cols;
Partition MetaD::rows;
Partition MetaD::lcm;

#endif //DISTRIBUTEDNE_TYPE_HPP
