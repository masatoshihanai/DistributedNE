/*
 * dgraph.hpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_DGRAPH_HPP
#define DISTRIBUTEDNE_DGRAPH_HPP

#include <algorithm>
#include <atomic>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <parallel/algorithm>
#include <unordered_set>
#include <sstream>
#include <stdint.h>
#include <sys/stat.h>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "partition.hpp"
#include "com.hpp"
#include "spinlock.hpp"
#include "type.hpp"
#include "profiler.hpp"

int RSEED = 1;

class DistributedGraph {
  /* Master Graph */
  std::vector<EdgeNum> offSet_;
  std::deque<VertT> neighbor_;
  std::vector<bool> isReverseNe_;

 public:
  std::vector<std::atomic<Score>> score_;

 private:
  SpinLocks vertLocks_;

  /* Allocated Graph */
  std::deque<std::atomic<Partition>> allocPartition_; /* max means non-allocated edge */
  std::vector<PartitionSet> allocPartitionSets_;

  /* Graph info */
  VertT numGlobalVertices_;
  VertT numCompressedVertices_ = 0;

  /* one edge is represented by two directed edges */
  EdgeNum numGlobalDirectedEdges_;
  std::vector<EdgeNum> numDirectedEdgesPerRank_;
  EdgeNum countRest_;
  std::vector<EdgeNum> numRestDirectedEdgesPerRank_;

  EdgeNum numGlobalReverseEdges_ = 0;
  std::vector<EdgeNum> numReverseEdges_;

  std::vector<EdgeNum> numLocalAllocDirectedEdgesPerRank_;
  std::vector<EdgeNum> numGlobalAllocDirectedEdgesPerRank_;

  /* Edge Buffer */
  EdgeBuf edgeBuf_;

  /* MPI conf */
  int rank_;
  int rankSize_;

  /* Thread conf */
  int numThr_ = 1;

  int iter = 0;

  EdgeNum edgeNumOrg_ = 0;

  VertT randLocalID;

 public:
  std::vector<std::vector<VertT>> syncVertBuf_;
  std::vector<bool> isSync_;

 public:
  void init(char* fileName, int rank, int rankSize, int numThr);

  void getRandom(bool*& recRequest, VPBuf& expandBuffer);

  void allocEdges(VPs& rankVertPair,
                  VPBuf& dstRankPair);

  void allocInBoundary(VPs& recDstRankPair,
                       VertScoresBuf& buf);

  void updateInfo();

  EdgeBuf& edgeBuf() {
    return edgeBuf_;
  };

  VertT numVertices() {
    return numGlobalVertices_;
  };

  VertT numUniqueVertices() {
    return numCompressedVertices_  == 0 ? numGlobalVertices_ : numCompressedVertices_;
  }

  EdgeNum numEdges() {
    return numGlobalDirectedEdges_ - numGlobalReverseEdges_;
  };

  EdgeNum numOrgEdges() {
    return edgeNumOrg_;
  };

  EdgeNum numDEdges() {
    return numGlobalDirectedEdges_;
  }

  EdgeNum allocDirectedEdges(int rank) {
    return numGlobalAllocDirectedEdgesPerRank_[rank];
  };

  EdgeNum restDirectedEdges(int rank) {
    return numRestDirectedEdgesPerRank_[rank];
  };

  EdgeNum sumRestDirectedEdges() {
    EdgeNum ret = 0;
    for (int r = 0; r < rankSize_; ++r) {
      ret += restDirectedEdges(r);
    }
    return ret;
  }

  VertItr beginNeighbors(VertT id) {
    VertT localID = MetaD::GtoL(id);
    return neighbor_.begin() +
           (localID == 0 ? 0 : offSet_[localID - 1]);
  };

  VertItr endNeighbors(VertT id) {
    VertT localID = MetaD::GtoL(id);
    return neighbor_.begin() + offSet_[localID];
  };

  VertItr beginNeighbors() {
    return neighbor_.begin();
  };

  VertItr endNeighbors() {
    return neighbor_.end();
  }

  EdgeNum indexOfNeighbor(VertItr it) {
    return it - neighbor_.begin();
  };

  VertItr bSearchDst(VertT src, VertT dst) {
    VertItr ret = std::lower_bound(beginNeighbors(src),
                                   endNeighbors(src),
                                   dst);
    if (*ret == dst) return ret;
  };

  EdgeNum numReplica() {
    std::vector<EdgeNum> tmp(numThr_);

    for (EdgeNum i = 0; i < allocPartitionSets_.size(); ++i) {
      tmp[omp_get_thread_num()] += allocPartitionSets_[i].countBit();
    }

    EdgeNum ret = 0;
    for (auto x: tmp) {
      ret += x;
    }
    return ret;
  };

  DistributedGraph() {};

 private:
  void constructGraphFromFile(char* fileName);

 private:
  DistributedGraph(const DistributedGraph&);
  void operator=(const DistributedGraph&);

  /* for VERBOSE */
 public:
  void debugInfo(int iteration) {
    if (iteration == 0) {
      std::cout << " Itr|"
                << std::setw(std::max((int)std::to_string(numEdges()).length(), 11)) << " RestEdges ";
      Profiler::ProfileInit();
      std::cout << std::endl;
    }
    std::cout << std::setw(4) << iteration;

    std::cout << "|" << std::setw(std::max((int)std::to_string(numEdges()).length(), 11)) << sumRestDirectedEdges();

    Profiler::Profile();

    std::cout << std::endl;
  }
};

void DistributedGraph::getRandom(bool*& randRequest, VPBuf& expandBuffer) {
  using PartiNumPair = Edge;
  EdgeNum numRest = 0;
  std::vector<Edge> numRestEdgesSort;
  EdgeNum numTotal = 0;
  for (int r = 0; r < rankSize_; ++r) {
    numRestEdgesSort.emplace_back((VertT) numRestDirectedEdgesPerRank_[r], r);
    numRest += numRestDirectedEdgesPerRank_[r];
    numTotal += numGlobalAllocDirectedEdgesPerRank_[r];
  }
  std::sort(numRestEdgesSort.begin(), numRestEdgesSort.end(), std::greater<Edge>());
//  EdgeNum averageAllocNum = numTotal / rankSize_;

  /* Get random vertex from rank with the most # of non-alloc edges */
  for (VertT r = 0, i = 0; r < rankSize_; ++r) {
    if (randRequest[r]) {
      if (numRestEdgesSort[i].to_ == rank_) {
        /* Get from this rank */
        std::unordered_set<VertT> expanded;
        for (VertT x = 0; x < offSet_.size(); ++x) {
          randLocalID = (randLocalID + 1) % offSet_.size();
          EdgeNum edgeIndx = randLocalID == 0 ? 0 : offSet_[randLocalID - 1];
          EdgeNum ldegree = (offSet_[randLocalID] - edgeIndx);
          if (iter < 50 && ldegree > 2) continue; // skip high degree
          if (score_[randLocalID] == 0) continue; // skip 0 degree

          if (expanded.count(MetaD::LtoG(randLocalID)) == 0) { // remove duplicated
            expanded.insert(MetaD::LtoG(randLocalID));
            VertParti ir(MetaD::LtoG(randLocalID), (Partition) r);
            for (auto itr = MetaD::beginPartitions(MetaD::LtoG(randLocalID));
                 itr != MetaD::endPartitions(MetaD::LtoG(randLocalID));
                 ++itr) {
              expandBuffer[*itr].push_back(ir);
            }
            break;
          }
        }
// TODO smarter algorithm
//
//        EdgeNum ranNum = 1;
//        if (iter > 50) {
//          ranNum = averageAllocNum < allocDirectedEdges(r) ? 1 :
//              (averageAllocNum - numGlobalAllocDirectedEdgesPerRank_[r]) * 0.1;
//        }
//
//        std::unordered_set<VertT> expanded;
//        EdgeNum n = 0;
//        for (VertT x = 0; x < offSet_.size(); ++x) {
//          randLocalID = (randLocalID + 1) % offSet_.size();
//          EdgeNum edgeIndx = randLocalID == 0 ? 0 : offSet_[randLocalID - 1];
//          EdgeNum ldegree = (offSet_[randLocalID] - edgeIndx);
//          if (iter < 50 && ldegree > 2) continue; // skip high degree
//          if (score_[randLocalID] == 0) continue; // skip 0 degree
//
//          if (expanded.count(MetaD::LtoG(randLocalID)) == 0) { // remove deplicated
//            expanded.insert(MetaD::LtoG(randLocalID));
//            VertParti ir(MetaD::LtoG(randLocalID), (Partition) r);
//            for (auto itr = MetaD::beginPartitions(MetaD::LtoG(randLocalID));
//                 itr != MetaD::endPartitions(MetaD::LtoG(randLocalID));
//                 ++itr) {
//              expandBuffer[*itr].push_back(ir);
//            }
//            n += score_[randLocalID] * (MetaD::rows + MetaD::cols);
//          }
//
//          if (n > ranNum) break;
//        }
      }
      ++i;
    }
  }
};

void DistributedGraph::allocEdges(VPs& rankVertPair,
                                  VPBuf& dstRankPair) {
  std::vector<EdgeNum> countAllocEdges(numThr_, 0);
  std::vector<std::atomic<EdgeNum>> localAlloc(rankSize_);

  std::sort(rankVertPair.begin(), rankVertPair.end());

  SpinLocks locks(rankSize_);
  #pragma omp parallel for
  for (VertT i = 0; i < rankVertPair.size(); ++i) {
    VertT src = rankVertPair[i].id_;
    Partition rank = rankVertPair[i].parti_;

    if (numGlobalAllocDirectedEdgesPerRank_[rank] + localAlloc[rank] * rankSize_
        >  numGlobalDirectedEdges_ / rankSize_ * BALANCE_FACTOR) {
      continue;
    }
    if (score_[MetaD::GtoL(src)] == 0) continue;

    for (auto srcToDst = beginNeighbors(src); srcToDst != endNeighbors(src);) {
      /* allocate srcToDst edge */
      Partition max = PART_MAX;
      if (allocPartition_[indexOfNeighbor(srcToDst)].compare_exchange_weak(max, rank)) {
        /* allocate dstToSrc edge */
        VertT dst = *srcToDst;
        VertItr dstToSrc = bSearchDst(dst, src);
        if (!allocPartition_[indexOfNeighbor(dstToSrc)].compare_exchange_weak(max, rank)) {
          /* dst-to-src has already allocated */
          if (src < dst) {
            allocPartition_[indexOfNeighbor(dstToSrc)] = rank;
          } else {
            ++srcToDst;
            continue;
          }
        }

        /* buffer send edges */
        if (!isReverseNe_[indexOfNeighbor(srcToDst)]) {
          locks[rank].lock();
          edgeBuf_[rank].emplace_back(src, *srcToDst, isReverseNe_[indexOfNeighbor(srcToDst)]);
          locks[rank].unlock();
        }

        if (!isReverseNe_[indexOfNeighbor(dstToSrc)]) {
          locks[rank].lock();
          edgeBuf_[rank].emplace_back(*srcToDst, src, isReverseNe_[indexOfNeighbor(dstToSrc)]);
          locks[rank].unlock();
        }

        ++localAlloc[rank];
        ++localAlloc[rank];

        countAllocEdges[omp_get_thread_num()] += 2;
        --score_[MetaD::GtoL(dst)];
        --score_[MetaD::GtoL(src)];

        /* set partitions */
        {
          vertLocks_[MetaD::GtoL(dst)].lock();
          isSync_[MetaD::GtoL(dst)] = true;
          allocPartitionSets_[MetaD::GtoL(dst)].set(rank);
          vertLocks_[MetaD::GtoL(dst)].unlock();
        }

        VertParti pair(dst, rank);
        for (auto itr = MetaD::beginPartitions(dst);
             itr != MetaD::endPartitions(dst);
             ++itr) {
          {
            locks[*itr].lock();
            dstRankPair[*itr].push_back(pair);
            locks[*itr].unlock();
          }
        }
        ++srcToDst;
      } else {
        /* skip directed edge */
        VertT dst = *srcToDst;
        ++srcToDst;
        if (srcToDst != endNeighbors(src) && *srcToDst == dst) {
          ++srcToDst;
        }
      }
    }
  } /* end #pragma omp parallel for */

  EdgeNum num = 0;
  for (int r = 0; r < rankSize_; ++r) {
    num += dstRankPair[r].size();
  }
  VPs().swap(rankVertPair);

  for (int thr = 0; thr < numThr_; ++thr) {
    countRest_ -= countAllocEdges[thr];
  }
  for (int r = 0; r < rankSize_; ++r) {
    numLocalAllocDirectedEdgesPerRank_[r] += localAlloc[r];
  }
}; /* DistributedGraph::allocEdges() */

void DistributedGraph::allocInBoundary(VPs& recDstRankPair, VertScoresBuf& buf) {
  std::sort(recDstRankPair.begin(), recDstRankPair.end());

  /* Sync partition set */
  #pragma omp parallel for
  for (VertT i = 0; i < recDstRankPair.size(); ++i) {
    VertT dst = recDstRankPair[i].id_;
    Partition parti = recDstRankPair[i].parti_;
    {
      vertLocks_[MetaD::GtoL(dst)].lock();
      isSync_[MetaD::GtoL(dst)] = true;
      allocPartitionSets_[MetaD::GtoL(dst)].set(parti);
      vertLocks_[MetaD::GtoL(dst)].unlock();
    }
  }

  auto func = [](const VertParti& r, const VertParti& l){ return r.id_ == l.id_;};
  VPs removeDeplicated(recDstRankPair.size());
  auto itr = std::unique_copy(recDstRankPair.begin(), recDstRankPair.end(), removeDeplicated.begin(), func);
  removeDeplicated.erase(itr, removeDeplicated.end());

  /* allocate in-boundary */
  std::vector<EdgeNum> countAllocEdges(numThr_, 0);
  std::vector<std::vector<EdgeNum>> localAlloc(numThr_, std::vector<EdgeNum>(rankSize_,0));
  SpinLocks locks(rankSize_);

  #pragma omp parallel for
  for (VertT i = 0; i < removeDeplicated.size(); ++i) {
    VertT dst = removeDeplicated[i].id_;
    if (score_[MetaD::GtoL(dst)] == 0) continue;
    for (auto twohopDstItr = beginNeighbors(dst); twohopDstItr != endNeighbors(dst);) {
      if (allocPartition_[indexOfNeighbor(twohopDstItr)].load() == PART_MAX) {
        PartitionSet parti = PartitionSet::InterSection(allocPartitionSets_[MetaD::GtoL(dst)],
                                                        allocPartitionSets_[MetaD::GtoL(*twohopDstItr)]);
        /* Get min partition */
        Partition sharedRank = PART_MAX;
        EdgeNum numAlloced = numGlobalDirectedEdges_ / rankSize_ * BALANCE_FACTOR;
        int p = parti.next();
        for (; p != -1; p = parti.next()) {
          EdgeNum val = numGlobalAllocDirectedEdgesPerRank_[p]
                        + localAlloc[omp_get_thread_num()][p] * numThr_ * rankSize_;
          if (numAlloced > val) {
            sharedRank = (Partition) p;
            numAlloced = val;
          }
        }

        if (sharedRank != PART_MAX && numAlloced < numGlobalDirectedEdges_ / rankSize_ * BALANCE_FACTOR) {
          /* alloc dst to twohopDst */
          Partition max = PART_MAX;
          if (allocPartition_[indexOfNeighbor(twohopDstItr)].compare_exchange_weak(max, sharedRank)) {
            /* twohopDst to dst */
            VertT twohopDst = *twohopDstItr;
            VertItr twohopDstToDstItr = bSearchDst(twohopDst, dst);
            if (!allocPartition_[indexOfNeighbor(twohopDstToDstItr)].compare_exchange_weak(max, sharedRank)) {
              if (dst < twohopDst) {
                allocPartition_[indexOfNeighbor(twohopDstToDstItr)] = sharedRank;
              } else {
                /* skip directed edge */
                VertT twohopDst = *twohopDstItr;
                ++twohopDstItr;
                if (twohopDstItr != endNeighbors(dst) && *twohopDstItr == twohopDst) {
                  ++twohopDstItr;
                }
                continue;
              }
            }

            /* buffer send edges */
            if (!isReverseNe_[indexOfNeighbor(twohopDstItr)]) {
              locks[sharedRank].lock();
              edgeBuf_[sharedRank].emplace_back(dst, twohopDst, isReverseNe_[indexOfNeighbor(twohopDstItr)]);
              locks[sharedRank].unlock();
            }

            if (!isReverseNe_[indexOfNeighbor(twohopDstToDstItr)]) {
              locks[sharedRank].lock();
              edgeBuf_[sharedRank].emplace_back(twohopDst, dst, isReverseNe_[indexOfNeighbor(twohopDstToDstItr)]);
              locks[sharedRank].unlock();
            }

            localAlloc[omp_get_thread_num()][sharedRank] += 2;
            countAllocEdges[omp_get_thread_num()] += 2;

            if (!isSync_[MetaD::GtoL(twohopDst)]) {
              for (auto itr = MetaD::beginPartitions(twohopDst);
                   itr != MetaD::endPartitions(twohopDst);
                   ++itr) {
                {
                  locks[*itr].lock();
                  syncVertBuf_[*itr].push_back(twohopDst);
                  locks[*itr].unlock();
                }
              }
            }

            --score_[MetaD::GtoL(dst)];
            --score_[MetaD::GtoL(twohopDst)];

            ++twohopDstItr;
            /* Directed edge */
            if (twohopDstItr != endNeighbors(dst) && *twohopDstItr == twohopDst) {
              allocPartition_[indexOfNeighbor(twohopDstItr)] = sharedRank;
              if (!isReverseNe_[indexOfNeighbor(twohopDstItr)]) {
                locks[sharedRank].lock();
                edgeBuf_[sharedRank].emplace_back(dst, twohopDst, isReverseNe_[indexOfNeighbor(twohopDstItr)]);
                locks[sharedRank].unlock();
              }

              ++twohopDstToDstItr;
              allocPartition_[indexOfNeighbor(twohopDstToDstItr)] = sharedRank;
              if (!isReverseNe_[indexOfNeighbor(twohopDstToDstItr)]) {
                locks[sharedRank].lock();
                edgeBuf_[sharedRank].emplace_back(twohopDst, dst, isReverseNe_[indexOfNeighbor(twohopDstToDstItr)]);
                locks[sharedRank].unlock();
              }
              ++twohopDstItr;
              localAlloc[omp_get_thread_num()][sharedRank] += 2;
              countAllocEdges[omp_get_thread_num()] += 2;
              --score_[MetaD::GtoL(dst)];
              --score_[MetaD::GtoL(twohopDst)];
            }
          } else {
            /* skip directed edge */
            VertT twohopDst = *twohopDstItr;
            ++twohopDstItr;
            if (twohopDstItr != endNeighbors(dst) && *twohopDstItr == twohopDst) {
              ++twohopDstItr;
            }
          }
        } else {
          /* skip directed edge */
          VertT twohopDst = *twohopDstItr;
          ++twohopDstItr;
          if (twohopDstItr != endNeighbors(dst) && *twohopDstItr == twohopDst) {
            ++twohopDstItr;
          }
        }
      } else {
        /* skip directed edge */
        VertT twohopDst = *twohopDstItr;
        ++twohopDstItr;
        if (twohopDstItr != endNeighbors(dst) && *twohopDstItr == twohopDst) {
          ++twohopDstItr;
        }
      }
    }
  } /* end #pragma omp parallel for */

  for (int thr = 0; thr < numThr_; ++thr) {
    countRest_ -= countAllocEdges[thr];
    for (int r = 0; r < rankSize_; ++r) {
      numLocalAllocDirectedEdgesPerRank_[r] += localAlloc[thr][r];
    }
  }

  /* Sync Score */
  auto f = [](const VertParti& r, const VertParti& l){ return r.id_ == l.id_ && r.parti_ == l.parti_;};
  VPs syncVert(recDstRankPair.size());
  auto i = std::unique_copy(recDstRankPair.begin(), recDstRankPair.end(), syncVert.begin(), f);
  syncVert.erase(i, syncVert.end());

  #pragma omp parallel for
  for (VertT i = 0; i < syncVert.size(); ++i) {
    VertT dst = syncVert[i].id_;
    Partition parti = syncVert[i].parti_;
    if (score_[MetaD::GtoL(dst)] == 0) continue;
    locks[parti].lock();
    VertScore vs(dst, score_[MetaD::GtoL(dst)]);
    buf[parti].push_back(vs);
    locks[parti].unlock();
  }

  VPs().swap(recDstRankPair);
}; /* DistributedGraph::allocInBoundary() */

void DistributedGraph::updateInfo() {
  for (auto& buf: edgeBuf_) {
    buf.clear();
  }

  MPI_Allreduce(&numLocalAllocDirectedEdgesPerRank_[0],
                &numGlobalAllocDirectedEdgesPerRank_[0],
                rankSize_,
                MPI_EDGE_NUM,
                MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allgather(&countRest_,
                1,
                MPI_EDGE_NUM,
                &numRestDirectedEdgesPerRank_[0],
                1,
                MPI_EDGE_NUM,
                MPI_COMM_WORLD);
  ++iter;

  /* compute rest */
  auto min = std::min_element(numGlobalAllocDirectedEdgesPerRank_.begin(),
                              numGlobalAllocDirectedEdgesPerRank_.end());
  EdgeNum sum = 0;
  for (auto& x: numRestDirectedEdgesPerRank_) { sum += x; }

  bool allocRest = true;
  for (auto it = numGlobalAllocDirectedEdgesPerRank_.begin();
       it != numGlobalAllocDirectedEdgesPerRank_.end(); ++it) {
    if (it != min && *it < (numGlobalDirectedEdges_ / rankSize_)) {
      allocRest = false;
    }
  }

  if (allocRest || sum + *min < (numGlobalDirectedEdges_ / rankSize_)) {
    int targetRank = min - numGlobalAllocDirectedEdgesPerRank_.begin();
    for (VertT localID = 0; localID < offSet_.size(); ++localID) {
      EdgeNum edgeItr = localID == 0 ? 0 : offSet_[localID - 1];
      while (edgeItr < offSet_[localID]) {
        if (allocPartition_[edgeItr].load() == PART_MAX) {
          allocPartition_[edgeItr] = targetRank;
          ++numLocalAllocDirectedEdgesPerRank_[targetRank];
          --countRest_;
          if (!isReverseNe_[edgeItr]) {
            edgeBuf_[targetRank].emplace_back(MetaD::LtoG(localID), neighbor_[edgeItr], isReverseNe_[edgeItr]);
          }
        }
        ++edgeItr;
      }
    }
  }
}; /* DistributedGraph::updateInfo() */

void DistributedGraph::init(char* fileName, int rank, int rankSize, int numThr) {
  rank_ = rank;
  rankSize_ = rankSize;
  numThr_ = numThr;

  if (!VERBOSE) std::ios::sync_with_stdio(false);
  constructGraphFromFile(fileName);

  /* Init allocated Partition */
  allocPartition_ = std::deque<std::atomic<Partition>>(neighbor_.size());
  std::fill(allocPartition_.begin(), allocPartition_.end(), PART_MAX);
  allocPartitionSets_.resize(offSet_.size());
  /* init score */
  score_ = std::vector<std::atomic<Score>>(offSet_.size());
  score_[0] = offSet_[0];
  for (VertT i = 1; i < offSet_.size(); ++i) {
    score_[i] = offSet_[i] - offSet_[i - 1];
  }

  /* init lock */
  vertLocks_ = SpinLocks(offSet_.size());
  isSync_ = std::vector<bool>(offSet_.size(), false);

  /* Sync num reverse edges */
  EdgeNum num = numReverseEdges_[rank_];
  MPI_Allgather(&num,
                1,
                MPI_EDGE_NUM,
                &numReverseEdges_[0],
                1,
                MPI_EDGE_NUM,
                MPI_COMM_WORLD);
  numGlobalReverseEdges_ = 0;
  for (auto x: numReverseEdges_) {
    numGlobalReverseEdges_ += x;
  }

  /* init global info */
  numDirectedEdgesPerRank_.resize(rankSize);
  EdgeNum localNumMasterEdges = neighbor_.size();
  MPI_Allgather(&localNumMasterEdges,
                1,
                MPI_EDGE_NUM,
                &numDirectedEdgesPerRank_[0],
                1,
                MPI_EDGE_NUM,
                MPI_COMM_WORLD);

  numGlobalDirectedEdges_ = 0;
  for (auto& x: numDirectedEdgesPerRank_) { numGlobalDirectedEdges_ += x; };

  countRest_ = neighbor_.size();

  numRestDirectedEdgesPerRank_.insert(numRestDirectedEdgesPerRank_.end(),
                                      numDirectedEdgesPerRank_.begin(),
                                      numDirectedEdgesPerRank_.end());

  numLocalAllocDirectedEdgesPerRank_.resize(rankSize);
  numGlobalAllocDirectedEdgesPerRank_.resize(rankSize);

  /* init buffers */
  edgeBuf_ = EdgeBuf(rankSize);

  syncVertBuf_ = std::vector<std::vector<VertT>>(rankSize);

  /* Random Seed */
  srand(RSEED);
  randLocalID = rank_ == numGlobalVertices_ ? rank_ : (rand() * rank_) % offSet_.size();

  /* Profile Mem Usage */
  if (VERBOSE && rank_ == 0) {
    std::cout << "---Mem Usage for Graph (BYTE) ---" << std::endl;
    uint64_t vMem = offSet_.size() * sizeof(VertT);
    uint64_t eMem = neighbor_.size() * sizeof(VertT);
    uint64_t sMem = score_.size() * sizeof(Score);
    uint64_t lMem = vertLocks_.size() * sizeof(SLock);
    uint64_t aMem = allocPartition_.size() * sizeof(Partition);
    uint64_t pMem = allocPartitionSets_.size() * sizeof(PartitionSet);
    int size = std::to_string(vMem + eMem + sMem + lMem + aMem + pMem).length();
    std::cout << "  Total       : " << std::setw(size) << (vMem + eMem + sMem + lMem + aMem + pMem) << std::endl;
    std::cout << "    OffSet    : " << std::setw(size) << vMem << std::endl;
    std::cout << "    Neighbor  : " << std::setw(size) << eMem << std::endl;
    std::cout << "    Score     : " << std::setw(size) << sMem << std::endl;
    std::cout << "    Locks     : " << std::setw(size) << lMem << std::endl;
    std::cout << "    Alloc edge: " << std::setw(size) << aMem << std::endl;
    std::cout << "    Alloc sets: " << std::setw(size) << pMem << std::endl;
    std::cout << "---------------------------------" << std::endl;
  }
  return;
}; /* DistributedGraph::init() */

void DistributedGraph::constructGraphFromFile(char* fileName) {
  int rank = rank_;
  int rankSize = rankSize_;
  char* fileData;
  uint64_t startOffset;
  if (rank_ == 0) std::cout << "  Open file: " << fileName << std::endl;
  /* Read File */
  {
    MPI_File file;
    int error = MPI_File_open(MPI_COMM_WORLD,
                              fileName,
                              MPI_MODE_RDONLY,
                              MPI_INFO_NULL,
                              &file);
    if (error) {
      std::cerr << "fail to open file: " << fileName << std::endl;
      exit(1);
    }
    MPI_Offset totalFileSize, localScanSize;
    MPI_File_get_size(file, &totalFileSize);

    /* Set start point of each rank */
    localScanSize = totalFileSize / rankSize_;
    MPI_Offset localStartOffset = rank_ * localScanSize;

    /* Set tail of each rank */
    const int MAX_LINE_NUM = 30; // todo tuning for large scale
    MPI_Offset localEndOffset;
    if (rank_ == rankSize_ - 1) {
      /* final rank will scan to the end */
      localEndOffset = totalFileSize;
    } else {
      /* other rank will scan to the end of line */
      /* the end of line ('\n') must be included in MAX_LINE_NUM */
      localEndOffset = localStartOffset + localScanSize - 1 + MAX_LINE_NUM;
    }

    MPI_Offset localSize = localEndOffset - localStartOffset + 1;
    fileData = new char[localSize + 1];

    /* Read file in parallel */
    MPI_File_close(&file);
    std::fstream filess;
    filess.open(fileName, std::fstream::in);
    filess.seekg(localStartOffset, filess.beg);
    filess.read(fileData, localEndOffset - localStartOffset);

    /* init start point */
    size_t startLines = 0;
    if (rank_ != 0) {
      while (fileData[startLines] != '\n') {
        ++startLines;
      }
      ++startLines;
    }

    /* init end point */
    size_t endLines = localSize;
    if (rank_ != rankSize_ - 1) {
      endLines -= MAX_LINE_NUM;
      while (fileData[endLines] != '\n') {
        ++endLines;
      }
    }
    fileData[endLines + 1] = '\0';

    startOffset = startLines;

    filess.close();
  } /* Read File */

  /* Scan Edges */
  EdgeBuf edgeSendBuf(rankSize);
  VertBuf vertexBuffer(rankSize);
  size_t fileLen = strlen(fileData);
  if (VERBOSE && rank_ == 0) {
    std::cout << "  Finish reading file: size (Byte)" << fileLen << std::endl;
  }

  for (auto& sndBuf: edgeSendBuf) {
    sndBuf.reserve(fileLen / 15 / rankSize);
  }

  VertT localNumEdges = 0;
  VertT scale = 14;
  VertT localMaxVertID = 0;

  EdgeNum lEdgeNumOrg = 0;
  char newline[] = "\n";
  for (char* p = strtok(&fileData[startOffset], newline);
       p != NULL;
       p = strtok(NULL, newline)) {
    if (p[0] == '#' || p[0] == '%') {
      continue;
    }
    ++lEdgeNumOrg;
    std::istringstream iss(p);
    VertT from, to;
    if (!(iss >> from >> to)) continue;
    if (from == to) continue;
    bool reverse = false;
    if (from > to) {
      VertT tmp = from;
      from = to;
      to = tmp;
      reverse = true;
    }

    if (localMaxVertID < to) {
      localMaxVertID = to;
    }

    Edge edgeline;
    edgeline.from_ = from;
    edgeline.to_ = to;
    edgeline.reverse_ = reverse;
    Partition partition = MetaD::getMaster(edgeline);
    edgeSendBuf[partition].push_back(edgeline);
    ++localNumEdges;

    if (VERBOSE) {
      if (rank_ == 0 && localNumEdges % (VertT) pow(2, scale) == 0) {
        std::cout << "\r  Read edges counts: " << localNumEdges << " / " << fileLen / 15 << std::flush;
        ++scale;
      }
    }
  }
  if (VERBOSE && rank == 0) {
    std::cout << std::endl;
    std::cout << "  Finish reading edges: " << localNumEdges << " at rank 0" << std::endl;
  }
  delete[] fileData;
  MPI_Allreduce(&lEdgeNumOrg, &edgeNumOrg_, 1, MPI_EDGE_NUM, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localMaxVertID, &numGlobalVertices_, 1, MPI_VERT_T, MPI_MAX, MPI_COMM_WORLD);

  /* shuffle edges */
  std::vector<Edge> recEdgeBuf;
  {
    Communicator com;
    com.init(rank_, rankSize_);
    com.AlltoAll(edgeSendBuf, recEdgeBuf);
  }
  for (int r = 0; r < rankSize_; ++r) { edgeSendBuf.at(r).clear(); }
  edgeSendBuf.clear();
  MPI_Barrier(MPI_COMM_WORLD);

  /* construct graph */
  std::sort(recEdgeBuf.begin(), recEdgeBuf.end());
  std::vector<VertT> lReverseCount(numThr_);

  #pragma omp parallel for
  for (uint64_t i = 0; i < recEdgeBuf.size(); ++i) {
    if (i != 0 && recEdgeBuf[i] == recEdgeBuf[i-1]) {
      ++lReverseCount[omp_get_thread_num()];
    }
  }

  VertT reverseCount = 0;
  for (auto x: lReverseCount) {
    reverseCount += x;
  }

  std::deque<Edge> edgeWithReverse;
  //edgeWithReverse.reserve((recEdgeBuf.size() - reverseCount)*2);
  for (uint64_t i = 0; i < recEdgeBuf.size(); ++i) {
    if (i != 0 && recEdgeBuf[i] == recEdgeBuf[i-1]) {
      /* update by original edge */
      edgeWithReverse.back().reverse_ = false;
    } else {
      if (recEdgeBuf[i].reverse_) {
        /* original edge */
        Edge e(recEdgeBuf[i].to_, recEdgeBuf[i].from_, false);
        edgeWithReverse.push_back(e);
        /* reverse edge */
        Edge re(recEdgeBuf[i].from_, recEdgeBuf[i].to_, true);
        edgeWithReverse.push_back(re);
      } else {
        /* original edge */
        Edge e(recEdgeBuf[i].from_, recEdgeBuf[i].to_, false);
        edgeWithReverse.push_back(e);
        /* reverse edge */
        Edge re(recEdgeBuf[i].to_, recEdgeBuf[i].from_, true);
        edgeWithReverse.push_back(re);
      }
    }
  }
  recEdgeBuf.clear();

  std::sort(edgeWithReverse.begin(), edgeWithReverse.end());

  neighbor_.resize(edgeWithReverse.size());
  isReverseNe_.resize(edgeWithReverse.size());
  std::vector<EdgeNum> localNum(numThr_);

  #pragma omp parallel for
  for (uint64_t i = 0; i < edgeWithReverse.size(); ++i) {
    neighbor_[i] = edgeWithReverse[i].to_;
    isReverseNe_[i] = edgeWithReverse[i].reverse_;
    if (edgeWithReverse[i].reverse_) ++localNum[omp_get_thread_num()];
  }
  numReverseEdges_.resize(rankSize);
  for (auto x: localNum) {
    numReverseEdges_[rank_] += x;
  }

  EdgeNum num = numReverseEdges_[rank_];
  MPI_Allgather(&num,
                1,
                MPI_EDGE_NUM,
                &numReverseEdges_[0],
                1,
                MPI_EDGE_NUM,
                MPI_COMM_WORLD);
  numGlobalReverseEdges_ = 0;
  for (auto x: numReverseEdges_) {
    numGlobalReverseEdges_ += x;
  }

  VertT prev = VERT_MAX;
  offSet_.resize(MetaD::GtoL(numGlobalVertices_) + 1);
  uint64_t lid = 0; /* for offset */
  for (Edge& edge: edgeWithReverse) {
    if (edge.from_ != prev) {
      /* move to the vertex */
      while(lid != MetaD::GtoL(edge.from_)) {
        offSet_[lid+1] = offSet_[lid];
        ++lid;
      }
      ++offSet_[MetaD::GtoL(edge.from_)];
    } else {
      ++offSet_[MetaD::GtoL(edge.from_)];
    }
    prev = edge.from_;
  }
  /* fill until end */
  while (lid != MetaD::GtoL(numGlobalVertices_)) {
    offSet_[lid+1] = offSet_[lid];
    ++lid;
  }

  edgeWithReverse.clear();
  edgeWithReverse.shrink_to_fit();

  if (VERBOSE && rank_ == 0) {
    std::cout << "  Finish construct graph" << std::endl;
  }
}/* DistributedGraph::constructGraphFromFile() */

#endif //DISTRIBUTEDNE_DGRAPH_HPP
