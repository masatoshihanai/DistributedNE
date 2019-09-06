/*
 * allocator.hpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_ALLOCATOR_HPP
#define DISTRIBUTEDNE_ALLOCATOR_HPP

#include <fstream>

#include "type.hpp"
#include "dgraph.hpp"
#include "com.hpp"
#include "profiler.hpp"

class DistributedAllocator {
  Communicator com_;
  DistributedGraph* graph_;

  int rank_;
  int rankSize_;

  VPBuf exbuf_;
  bool expandRandom_ = false;
  bool* recRandomRequest_;

  VPs recVertParti;

  VPBuf dstRankPair_;
  VPs recDstRankPair_;

  VertScoresBuf sndScores_;
  VertScores recScores_;

 public:
  std::deque<Edge> edges_;

 public:
  void init(int rank, int rankSize, DistributedGraph* graph);
  void expandRandom();
  void bufferExpandVertex(VertT id);
  void expandBoundary();
  void updateBoundary(BoundaryQueue& bounday,
                      std::unordered_set<VertT>& replications);
  void clearBuffer();

  bool checkResult();
  void outputResult(const std::string& fileName);
  EdgeNum numReplica();
  EdgeNum numLocalReplica();
  EdgeNum numInitReplica();
  double balanceQuality();

  DistributedAllocator(){};
  ~DistributedAllocator() { delete[] recRandomRequest_; };
 private:
  DistributedAllocator(const DistributedAllocator&);
  void operator=(const DistributedAllocator&);
};

void DistributedAllocator::init(int rank, int rankSize, DistributedGraph* graph) {
  rank_ = rank;
  rankSize_ = rankSize;
  graph_ = graph;
  com_.init(rank, rankSize);
  recRandomRequest_ = new bool[rankSize];

  /* Reserve Buffer */
  exbuf_.resize(rankSize);
  dstRankPair_.resize(rankSize);
  sndScores_.resize(rankSize);
};

void DistributedAllocator::expandRandom() {
  expandRandom_ = true;
};

void DistributedAllocator::bufferExpandVertex(VertT id) {
  for (auto itr = MetaD::beginPartitions(id);
       itr != MetaD::endPartitions(id);
       ++itr) {
    exbuf_[*itr].emplace_back(id, rank_);
  }
};

void DistributedAllocator::expandBoundary() {
  /* Send expand vertices */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  com_.allGatherBool(expandRandom_, recRandomRequest_);
  graph_->getRandom(recRandomRequest_, exbuf_);
  com_.AlltoAll<VertParti>(exbuf_, recVertParti);
  if (VERBOSE) {
    Profiler::Expand += (MPI_Wtime() - Profiler::Start);
  }

  /* Alloc edges  */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  graph_->allocEdges(recVertParti, dstRankPair_);
  if (VERBOSE) {
    MPI_Barrier(MPI_COMM_WORLD);
    Profiler::Alloc += (MPI_Wtime() - Profiler::Start);
  }

  /* Rank sets synchronization */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  com_.AlltoAll<VertParti>(dstRankPair_, recDstRankPair_);
  if (VERBOSE) Profiler::ComRanks += (MPI_Wtime() - Profiler::Start);

  /* Alloc in-boundary edges */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  VertScoresBuf syncVertScore(rankSize_);
  graph_->allocInBoundary(recDstRankPair_, syncVertScore);
  if (VERBOSE) {
    MPI_Barrier(MPI_COMM_WORLD);
    Profiler::AllocInBoundary += (MPI_Wtime() - Profiler::Start);
  }

  /* Compute Score */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  com_.AlltoAll<VertScore>(syncVertScore, recScores_);
  if (VERBOSE) {
    Profiler::ComScore += (MPI_Wtime() - Profiler::Start);
  }

  /* Communicate Edges */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  std::vector<Edge> rec;
  com_.AlltoAll<Edge>(graph_->edgeBuf(), rec);
  edges_.insert(edges_.end(), rec.begin(), rec.end());
  if (VERBOSE) {
    MPI_Barrier(MPI_COMM_WORLD);
    Profiler::ComEdges += (MPI_Wtime() - Profiler::Start);
  }


  /* Synchronize Global Info */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  graph_->updateInfo();
  if (VERBOSE) Profiler::UpdateInf += (MPI_Wtime() - Profiler::Start);

  /* clear buffer */
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  clearBuffer();
  if (VERBOSE) Profiler::ClearBuf += (MPI_Wtime() - Profiler::Start);
}; /* DistributedAllocator::expandBoundary() */

void DistributedAllocator::updateBoundary(BoundaryQueue& boundary,
                                          std::unordered_set<VertT>& replications) {
  if (VERBOSE) {
    MPI_Barrier(MPI_COMM_WORLD);
    Profiler::Start = MPI_Wtime();
  }
  if (recScores_.size() == 0) {
    if (VERBOSE) {
      MPI_Barrier(MPI_COMM_WORLD);
      Profiler::UpdateBoundary += (MPI_Wtime() - Profiler::Start);
    }
    return;
  }
  std::sort(recScores_.begin(), recScores_.end());

  VertT prevId = recScores_[0].id_;
  Score scr = recScores_[0].score_;
  if (VERBOSE || VERTEX_BALANCE) replications.insert(prevId);
  std::vector<VertScore> vec;
  vec.reserve(recScores_.size());
  for (VertT i = 1; i < recScores_.size(); ++i) {
    if (recScores_[i].id_ == prevId) {
      scr += recScores_[i].score_;
    } else {
      vec.push_back(VertScore(prevId, scr));
      prevId = recScores_[i].id_;
      scr = recScores_[i].score_;
      if (VERBOSE || VERTEX_BALANCE) replications.insert(prevId);
    }
  }
  vec.push_back(VertScore(prevId, scr));
  if (VERBOSE || VERTEX_BALANCE) replications.insert(prevId);

  boundary.insert(vec);
  VertScores().swap(recScores_);
  if (VERBOSE) {
    MPI_Barrier(MPI_COMM_WORLD);
    Profiler::UpdateBoundary += (MPI_Wtime() - Profiler::Start);
  }
};

void DistributedAllocator::clearBuffer() {
  expandRandom_ = false;
  for (int r = 0; r < rankSize_; ++r) {
    exbuf_[r].clear();
    recRandomRequest_[r] = false;
  }

  for (int r = 0; r < rankSize_; ++r) {
    dstRankPair_[r].clear();
    sndScores_[r].clear();
  }

  recDstRankPair_.clear();
  recVertParti.clear();
};

bool DistributedAllocator::checkResult() {
  EdgeNum total;
  EdgeNum edgeCount = 0;
  edgeCount = edges_.size();

  MPI_Allreduce(&edgeCount, &total, 1, MPI_EDGE_NUM, MPI_SUM, MPI_COMM_WORLD);

  if (total != graph_->numEdges()) {
    if (rank_ == 0) {
      std::cerr << "Error total edges: " << total << " actual: " << graph_->numEdges() << std::endl;
    }
    return false;
  }
  // TODO more checks
  return true;
}

void DistributedAllocator::outputResult(const std::string& fName) {
  if (VERBOSE) Profiler::Start = MPI_Wtime();
  /* use std::io */
  std::string outfile_name = fName + "." + std::to_string(rankSize_) + ".pedges";
  if (rank_ == 0) {
    std::cout << "<<<< Start to Output Result to  " << outfile_name << " >>>>"<< std::endl;
    std::cout << " Rank ... ";
  }
  for (int r = 0; r < rankSize_; ++r) {
    if (r == rank_) {
      std::cout << " " << r << std::flush;
      std::fstream oFile;
      if (r == 0) {
        oFile.open(outfile_name, std::ios::out);
        oFile << "# SrcID  DstID  PartitionID \n";
      } else {
        oFile.open(outfile_name, std::ios::app);
      }
      for (auto it = edges_.begin(); it != edges_.end(); ++it) {
        oFile << it->from_ << " " << it->to_ << " " << rank_ << "\n";
      }
      oFile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (rank_ == 0) {
    std::cout << " ... Finish. " << std::endl;
    std::cout << "<<<< The Partitioned Graph is in " << outfile_name << " >>>>"<< std::endl;
  }
  if (VERBOSE) Profiler::Output = (MPI_Wtime() - Profiler::Start);
};

EdgeNum DistributedAllocator::numReplica() {
  std::unordered_set<VertT> verts;
  verts.reserve(graph_->numVertices() / rankSize_);

  for (size_t i = 0; i < edges_.size(); ++i) {
    verts.insert(edges_[i].from_);
    verts.insert(edges_[i].to_);
  }

  EdgeNum localReplica = verts.size();
  EdgeNum numReplica;
  MPI_Allreduce(&localReplica, &numReplica, 1, MPI_EDGE_NUM, MPI_SUM, MPI_COMM_WORLD);

  return numReplica;
};

EdgeNum DistributedAllocator::numLocalReplica() {
  std::unordered_set<VertT> verts;
  verts.reserve(graph_->numVertices() / rankSize_);
  for (size_t i = 0; i < edges_.size(); ++i) {
    verts.insert(edges_[i].from_);
    verts.insert(edges_[i].to_);
  }
  return verts.size();
};

EdgeNum DistributedAllocator::numInitReplica() {
  std::unordered_set<VertT> verts;
  verts.reserve(graph_->numVertices() / rankSize_);
  for (auto it = graph_->beginNeighbors(); it != graph_->endNeighbors(); ++it) {
    verts.insert(*it);
  }
  EdgeNum localReplica = verts.size();
  EdgeNum numReplica;
  MPI_Allreduce(&localReplica, &numReplica, 1, MPI_EDGE_NUM, MPI_SUM, MPI_COMM_WORLD);

  return numReplica;
};

double DistributedAllocator::balanceQuality() {
  EdgeNum maxEdge = 0;
  EdgeNum local =  graph_->allocDirectedEdges(rank_);
  MPI_Allreduce(&local, &maxEdge, 1, MPI_EDGE_NUM, MPI_MAX, MPI_COMM_WORLD);

  EdgeNum average = 0;
  for (int r = 0; r < rankSize_; ++r) {
    average += graph_->allocDirectedEdges(r);
  }
  average = average / rankSize_;

  return (double) maxEdge / average;
};

#endif //DISTRIBUTEDNE_ALLOCATOR_HPP
