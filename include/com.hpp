/*
 * com.hpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_COM_HPP
#define DISTRIBUTEDNE_COM_HPP

#include "type.hpp"

class Communicator {
  int rank_;
  int rankSize_;

 public:
  template<class T>
  void AlltoAll(std::vector<std::vector<T>>& snd,
                std::vector<T>& rec);


  void allGatherBool(bool requestRandom,
                     bool*& recRequest);

  void init (int rank, int rankSize) {
    rank_ = rank;
    rankSize_ = rankSize;
  }

  Communicator(){};
 private:
  Communicator(const Communicator&);
  void operator=(const Communicator&);
};

template<class T>
void Communicator::AlltoAll(std::vector<std::vector<T>>& snd,
                            std::vector<T>& rec) {
  const static long MAX_BUFFER = (long) std::numeric_limits<int>::max() / rankSize_;
  /* Check number of messages */
  std::vector<uint64_t> sndNumMsg(rankSize_);
  uint64_t maxSndNum = 0;
  for (int r = 0; r < rankSize_; ++r) {
    sndNumMsg[r] = snd[r].size();
    maxSndNum = std::max(maxSndNum, (uint64_t) snd[r].size());
  }

  std::vector<uint64_t> recNumMsg(rankSize_);
  MPI_Alltoall(&sndNumMsg[0],
               1,
               MPI_UINT64_T,
               &recNumMsg[0],
               1,
               MPI_UINT64_T,
               MPI_COMM_WORLD);

  uint64_t sumMsg = 0;
  for (int r = 0; r < rankSize_; ++r) {
    sumMsg += recNumMsg[r];
  }

  uint64_t maxNum = 0;
  MPI_Allreduce(&sumMsg, &maxNum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

  int numItr = maxNum / (MAX_BUFFER / sizeof(T)) + 1;
  if (maxNum > (MAX_BUFFER / sizeof(T))) {
    /* Multi AlltoAll */
    rec.resize(sumMsg);
    std::vector<uint64_t> sOffset(rankSize_);
    uint64_t rOffset = 0;
    for (int i = 0; i < numItr; ++i) {
      std::vector<int> sndMsgByte(rankSize_);
      for (int r = 0; r < rankSize_; ++r) {
        if ((sndNumMsg[r] - sOffset[r]) * sizeof(T) > MAX_BUFFER) {
          sndMsgByte[r] = (int) (MAX_BUFFER / sizeof(T)) * sizeof(T);
        } else {
          sndMsgByte[r] = (int) (sndNumMsg[r] - sOffset[r]) * sizeof(T);
        }
      }

      std::vector<int> recMsgByte(rankSize_);
      MPI_Alltoall(&sndMsgByte[0],
                   1,
                   MPI_INT,
                   &recMsgByte[0],
                   1,
                   MPI_INT,
                   MPI_COMM_WORLD);
      int rbyte = 0;
      std::vector<int> recMsgDips(rankSize_);
      for (int r = 0; r < rankSize_; ++r) {
        recMsgDips[r] = rbyte;
        rbyte += recMsgByte[r];
      }

      #pragma omp parallel for
      for (int r = 0; r < rankSize_; ++r) {
        MPI_Gatherv(&(snd[r][sOffset[r]]),
                    sndMsgByte[r],
                    MPI_BYTE,
                    &rec[rOffset],
                    &recMsgByte[0],
                    &recMsgDips[0],
                    MPI_BYTE,
                    r,
                    MPI_COMM_WORLD);
        sOffset[r] += ((long) sndMsgByte[r]) / sizeof(T);
      }
      rOffset += ((long) rbyte) / sizeof(T);
    }

    for (int r = 0; r < rankSize_; ++r) {
      std::vector<T>().swap(snd[r]);
    }

  } else {
    /* Single AlltoAll */
    std::vector<int> sndMsgByte(rankSize_);
    for (int r = 0; r < rankSize_; ++r) {
      sndMsgByte[r] = (int) sndNumMsg[r] * sizeof(T);
    }

    int rbyte = 0;
    std::vector<int> recMsgByte(rankSize_);
    std::vector<int> recMsgDips(rankSize_);
    for (int r = 0; r < rankSize_; ++r) {
      recMsgDips[r] = rbyte;
      recMsgByte[r] = (int) recNumMsg[r] * sizeof(T);
      rbyte += recMsgByte[r];
    }

    rec.resize(rbyte / sizeof(T));


    #pragma omp parallel for
    for (int r = 0; r < rankSize_; ++r) {
      MPI_Gatherv(&(snd[r][0]),
                  sndMsgByte[r],
                  MPI_BYTE,
                  &rec[0],
                  &recMsgByte[0],
                  &recMsgDips[0],
                  MPI_BYTE,
                  r,
                  MPI_COMM_WORLD);
    }

    for (int r = 0; r < rankSize_; ++r) {
      std::vector<T>().swap(snd[r]);
    }
  }
};

void Communicator::allGatherBool(bool requestRandom, bool*& recRequest) {
  MPI_Allgather(&requestRandom,
                1,
                MPI_CXX_BOOL,
                recRequest,
                1,
                MPI_CXX_BOOL,
                MPI_COMM_WORLD);
};

#endif //DISTRIBUTEDNE_COM_HPP
