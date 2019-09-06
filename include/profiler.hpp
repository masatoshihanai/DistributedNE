/*
 * profiler.hpp 
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef DISTRIBUTEDNE_PROFILER_HPP
#define DISTRIBUTEDNE_PROFILER_HPP

struct Profiler {
  static double Start;
  static double Expand;
  static double Alloc;
  static double ComRanks;
  static double AllocInBoundary;
  static double SyncScore;
  static double ComScore;
  static double ComEdges;
  static double UpdateInf;
  static double ClearBuf;
  static double UpdateBoundary;
  static double Output;

  static void ProfileInit() {
    std::cout << "| Total Time | Expand | Alc 1hop | Com Parti. | Alc 2hop | Com Score | Com Edges | Updt Bndry |";
  }
  static void Profile() {
    int size = std::max((int) std::to_string(Expand + Alloc + ComRanks + AllocInBoundary + SyncScore + ComScore + ComEdges + UpdateInf + ClearBuf + UpdateBoundary).length(), 12);
    std::cout << "|" << std::setw(size) << Expand + Alloc + ComRanks + AllocInBoundary + SyncScore + ComScore + ComEdges + UpdateInf + ClearBuf + UpdateBoundary;

    size = std::max(8, (int)std::to_string(Expand).length());
    std::cout << "|" << std::setw(size) << Expand;

    size = std::max(10, (int)std::to_string(Alloc).length());
    std::cout << "|" << std::setw(size) << Alloc;

    size = std::max(12, (int)std::to_string(ComRanks).length());
    std::cout << "|" << std::setw(size) << ComRanks;

    size = std::max(10, (int)std::to_string(AllocInBoundary).length());
    std::cout << "|" << std::setw(size) << AllocInBoundary;

    size = std::max(11, (int)std::to_string(ComScore).length());
    std::cout << "|" << std::setw(size) << ComScore;

    size = std::max(11, (int)std::to_string(ComEdges).length());
    std::cout << "|" << std::setw(size) << ComEdges;

    size = std::max(12, (int)std::to_string(UpdateBoundary).length());
    std::cout << "|" << std::setw(size) << UpdateBoundary;

    std::cout << "|";
  };
 private:
  Profiler(){};
  ~Profiler(){};
};

double Profiler::Start = 0.0;
double Profiler::Expand = 0.0;
double Profiler::Alloc = 0.0;
double Profiler::ComRanks = 0.0;
double Profiler::AllocInBoundary = 0.0;
double Profiler::SyncScore = 0.0;
double Profiler::ComScore = 0.0;
double Profiler::ComEdges = 0.0;
double Profiler::UpdateInf = 0.0;
double Profiler::ClearBuf = 0.0;
double Profiler::UpdateBoundary = 0.0;
double Profiler::Output = 0.0;

#endif //DISTRIBUTEDNE_PROFILER_HPP
