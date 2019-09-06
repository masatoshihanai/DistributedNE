/*
 * main.cpp
 * 
 *  Copyright (c) 2019 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#include <iomanip>
#include <queue>
#include <unordered_set>
#include <unistd.h>

#include "boundary_q.hpp"
#include "allocator.hpp"
#include "dgraph.hpp"
#include "profiler.hpp"
#include "type.hpp"

int main(int argc, char** argv) {
  /* Init MPI configuration */
  int THR_SUPPORT = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &THR_SUPPORT);
  int rankSize, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &rankSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* Init Threads */
  int numThr = 1;
  if (getenv("OMP_NUM_THREADS")) {
    numThr = std::atoi(getenv("OMP_NUM_THREADS"));
  }
  omp_set_num_threads(numThr);

  /* Get Options */
  double expandRatio = 0.1;
  int opt;
  bool DRYRUN = false;
  bool waitForEnter = false;
  while ((opt = getopt(argc, argv, "hvde:b:ms:w")) != -1) {
    switch (opt) {
      case 'h':
        if (rank == 0) {
          std::cout << "==========================================================================" << std::endl;
          std::cout << "Distributed NE: A Scalable Parallel and Distributed Graph Edge Partitioner" << std::endl;
          std::cout << "==========================================================================\n" << std::endl;
          std::cout << "Usage: \n"
                       "  $ mpirun -n <# process> " << argv[0]
                    << " [option] <graph-file> <# partition> \n" << std::endl;
          std::cout << "Options:\n"
                    << "  -h: Help \n"
                    << "  -v: Output logs and performance details \n"
                    << "  -d: DryRun without outputting partitioned edges \n"
                    << "  -e <expansion ratio>: Expansion ratio (default=0.1) \n"
                    << "  -b <balance factor>: Balance factor. (default=1.01) \n"
                    //<< "      (max # of edges) / (average # of edges) becomes less than the factor. \n"
                    //<< "  -m: Multi objective mode. # of vertices is also balanced. \n"
                    << "  -s <seed>: Random seed (default=1) \n"
                    //<< "  -w: wait for attach process to debug \n"
                    << std::endl;
          std::cout << "Example to use: \n"
                    << "  $ mpirun -n 4 " << argv[0] << " ../data/Slashdot.edges 4\n"
                    << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
      case 'v':
        VERBOSE = true;
        break;
      case 'd':
        DRYRUN = true;
        break;
      case 'e':
        expandRatio = std::atof(optarg);
        break;
      case 'b':
        BALANCE_FACTOR = std::atof(optarg);
        break;
      case 'm':
        VERTEX_BALANCE = true;
        break;
      case 's':
        RSEED = std::atoi(optarg);
        break;
      case 'w':
        waitForEnter = true;
        break;
      default:
        break;
    }
  }

  if (waitForEnter) {
    if (rank == 0) {
      std::cout << "Enter to Start " << std::endl;
      while(getchar() != '\n');
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  /* Get file name and # partition */
  char* inputFile = argv[optind];

  /* Check whether options are vaild */
  if (rankSize != std::atoi(argv[optind + 1])) {
    if (rank == 0) {
      std::cerr << "!!! # of partitions must equal to # of MPI processes !!!" << std::endl;
      std::cerr << "# of partitions = " << argv[optind + 1] << "; # of MPI processes = " << rankSize << std::endl;
      exit(1);
    }
  }

  if (rankSize > P_BITS_SIZE*sizeof(uint64_t)*8) {
    if (rank == 0) {
      std::cerr << "!!! Partition size is too Large. Max Partition Size is !!!!!" << P_BITS_SIZE * sizeof(VertT) * 8
                << std::endl;
      std::cerr << "!!! Increase the size of the partition set (See P_BITS_SIZE in CMakeLists.txt).!!!!" << std::endl;
      exit(1);
    }
  }

  if (expandRatio > 1.0) {
    if (rank == 0) {
      std::cerr << "!!! Warning: Expand ratio must be less than 1.0 !!! " << std::endl;
      std::cerr << "!!! Warning: Run with Expand Ratio = 1.0 !!! " << std::endl;
    }
    expandRatio = 1.0;
  }

  if (rank == 0) std::cout << "Start " << argv[0] << std::endl;

  /* Init Distributed Graph */
  double initStart = MPI_Wtime();
  if (rank == 0) std::cout << "Initialization Start." << std::endl;
  MetaD::init(rankSize, rank);
  DistributedGraph graph;
  graph.init(inputFile, rank, rankSize, numThr);

  DistributedAllocator alloc;
  alloc.init(rank, rankSize, &graph);
  double initEnd = MPI_Wtime();
  if (rank == 0) std::cout << "Initialization Finish. Time (sec): " << (initEnd - initStart) << std::endl;

  if (rank == 0) {
    std::cout << "== Initialization Information ==" << std::endl;
    std::cout << "  Graph: " << inputFile << std::endl;
    int size = std::to_string(graph.numOrgEdges()).length();
    std::cout << "  # of Verts: " << std::setw(size) << graph.numVertices() << std::endl;
    std::cout << "  # of Edges: " << std::setw(size) << graph.numOrgEdges() << std::endl;
    std::cout << "  # of Partition: " << std::setw(size-5) << rankSize << std::endl;
    std::cout << "  Balance Factor: " << std::setw(size-5) << BALANCE_FACTOR << std::endl;
    if (VERBOSE) {
      std::cout << "    # of Threads   : " << numThr << std::endl;
      std::cout << "    Expansion Ratio: " << expandRatio << std::endl;
      //std::cout << "    Vertex Balance Mode: " << std::boolalpha << VERTEX_BALANCE << std::endl;
      std::cout << "    Random Seed    : " << RSEED << std::endl;
      std::cout << "    Vertex ID type : " << sizeof(VertT) * 8 << "-bit " << std::endl;
      std::cout << "    Max Partitions : " << MAX_NUM_PARTS << std::endl;
    }
    std::cout << "================================" << std::endl;
  }

  /* Set thresholds */
  EdgeNum localEdgeNumThreshold = (graph.numDEdges() / rankSize) * BALANCE_FACTOR;
  EdgeNum addEdgeNum = 0;

  std::unordered_set<VertT> replcaVert;
  VertT sumReplica;
  if (VERBOSE) {
    replcaVert.reserve(graph.numVertices() * 2.0 / rankSize);
  }

  /* Iteration */
  BoundaryQueue boundary;
  double start = MPI_Wtime();
  if (rank == 0) {
    std::cout << "!!!!!!!!!! Compute Partition !!!!!!!!!!!!!" << std::endl;
  }
  int iteration = 0;
  while (true) {
    /* Select expand vertices */
    VertT numExpand = expandRatio * boundary.size() < 1 ? 1 : expandRatio * boundary.size();

    if (boundary.size() == 0) {
      if (graph.allocDirectedEdges(rank) < localEdgeNumThreshold) {
        alloc.expandRandom();
      }
    } else {
      for (VertT i = 0; i < numExpand; ++i) {
        /* For edge balance */
        if (graph.allocDirectedEdges(rank) + addEdgeNum > localEdgeNumThreshold) {
          break;
        }

        /* Get from boundary */
        if (boundary.size() == 0) { break; }
        VertT x = boundary.smallestValuedID();
        boundary.popSmallest();
        alloc.bufferExpandVertex(x);
      }
    }

    /* Expand Boundary */
    alloc.expandBoundary();

    /* Update new boundaries */
    alloc.updateBoundary(boundary, replcaVert);

    /* Check termination */
    EdgeNum sumRestEdges = graph.sumRestDirectedEdges();
    if (sumRestEdges == 0 || sumRestEdges > graph.numEdges() * 2) {
      break;
    }

    /* Output Info */
    if (VERBOSE) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0) {
        graph.debugInfo(iteration);
      }
    } else {
      if (rank == 0) {
        std::cout << "\r # of Iteration: " << iteration
                  << ". # of Processed Edges: "
                  << std::setw(std::to_string(graph.numEdges()).length())
                  << std::right << graph.numEdges() - (sumRestEdges / 2)
                  << " / " << graph.numEdges() << std::flush;
      }
    }

    ++iteration;
  } /* while */

  if (rank == 0) {
    if (!VERBOSE) std::cout << std::endl;
    std::cout << "!!!!!!!!!! Finish Partition !!!!!!!!!!!!!" << std::endl;
  }
  double end = MPI_Wtime();

  /* Analyze Result */
  EdgeNum numReplica = alloc.numReplica();
  EdgeNum numInitReplica = 0;
  VertT maxReplica;
  VertT localReplica;
  double balanceScore;
  if (VERBOSE) {
    numInitReplica = alloc.numInitReplica();
    localReplica = alloc.numLocalReplica();
    MPI_Allreduce(&localReplica, &sumReplica, 1, MPI_VERT_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localReplica, &maxReplica, 1, MPI_VERT_T, MPI_MAX, MPI_COMM_WORLD);
    balanceScore = alloc.balanceQuality();
  }
  if (rank == 0) {
    std::cout << "======== Result Information ========" << std::endl;
    std::cout << "  Replication factor: " << ((double) numReplica / (double) graph.numUniqueVertices()) << std::endl;
    std::cout << "  Compute time (sec): " << (end - start) << std::endl;
    if (VERBOSE) {
      std::cout << "    # of Replication : " << numReplica << std::endl;
      std::cout << "    # of Iterations  : " << iteration << std::endl;
      std::cout << "    Edge Balance     : " << balanceScore << std::endl;
      std::cout << "    Vertex Balance   : " << ((double) maxReplica * (double) rankSize / (double) sumReplica) << std::endl;
      std::cout << "    Init. RF(2D-hash): " << ((double) numInitReplica / (double) graph.numUniqueVertices()) << std::endl;
    }
    std::cout << "====================================" << std::endl;
  }

  /* Output Result */
  if (!DRYRUN) alloc.outputResult(inputFile);

  /* Check result */
  if (!alloc.checkResult()) {
    if (rank == 0) { std::cerr << "Fail to compute " << std::endl; }
    MPI_Finalize();
    return 1;
  };

  MPI_Finalize();
  return 0;
}