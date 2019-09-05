![DNE](pics/DNE_logo.png)

__*Distributed Neighbor Expansion (Distributed NE)*__ : A scalable parallel and distributed graph partitioner for high-quality edge partitioning, which has several key features:

- __Low-memory and high-performance implementation__
- __Provide high-quality edge partitionins__
- __Scalablable to trillion-edge graphs on 200+ machines__

The algorithm is based on  _parallel expansion_, where the edge parts are greedily expanded in parallel in such a way that the increase of vertex-cuts in each part becomes the local minimum. 

![expansion](pics/ParallelExpansion.png)

*Reference*    

> M. Hanai, et.al. *"Distributed Edge Partitioning for Trillion-edge Graphs"* (PVLDB 2020, [Paper](https://arxiv.org/pdf/1908.05855.pdf))

## Quick Start {#q-start} 

 (Require `MPI` and `C/C++` )

```bash
$ git clone git@github.com:masatoshihanai/DistributedNE.git
$ cd DistributedNE; mkdir build; cd build
$ cmake ..; make
$ mpirun -n 4 ./DistributedNE ../data/LJ.edges 4
```

## Outline of the Page

1. [__How to compile__](#compile)
2. [__How to run__](#run)
3. [__Graph data format__](#data)
4. [__Related tools__](#related)
5. [__Acknowledgement / Contact__](#ack)

---

### 1. How to compile {#compile}

First, prepare the standard development tools for C/C++ on MPI

__Requirements__: `MPI` `C++`  `Make` `CMake` `Git`

- __MPI__ : [OpenMPI](https://www.open-mpi.org/) , [MPICH](https://www.mpich.org/), or etc.
- __C++ compiler__: [GCC](https://gcc.gnu.org/install/), [IntelC++](https://software.intel.com/en-us/c-compilers), [Clang/LLVM](https://clang.llvm.org/index.html), or etc.
- [Make](https://www.gnu.org/software/make/), [CMake](https://cmake.org/), [Git](https://git-scm.com/)

Next, get the code from github

```bash
### get the code from git repository
$ git clone git@github.com:masatoshihanai/DistributedNE.git
```

Configure and compile

```bash
### make new build directory
$ cd DistributedNE
$ mkdir build
$ cd build
```

```bash
### configure and compile
$ cmake ..
$ make
```

We have successed to compile and run the program in these environments:

##### Tested Environments

- MPI:  `OpenMPI 2.1` `OpenMPI 4.0` `Intel MPI 5.1`
- Compiler: `GCC 4.9` `GCC 5.4` `GCC 7.4`
- OS: `Red Hat Server 6.9` `Ubuntu 18.04`

### 2. How to run {#run}

```bash
$
```



### 3. Graph data format {#data}

Distributed NE supports the de-facto standard __edge-list format__ by [SNAP](https://snap.stanford.edu/data/index.html),  where each line represents one directed edge (source id and destination id).

For example, <img src="pics/Graph.png" alt="Graph" style="zoom:50%;" />is  represented as follows:

 ```bash
# src dst
0    1
1    0
1    2
1    3
2    1
3    2
 ```



### 4. Related tools {#related}

>_PowerGraph_ - Edge-partitioned distributed graph processing system ([Code](https://github.com/jegonzal/PowerGraph))
>
>_GraphX_ - Spark-based Edge-partitioned distributed graph processing system ([Code](https://spark.apache.org/graphx/))
>
>_NE_ - Graph edge partitioner based on sequential neighbor expansion ([Code](https://www.kdd.org/kdd2017/papers/view/graph-edge-partitioning-via-neighborhood-heuristic)) 
>
>_METIS_ - Standart vertex patitoining tools ([Code](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview))

### 5. Acknowledgement {#ack}

![SUStech](pics/sustech.png)  ![IBM](pics/IBM.png) ![NTU](pics/NTU.png) 

###### Contact: 

_Masatoshi Hanai_ (mhanai at acm.org)