![DNE](pics/DNE_logo.png)

__*Distributed Neighbor Expansion (Distributed NE)*__ : A scalable parallel and distributed graph partitioner for high-quality edge partitioning, which has several key features:

- __Low-memory and high-performance implementation__
- __Provide high-quality edge partitionins__
- __Scalable to trillion-edge graphs on 200+ machines__ (1 trillion = 1,000,000,000,000)

The algorithm is based on  _parallel expansion_, where the edge parts are greedily expanded in parallel in such a way that the increase of vertex-cuts in each part becomes the local minimum. 

![expansion](pics/ParallelExpansion.png)

*Reference*    

> M. Hanai, et.al. *"Distributed Edge Partitioning for Trillion-edge Graphs"* (PVLDB 2019, [Paper](https://arxiv.org/pdf/1908.05855.pdf))

## Quick Start {#q-start} 

```bash
$ git clone git@github.com:masatoshihanai/DistributedNE.git
$ cd DistributedNE; mkdir build; cd build
$ cmake ..; make ## Require MPI 
$ mpirun -n 4 ./DistributedNE ../data/Slashdot.edges 4
```

## Outline of the Page

1. [__Graph data format__](#data)
2. [__How to compile__](#compile)
3. [__How to run__](#run)
4. [__Related tools__](#related)
5. [__Acknowledgement / Contact__](#ack)

---

### 1. Graph data format {#data}

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

### 2. How to compile {#compile}

First, prepare the standard development tools for C/C++ on MPI

__Requirements__: `MPI` `C/C++ compiler`  `Make` `CMake` `Git`

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

- MPI:  `OpenMPI 2.1` `OpenMPI 4.0` `IntelMPI 5.1`
- Compiler: `GCC 4.9` `GCC 5.4` `GCC 7.4`
- OS: `RedHat 6.9` `Ubuntu 18.04`

### 3. How to run {#run}

```bash
Usage: 
  $ mpirun -n <\# process> ./DistributedNE [option] <graph-file> <\# partition> 
```

The program uses `MPI` to execute on the distribtued environment. Thus, for example, if you want to partition `../data/Slashdot.edges` into 4 parts on 4 MPI processes,  the command is as follows:

```bash
$ mpirun -n 4 ./DistributedNE ../data/Slashdot.edges 4
```

Then, the program outputs the partitioned graph to `../data/Slashdot.edges.4.pedges` , where each line represents one edge and its partition ID including  \<Source ID\>,  \<Destination ID\>, \<Partition ID\> like that: 

```bash
$ head ../data/Slashdot.edges.4.pedges
# SrcID  DstID  PartitionID 
16808 3450 0
3450 16808 0
16808 4986 0
4986 16808 0
16808 5102 0
5102 16808 0
16808 5490 0
16808 5568 0
5568 16808 0
```

There are several useful options as you can see by running with  `-h` option

```bash
$ ./DistributedNE -h
==========================================================================
Distributed NE: A Scalable Parallel and Distributed Graph Edge Partitioner
==========================================================================

Usage: 
  $ mpirun -n <\# process> ./DistributedNE [option] <graph-file> <\# partition> 

Options:
  -h: Help 
  -v: Output logs and performance details 
  -d: DryRun without outputting partitioned edges 
  -e <expansion ratio>: Expansion ratio (default=0.1) 
  -b <balance factor>: Balance factor. (default=1.01) 
  -s <seed>: Random seed (default=1) 

Example to use: 
  $ mpirun -n 4 ./DistributedNE ../data/Slashdot.edges 4
  
```

### 4. Related tools {#related}

>_PowerGraph_ - Edge-partitioned distributed graph processing system ([Code](https://github.com/jegonzal/PowerGraph))
>
>_GraphX_ - Spark-based edge-partitioned distributed graph processing system ([Code](https://spark.apache.org/graphx/))
>
>_NE_ - Graph edge partitioner based on sequential neighbor expansion ([Code](https://www.kdd.org/kdd2017/papers/view/graph-edge-partitioning-via-neighborhood-heuristic)) 
>
>_METIS_ - Standard vertex partitioner ([Code](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview))

### 5. Acknowledgement {#ack}

![SUStech](pics/sustech.png)  ![IBM](pics/IBM.png) ![NTU](pics/NTU.png) 

#### Contact 

_Masatoshi Hanai_ (mhanai at acm.org)
