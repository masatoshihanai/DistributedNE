![DNE](pics/DNE_logo.png)

__*Distributed Neighbor Expansion (Distributed NE)*__ : A scalable parallel and distributed graph partitioner for high-quality edge partitioning, which has several key features:

- __Low-memory and high-performance implementation__
- __Provide high-quality edge partitionins__

- __Scalablable to trillion-edge graphs on 200+ Machines__

The algorithm is based on  _parallel expansion_, where the edge parts are greedily expanded in parallel in such way that the increase of vertex-cuts in each part becomes locally minimal as the figure. 

![expansion](pics/ParallelExpansion.png)    

### Reference.

M. Hanai, et.al. *"Distributed Edge Partitioning for Trillion-edge Graphs"* (PVLDB 2020, [Paper](https://arxiv.org/pdf/1908.05855.pdf))

## Quick Start

 (Require `MPI` and `C/C++` )

```bash
$ git clone git@github.com:masatoshihanai/DistributedNE.git
$ cd DistributedNE; mkdir build; cd build
$ cmake ..; make
$ mpirun -n 4 ./DistributedNE ../data/LJ.edges 4
```

### How to compile

First, prepare the standard development tools for C/C++ on MPI

__Requirements__: `MPI` `C++`  `Make` `CMake` `Git`

- __MPI__ : [OpenMPI](https://www.open-mpi.org/) , [MPICH](https://www.mpich.org/)
- __C++ compiler__: [GCC](https://gcc.gnu.org/install/), [IntelC++](https://software.intel.com/en-us/c-compilers), [Clang/LLVM](https://clang.llvm.org/index.html)
- [Make](https://www.gnu.org/software/make/), [CMake](https://cmake.org/) [Git](https://git-scm.com/)

Next, get the code from github.

```bash
### get the code from git repository
$ git clone git@github.com:masatoshihanai/DistributedNE.git
```

Build via CMake

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

We have successed to compile & run the program in these environments:

##### Tested Environments

- MPI:  `Intel MPI 5.1.2`
- Compiler: `GCC 5.4`
- OS: `Red Hat Server 6.9`

### How to Run Tools

```bash
$
```



### Supported data format

``` bash
$
```



#### Aknowledgement

![SUStech](pics/sustech.png)  ![IBM](pics/IBM.png) ![NTU](pics/NTU.png) 

###### Contact: 

_Masatoshi Hanai_ (mhanai at acm.org)