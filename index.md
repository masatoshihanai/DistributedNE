![DNE](pics/DNE_logo.png)

__*Distributed Neighbor Expansion (Distributed NE)*__ : A scalable parallel and distributed graph partitionor for high-quality edge partitioning, which has several key important features:

- __Scalablable to trillion-edge graphs__
- __Provide high-quality edge partitionins__

- __Theoretical gurantee of partitioning quality__

The algorithm is based on  _parallel expansion_, where the edge parts are greedy expanded in parallel in such way that the increase of vertex-cuts in each part becomes locally minimal as the figure.

![expansion](pics/ParallelExpansion.png)    

### Reference.

M. Hanai, et.al. *"Distributed Edge Partitioning for Trillion-edge Graphs"* (PVLDB 2020, [Paper](https://arxiv.org/abs/1908.05855))



## Quick Start

 (Require `MPI` and `C++11`)

```bash
$ git clone git@github.com:masatoshihanai/DistributedNE.git
$ cd DistributedNE; mkdir build; cd build;
$ 
```

#### Tested Environments

- MPI:  `Intel MPI 1.4.4``OpenMPI`
- Compiler: `g++xxx`
- OS: `Ubuntsu`xxxxx 

### How to compile

```bash
$

```



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