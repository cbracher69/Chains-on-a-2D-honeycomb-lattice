Chains on a 2D honeycomb lattice
================================

This is a computing project started as a tutorial with two students at Bard College, Yan Chu and Yuexi Ma.  Although the challenge has interesting applications in math (graphs on a honeycomb lattice) and statistical physics (where the chains are a crude representation of bonds in a polymer molecule), our main motivation is to explore high-performance computing on a desktop computer, and learning how to use parallel computation on a GPU (NVidia CUDA):  How far can we go exploring a problem with exponentially increasing complexity?

### Introduction:  Chains on a honeycomb lattice

A [*honeycomb lattice*](http://en.wikipedia.org/wiki/Hexagonal_tiling) in two dimensions is the familiar tiling of the plane with regular hexagons.  In this lattice, every vertex is connected by edges to its three nearest neighbors.

A *chain* with N segments on this lattice consists of N such edges.  It is created by N sequential moves from a vertex to a neighboring vertex, such that the chain never immediately reverses direction (does not jump back to the vertex it came from in the immediately preceding step).  This leaves two possible directions at each vertex - the chain can continue turning either "left" or "right."  In a chain with N segments, N-2 such decisions are made.  Therefore, the number of different chains with N segments is 2^(N-2), and they can be encoded in binary form:

    011010111000 ... 1010
    
where 0 stands for "turn left" and 1 for "turn right."

In this project, we are only interested in non-overlapping ([self-avoiding](http://en.wikipedia.org/wiki/Self-avoiding_walk)) open chains and polygons.  An open chain is simply one that does not fold back on itself - every vertex is visited only once.  A *polygon* on a lattice is a chain that returns to its starting point.  In a self-avoiding polygon, the chain does not cross itself.  Finding all self-avoiding chains and polygons of a given length is considered an NP-hard problem.

##### Embedding of the honeycomb lattice

Mathematically, it can be insightful not to consider the honeycomb lattice as a fundamental pattern, but rather as a subspace of a higher-dimensional lattice with simpler structure.  The honeycomb lattice naturally arises in the three-dimensional simple cubic grid when it is "sliced" perpendicular to the space diagonal of the cube (in the (111) direction).  Using the 3D representation, it is easy to find a natural distance metric in the honeycomb grid (it equals the L1 metric, or ["Manhattan distance,"](http://en.wikipedia.org/wiki/Manhattan_distance) in the cubical lattice.

##### Outlook:  Higher-dimensional honeycomb grids

Using the embedding method, it becomes possible to define "honeycomb grids" in three- and higher-dimensional spaces, where each vertex has d+1 equlvalent nearest neighbors (d is the dimension of the space).  They are represented as hyperplanes in the (d+1)-dimensional simple cubic lattice (cartesian coordinate grid), and have the same underlying metric.  For instance, the 3D equivalent of the honeycomb lattice is the [*diamond lattice*](http://en.wikipedia.org/wiki/Diamond_lattice) commonly encountered in semiconductor physics.

### Contents of this repository

*   2D Chain Deterministic CPU.cpp

    This is the latest version (3.6) of the program designed to run on the CPU only.  The program takes the number of segments as an input, and then examines the whole set of possible chains.  It will collect and display summary information about the number of self-avoiding chains and polygons, the number of unique such polygons (i.e., polygons that are not congruent to each other - that cannot be mapped by a simple rotation or translation of the polygon on the grid), and their distribution among symmetry classes: Polygons can be sorted out according to their rotational symmetry (1-, 2-, 3-, or 6-fold rotation axis), and their symmetry properties (polygons are either mirror images of themselves, or they have no such symmetry).
    
    Version 3.6 uses a number of tricks to keep the numerical expense down, such as reusing parts of chains, skipping chains known to start with an overlapping sequence, and using the distance metric to minimize overlap checks, resulting in only about 10% overhead - most overlapping chains and polygons are eliminated before they are ever tested.  The program has been used to find all self-avoiding chains and polygons up to n=42 segments (see info sheet).
    
*   2D Chain Deterministic CUDA.cu

    This is a parallelized version (1.1) that uses a NVidia CUDA-enabled GPU to perform the calculation. (Note:  This requires installing CUDA (free download), available at the [NVidia website](https://developer.nvidia.com/about-cuda).)  The code employs parallel threads on the GPU to perform the overlap check for each chain, without any of the more sophisticated algorithms employed in the CPU version to eliminate "bad" chains without check.  For a system with a Tesla C2075 GPU and Core i7 CPU, parallelization led to a speed-up of a factor 20 compared to a one-thread CPU version of the same code. (The more sophisticated algorithm used in the current CPU code cuts that advantage down by a factor of 10.)
    
*   Honeycomb Chain Data
   
    Summary information about the number and distribution of self-avoiding chains and polygons on the 2D honeycomb lattice, for lengths up to 42 segments.

### Future directions

*   Combine the advantages of the CPU approach (chain elimination) and the GPU approach (brute force analysis of chains for overlaps).  Develop a divide-and-conquer model that splits the chains in half, collects "viable" segments using the CPU, and links pairs of them back together to a chain of the original length, with the results to be analyzed in parallel on the GPU.  The method should harvest the speed advantages of parallel computation, without increasing the overhead dramatically, leading to an expected 10-fold acceleration for the configuration noted above.
*   Improve memory handling to avoid caching/paging issues for large chain sizes.
*   Increase the length of tractable chains.  The target is n=50.  (Currently, the longest chains analyzed have n=42.)



