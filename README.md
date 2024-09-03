# Heuristics for Adjacency Graph Packing Considering Intergenic Regions

This repository has the implementation of two heuristics for the Adjacency Graph Packing problem considering intergenic regions. One (rand) tries to solve the problem by randomly generating cycle packings of the Adjacency Graph, constructed with best (there is an alternative decomposition implemented that randomly selects edges for the decomposition). The other is based on the Genetic Algorithm metaheuristic (ga).

We also have the implementation of an 4/3-approximation algorithm for the dcj distance problem with intergenic regions for genomes with a single copy of each gene [[1]](#1).

## Usage

Compile the code by running `make` and see the running options with `./dec --help`.

## Database

The folder `database` has files with simulated genomes used to test the heuristics.

## References

<a id="1">[1]</a> 
Guillaume Fertin, Géraldine Jean and Eric Tannier. Algorithms for computing the double cut and join distance on both gene order and intergenic sizes. Algorithms for Molecular Biology 2017;12(1):16.
