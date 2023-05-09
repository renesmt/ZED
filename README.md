# ZED
## Before use

The code depends on Blossom V(https://pub.ista.ac.at/~vnk/papers/blossom5.pdf) to solve the ground states of planar graphs. One should try to compile a MATLAB version of Blossom V algorithm. 
The Blossom V should be implemented in this way:


`blossom(number of nodes, number of edges, weights)`: for complete graphs, the weights is sorted from left to right, from top to down in the adjacency matrix of the graph. For example,
$\left(\begin{array}{cc} 
0.8944272 & 0.4472136\\
-0.4472136 & -0.8944272
\end{array}\right)$

`blossom2(number of nodes, number of edges, weights, edge pairs)`

This code also includes a random regular graph generator.(https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator)

You can also contact the author(smutian@wustl.edu) for a compiled UNIX matlab function.
## Solve single instance
  1. Initialize the bond configuration generator by `GGenerator.m`.
  2. Initialize the solver by `SSolver.m`.
  3. Solve the bond configuartion $J_{ij}$ by initialized sovler.
## Solve multiple instances and compute the zero-energy droplets
Just set the parameters in `CritDroplet.m`, including lattice type, system size, number of bond configurations and additional parameters that is necessary for some specific kinds of lattices such as random regular graph.
