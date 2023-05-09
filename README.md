# ZED
## Before use

The code depends on Blossom V (https://pub.ista.ac.at/~vnk/papers/blossom5.pdf) to solve the ground states of planar graphs. One should try to compile a MATLAB version of Blossom V algorithm. 
The Blossom V should be implemented in this way:


`blossom(number of nodes, number of edges, weights)`: for complete graphs, the weights is sorted from left to right, from top to down in the adjacency matrix of the graph. For example, for such kind of adjacency matrix:
```math
\left(\begin{array}{cccc} 
0 & a_1 & a_2 & a_3\\ 
a_1 & 0 & a_4 & a_5\\ 
a_2 & a_4 & 0 & a_6\\ 
a_3 & a_5 & a_6 & 0
\end{array}\right)
```
The value of weights should be taken in the order of $a_1,a_2,\cdots,a_6$.

`blossom2(number of nodes, number of edges, weights, edge pairs)`: for general graph. the `edge pairs` is in the form of $N_{\rm edge}\times 2$ matrix. For each row, the two elements correspond to the node index, from $0$ to $N_{\rm nodes}-1$.

This code also includes a random regular graph generator.(https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator)

You can also contact the author(smutian@wustl.edu) for a compiled UNIX matlab function.
## Solve single instance
  1. Initialize the bond configuration generator by `GGenerator.m`.
  2. Initialize the solver by `SSolver.m`.
  3. Solve the bond configuartion $J_{ij}$ by initialized sovler.
## Solve multiple instances and compute the zero-energy droplets
Just set the parameters in `CritDroplet.m`, including lattice type, system size, number of bond configurations and additional parameters that is necessary for some specific kinds of lattices such as random regular graph.
