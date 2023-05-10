# ZED
This package of code could be used to compute the zero-energy droplets of spin glasses, and compute some basic statistical properties of these droplets, such as the average droplet size or the distribution of critical threshold values.
## Before use
The author is keeping improving the readability of this code: now it is terrible :(. Please do not hesitate to contact the author (smutian@wustl.edu) if you are confused by anything.

The code depends on Blossom V (https://pub.ista.ac.at/~vnk/papers/blossom5.pdf) to solve the ground states of planar graphs. One should try to compile a MATLAB version of Blossom V algorithm. 
The Blossom V should be implemented in this way:

`blossom(number of nodes, number of edges, weights)`: for complete graphs, or fully-connected models. The weights are sorted from left to right, from top to down in the adjacency matrix of the graph. For example, for such kind of adjacency matrix:
```math
\left(\begin{array}{cccc} 
0 & a_1 & a_2 & a_3\\ 
a_1 & 0 & a_4 & a_5\\ 
a_2 & a_4 & 0 & a_6\\ 
a_3 & a_5 & a_6 & 0
\end{array}\right)
```
The `weights` should be  in the order of $[a_1,a_2,\cdots,a_6]$.

`blossom2(number of nodes, number of edges, weights, edge pairs)`: for general graphs. the `edge pairs` is in the form of $N_{\rm edge}\times 2$ matrix. For each row, the two elements correspond to the node index, from $0$ to $N_{\rm nodes}-1$.

An exact solver is also necessary. Here, Gurobi(www.gurobi.com) is adopted. Other exact solvers are also okay, as long as they are implemented in this format:  
`ground_state = GurobiComputeGround(JJJ)`,   
where `JJJ` is the adjacency matrix of the bond configuration $J_{ij}$.

You can also contact the author(smutian@wustl.edu) for a compiled UNIX matlab function.

This code also includes a random regular graph generator(https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator), which is already included in `GGenerator.m`.

## Solve single instance
  1. Initialize the bond configuration generator by `GGenerator.m`. For example,     
  `generator = GGenerator("square");`
  2. Generate the bond configuartion. For example, generate a $4\times 5$ (number of plaquettes) system with seed $12321$ by:   
   `[JJJ, i1, i2, nodes, c1, c2] = generator(4,5,12321);`   
  3. Initialize the solver by `SSolver.m`. For example,    
   `solver = SSolver("square");`
  4. Solve the bond configuartion $J_{ij}$ by initialized solver.   
   `ground_state = solver(nodes);`   
  For non-planar graph:   
  `ground_state = solver(JJJ);`
## Solve multiple instances and compute the zero-energy droplets
Just set the parameters in `CritDroplet.m`, including lattice type, system size, number of bond configurations, random number seed and additional parameters that is necessary for some specific kinds of lattices such as random regular graph. For example, to compute $10^5$ square instances with $L=128$ with random number seed $12321$:     
`CritDroplet("square",128,1e5,12321,0)`.
