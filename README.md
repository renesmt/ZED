# ZED
## Before use
The code depends on Blossom V(https://pub.ista.ac.at/~vnk/papers/blossom5.pdf) to solve the ground states of planar graphs. One should try to compile a MATLAB version of Blossom V algorithm. 
The Blossom V should be implemented in this way:
``
blossom($N_{node}$,$N_{edge}$)
``
This code also includes a random regular graph generator.(https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator)
You can also contact the author(smutian@wustl.edu) for a compiled UNIX matlab function.
## Solve single instance
  1. Initialize the bond configuration generator by `GGenerator.m`.
  2. Initialize the solver by `SSolver.m`.
  3. Solve the bond configuartion $J_{ij}$ by initialized sovler.
## Solve multiple instances and compute the zero-energy droplets
Just set the parameters in `CritDroplet.m`, including lattice type, system size, number of bond configurations and additional parameters that is necessary for some specific kinds of lattices such as random regular graph.
