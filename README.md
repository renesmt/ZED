# ZED
Usage:
## Solve single instance
  1. Initialize the bond configuration generator by `GGenerator.m`.
  2. Initialize the solver by `SSolver.m`.
  3. Solve the bond configuartion $J_{ij}$ by initialized sovler.
## Solve multiple instances and compute the zero-energy droplets
Just set the parameters in `CritDroplet.m`, including lattice type, system size, number of bond configurations and additional parameters that is necessary for some specific kinds of lattices such as random regular graph.
