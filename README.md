# DiCARP 1.0: Distributed Consensus Algorithm with Radial Partitioning  

This code implements a distributed consensus algorithm (DCA) for optimial power flow (OPF) problems using a special *Radial Partitioning* and adaptive *spectral parameter selection*. The code uses the python modelling package for optimization Pyomo, using the solver Ipopt. 

Some case instances from the MATPOWER library that are trasformed into Pyomo data format are given in the `Data` folder. There is MATLAB m file `opf_data_file.m` in the `Data` folder that transforms MATPOWER case files into the Pyomo DAT file for DiCARP. 

# Using DiCARP

For using DiCARP, Python, Pyomo, and the required modules must be installed. This code is using Ipopt as the underlying solver, which must be installed as well. The easiest way to do is downloading `ipopt.exe` from this [site](https://www.coin-or.org/download/binary/Ipopt/) and putting it in the required forlder. 

The `main` file has a code to run DiCARP for a given problem instance. The parameters to set are:

+ `problem`: A string with the name of the problem, for example `"case9"`.
+ `rho`: The parameter of the distributed algorithm. 
+ `tol`: The tolerance for the stopping criteria. 
+ `mat_iter`: Maximum number of iterations set for the distributed algorithm.  

The main function is `DCA_algorithm` which returns a python dictionary with the results and statistics of the distributed algorithm. 

# Citation

DiCARP is based on the following manuscript:

Mehdi Karimi. **Radial Partitioning with Spectral Penalty Parameter Selection in Distributed Optimization for Power Systems**, [https://arxiv.org/abs/2305.01032](https://arxiv.org/abs/2305.01032), 2024. 
