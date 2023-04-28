# DiCARP 1.0: Distributed Consensus Algorithm with Radial Partitioning  

This code implements a distributed consensus algorithm (DCA) for optimial power flow (OPF) problems using a special *Radial Partitioning*. The code uses the python modelling package for optimization Pyomo, using the solver Ipopt. 

Some case instances from the MATPOWER library that are trasformed into Pyomo data format are given in the `Data` folder. There is MATLAB m file `opf_data_file.m` that transforms MATPOWER case files into the Pyomo DAT file for DiCARP. 

# Using DiCARP

For using DiCARP, Python, DiCARP, and the required modules must be installed. This code is using Ipopt as the underlying solver, which must be installed. The easiest way to do is dwnloading `ipopt.exe` from this [site](https://www.coin-or.org/download/binary/Ipopt/) and putting it in the required forlder. 
