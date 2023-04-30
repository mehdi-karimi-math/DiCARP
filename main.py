# This code runs DiCARP (Distributed Consensus Algorithm with Radial Partitioning)
# to a problem given as a Pyomo data file.
#
# Copyright (c) 2023, by 
# Mehdi Karimi

from dist_opt_main_2 import *
import pickle


problem = "case9"  # The name of the problem
rho = 200 # Parameter of the DiCARP algorithm 
tol = .0001 # Tolerance of the algorithm for stopping criteria
max_iter = 2000 # maximum number of iterations to run the algorihtm 


result = DCA_algorithm(problem, tol , rho, max_iter)

print("The objective value of the DiCA algorithm = {} ".format(result["P_da"]))
print("The objective value of the Ipopt = {} ".format(result["P_ipm"]))
print("GAP between the DiCA and Ipopt = {} ".format(result["GAP"]))

