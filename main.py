# This code runs DiCARP (Distributed Consensus Algorithm with Radial Partitioning)
# to a problem given as a Pyomo data file.
#
# Copyright (c) 2024, by 
# Mehdi Karimi

from dist_opt_main_2 import *
import pickle
from radial_decom import *
folder = "Data/"

num_nodes = 9  # Number of nodes in the network
problem = "case9" # The name of the problem
rhov = 10000 # Parameter of the DiCARP algorithm 
rhop = 1000
if num_nodes  >= 1000:
    rhov = 100 # Parameter of the DiCARP algorithm 
    rhop = 10
tol = .0001 # Tolerance of the algorithm for stopping criteria
max_iter = 4000 # maximum number of iterations to run the algorihtm 
adaptive = True
component = False


result = DCA_algorithm(problem, tol , rhov, rhop, max_iter, adaptive, component)
print("The objective value of the DiCA algorithm = {} ".format(result["P_da"]))
print("The objective value of the Ipopt = {} ".format(result["P_ipm"]))
print("GAP between the DiCA and Ipopt = {} ".format(result["GAP"]))


# if component:
#     regions = create_regions_2_component(folder+problem)
# else:
#     regions = create_regions_2(folder+problem)

# # print(regions)
# num_com = 0
# for m in range(len(regions)):
#     num_com = num_com + len(regions[m+1]['vi_neigh'])
#     num_com = num_com + len(regions[m+1]['ei_neigh'])
#     num_com = num_com + len(regions[m+1]['gi_neigh'])
# print(num_com)