# The main method applies Radial Partitioning to a given connected graph. 
# The graph must be given as a Pyomo DAT file.
#
# Copyright (c) 2023, by 
# Mehdi Karimi

from OPF_opt import *
from collections import defaultdict
import random


def create_regions(problem):
    instance=mo.create_instance(problem+".dat")

    ver = list(instance.vi)
    edg = list(instance.ei)
    gen = list(instance.gi)

    fr = list(instance.fr.values())
    to = list(instance.to.values())

    adj_l = defaultdict(list)


    for i in range(len(edg)):
        adj_l[fr[i]].append(to[i])
        adj_l[to[i]].append(fr[i])

    num_reg = 1
    reg_n = {}
    ver_c = ver[:]
    ver_f = [0 for _ in range(len(ver_c))]

    while len(ver_c):
        root = ver_c[0] # random.choice(ver_c) #
        ver_f[root-1] = num_reg
        se = set([root])
        se_ext = set()
        ver_c.remove(root)
        stack=[(i,root) for i in adj_l[root] if ver_f[i-1]==0]
        while len(stack)>0:
            cur, par = stack[0]
            stack = stack[1:]
            temp = [i for i in adj_l[cur]]
            temp.remove(par)
            temp2 = temp[:]
            if len(se.intersection(set(temp))) == 0:
                se.add(cur)
                ver_c.remove(cur)
                ver_f[cur-1] = num_reg
                stack = stack + [(i,cur) for i in temp if ver_f[i-1]==0]
        reg_n[num_reg] = se
        num_reg = num_reg+1

    num_reg = num_reg-1
    regions = {}
    for i in range(num_reg):
        regions[i+1] = {"vi_region_all" : reg_n[i+1].copy(), "vi_region" : reg_n[i+1].copy(),
        "vi_neigh" : set(), "ei_region" : set(), "ei_neigh" : set(), "gi_region" : set()} 

    for i in edg:
        e = i-1
        if ver_f[fr[e]-1]==ver_f[to[e]-1]:
            re = ver_f[fr[e]-1]
            regions[re]["ei_region"].add(i)
        else:
            re1 = ver_f[fr[e]-1]
            re2 = ver_f[to[e]-1]
            regions[re1]["ei_region"].add(i)
            regions[re2]["ei_region"].add(i)
            regions[re1]["ei_neigh"].add(i)
            regions[re2]["ei_neigh"].add(i)
            regions[re1]["vi_region_all"].add(to[e])
            regions[re1]["vi_neigh"].update([to[e],fr[e]])
            regions[re2]["vi_region_all"].add(fr[e])
            regions[re2]["vi_neigh"].update([to[e],fr[e]])

    for g in gen:
        regions[ver_f[g[0]-1]]["gi_region"].add(g)
    
    return(regions)






