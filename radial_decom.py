# The main method applies Radial Partitioning to a given connected graph. 
# The graph must be given as a Pyomo DAT file.
#
# Copyright (c) 2024, by 
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
    
    reg_bel = defaultdict(list)

    while len(ver_c):
        root = ver_c[0] # 
        # root = random.choice(ver_c) #
        ver_f[root-1] = num_reg
        se = set([root])
        # reg_bel[root].append(num_reg)
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
                # reg_bel[cur].append(num_reg)
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
            # reg_bel[to[e]].append(re1)
            regions[re1]["vi_neigh"].update([to[e],fr[e]])
            regions[re2]["vi_region_all"].add(fr[e])
            # reg_bel[fr[e]].append(re2)
            regions[re2]["vi_neigh"].update([to[e],fr[e]])
    
    for i in range(num_reg):
            
        for nod in list(regions[i+1]["vi_neigh"]):
            reg_bel[nod].append(i+1)
            
    for g in gen:
        # for x in reg_bel[g[0]]:
        #     regions[x]["gi_region"].add(g[0])
        # if g[0]-1 in reg_bel:
        #     for i in reg_bel[g[0]-1]:
        #         regions[i]["gi_region"].add(g)
        # else:
            regions[ver_f[g[0]-1]]["gi_region"].add(g)
    return(regions)


def create_regions_2(problem):
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
        root = ver_c[0] # 
        # root = random.choice(ver_c) #
        ver_f[root-1] = num_reg
        se = set([root])
        # reg_bel[root].append(num_reg)
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
                # reg_bel[cur].append(num_reg)
                stack = stack + [(i,cur) for i in temp if ver_f[i-1]==0]
        reg_n[num_reg] = se
        num_reg = num_reg+1

    num_reg = num_reg-1
    regions = {}
    for i in range(num_reg):
        regions[i+1] = {"vi_region_all" : reg_n[i+1].copy(), "vi_region" : reg_n[i+1].copy(),
        "vi_neigh" : set(), "ei_region" : set(), "ei_neigh" : set(), "gi_region" : set(), "gi_neigh" : set()} 

    for i in edg:
        e = i-1
        if ver_f[fr[e]-1]!=ver_f[to[e]-1]:
            re1 = ver_f[fr[e]-1]
            re2 = ver_f[to[e]-1]
            regions[re1]["vi_region_all"].add(to[e])
            # reg_bel[to[e]].append(re1)
            regions[re1]["vi_neigh"].update([to[e],fr[e]])
            regions[re2]["vi_region_all"].add(fr[e])
            # reg_bel[fr[e]].append(re2)
            regions[re2]["vi_neigh"].update([to[e],fr[e]])
    
    reg_bel_all = defaultdict(set)
    reg_bel = defaultdict(set)
    for i in range(num_reg):
        for nod in list(regions[i+1]["vi_region_all"]):
            reg_bel_all[nod].add(i+1)
        for nod in list(regions[i+1]["vi_neigh"]):
            reg_bel[nod].add(i+1)

    for i in edg:
        e = i-1
        n1 = fr[e]
        n2 = to[e]
        set1 = reg_bel_all[n1]
        set2 = reg_bel_all[n2]
        int_set = set1.intersection(set2)
        if len(int_set) == 1:
            for val in int_set:
                regions[val]["ei_region"].add(i)
        else:
            for val in int_set:
                regions[val]["ei_region"].add(i)
                regions[val]["ei_neigh"].add(i)
            


    for g in gen:
        # for x in reg_bel[g[0]]:
        #     regions[x]["gi_region"].add(g[0])
        # if g[0]-1 in reg_bel:
        #     for i in reg_bel[g[0]-1]:
        #         regions[i]["gi_region"].add(g)
        # else:
        regions[ver_f[g[0]-1]]["gi_region"].add(g)
        if g[0] in reg_bel:
            for k in reg_bel[g[0]]:
                # print(k)
                regions[k]["gi_neigh"].add(g)
        
    return(regions)



