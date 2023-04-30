# The main method that implements the DiCARP algorithm and returns the results
# and statistics of the algorithm. 
#
# Copyright (c) 2023, by 
# Mehdi Karimi


from opt_agent_2 import *
from OPF_opt import *
from collections import defaultdict
from radial_decom import *
import matplotlib.pyplot as plt
import numpy


problem = "case89"
folder = "Data/"
tol = .0001
rho = 50


def DCA_algorithm(problem, tol = .0001, rho = 50, iter = 3000):
    regions = create_regions(folder+problem)
    # print(regions)
    num_R = len(regions)
    # print(num_R)

    result = {}
    result["num_R"] = num_R
    result["regions"] = regions

    inst = {}
    reg_bel = defaultdict(list)
    reg_bel_eg = defaultdict(list)
    for i in range(num_R):
        data = DataPortal()
        data["vi_region_all"] = regions[i+1]["vi_region_all"]
        data["vi_region"] = regions[i+1]["vi_region"]
        data["vi_neigh"] = regions[i+1]["vi_neigh"]
        data["ei_region"] = regions[i+1]["ei_region"]
        data["ei_neigh"] = regions[i+1]["ei_neigh"]
        data["gi_region"] = regions[i+1]["gi_region"]

        # Paramter rho for the consensus algorithm
        data["rho"] = {None: rho}
        
        data.load(filename= folder+problem+".dat", model=mo_agent_2)
        inst[i+1] = mo_agent_2.create_instance(data, name= "instance--"+str(i+1))
        for nod in list(inst[i+1].vi_neigh):
            reg_bel[nod].append(i+1)
        for eg in list(inst[i+1].ei_neigh):
            reg_bel_eg[eg].append(i+1)


    beta_copy_m = {}
    beta_copy_a = {}
    beta_copy_pt = {}
    beta_copy_pf = {}
    beta_copy_qt = {}
    beta_copy_qf = {}
    res_progress ={i+1:[] for i in range(num_R)}
    min_res_progress = []

    for i in range(iter):
        for m in inst:
            SolverFactory('ipopt').solve(inst[m], tee=False)
        
        for j in reg_bel:
            # temp_vm = sum([inst[k].v_vm[j]+(1/inst[k].dm[j])*inst[k].ym[j] for k in reg_bel[j]]) / len(reg_bel[j])
            # temp_va = sum([inst[k].v_va[j]+(1/inst[k].da[j])*inst[k].ya[j] for k in reg_bel[j]]) / len(reg_bel[j])
            temp_vm = sum([inst[k].dm[j]*inst[k].v_vm[j]+inst[k].ym[j] for k in reg_bel[j]]) / sum([inst[k].dm[j] for k in reg_bel[j]])
            temp_va = sum([inst[k].da[j]*inst[k].v_va[j]+inst[k].ya[j] for k in reg_bel[j]]) / sum([inst[k].da[j] for k in reg_bel[j]])
            for k in reg_bel[j]:
                beta_copy_m[(k,j)] = value(inst[k].beta_vm[j])
                beta_copy_a[(k,j)] = value(inst[k].beta_va[j])
                inst[k].beta_vm[j] = temp_vm
                inst[k].beta_va[j] = temp_va
        for j in reg_bel:
            for k in reg_bel[j]:
                inst[k].ym[j] = inst[k].ym[j]+inst[k].dm[j]*(inst[k].v_vm[j]-inst[k].beta_vm[j])
                inst[k].ya[j] = inst[k].ya[j]+inst[k].da[j]*(inst[k].v_va[j]-inst[k].beta_va[j])
        
        for j in reg_bel_eg:
            temp_pt = sum([inst[k].dpq[j]*inst[k].v_pt[j]+inst[k].ypt[j] for k in reg_bel_eg[j]]) / sum([inst[k].dpq[j] for k in reg_bel_eg[j]])
            temp_pf = sum([inst[k].dpq[j]*inst[k].v_pf[j]+inst[k].ypf[j] for k in reg_bel_eg[j]]) / sum([inst[k].dpq[j] for k in reg_bel_eg[j]])
            temp_qt = sum([inst[k].dpq[j]*inst[k].v_qt[j]+inst[k].yqt[j] for k in reg_bel_eg[j]]) / sum([inst[k].dpq[j] for k in reg_bel_eg[j]])
            temp_qf = sum([inst[k].dpq[j]*inst[k].v_qf[j]+inst[k].yqf[j] for k in reg_bel_eg[j]]) / sum([inst[k].dpq[j] for k in reg_bel_eg[j]])
            for k in reg_bel_eg[j]:
                beta_copy_pt[(k,j)] = value(inst[k].beta_pt[j])
                beta_copy_pf[(k,j)] = value(inst[k].beta_pf[j])
                beta_copy_qt[(k,j)] = value(inst[k].beta_qt[j])
                beta_copy_qf[(k,j)] = value(inst[k].beta_qf[j])
                inst[k].beta_pt[j] = temp_pt
                inst[k].beta_pf[j] = temp_pf
                inst[k].beta_qt[j] = temp_qt
                inst[k].beta_qf[j] = temp_qf
        for j in reg_bel_eg:
            for k in reg_bel_eg[j]:
                inst[k].ypt[j] = inst[k].ypt[j]+inst[k].dpq[j]*(inst[k].v_pt[j]-inst[k].beta_pt[j])
                inst[k].ypf[j] = inst[k].ypf[j]+inst[k].dpq[j]*(inst[k].v_pf[j]-inst[k].beta_pf[j])
                inst[k].yqt[j] = inst[k].yqt[j]+inst[k].dpq[j]*(inst[k].v_qt[j]-inst[k].beta_qt[j])
                inst[k].yqf[j] = inst[k].yqf[j]+inst[k].dpq[j]*(inst[k].v_qf[j]-inst[k].beta_qf[j])
                

        temp_d = {}
        
        parts = ["mp", "md", "ap", "ad", "ptp", "ptd", "pfp", "pfd", "qtp", "qtd", "qfp", "qfd"]
        res = {}
        res_norm = {}
        for p in parts:
            res[p] = []
            res_norm[p] = []

        for k in range(1,num_R+1,1):
            temp_d_m = value(list(inst[k].dm.values())[0])
            temp_d_a = value(list(inst[k].da.values())[0])
            temp_d_pq = value(list(inst[k].dpq.values())[0])


            def  calculate_residuals(inst_v, inst_beta, inst_y, lis, inst_d, i_copy, k, res_p, res_d, res_n_p, res_n_d):
                f2 = sum([(inst_v[i]-inst_beta[i])**2 for i in lis])
                f3 = sum([(inst_d[i]*(i_copy[(k,i)]-inst_beta[i]))**2 for i in lis])
                n1 = sum([(inst_beta[i])**2 for i in lis])
                n2 = sum([(inst_v[i])**2 for i in lis])
                n3 = sum([(inst_y[i])**2 for i in lis])
                res_p.append(value(f2))
                res_d.append(value(f3))
                res_n_p.append((value(n1),value(n2)))
                res_n_d.append(value(n3))
                
        
            
            calculate_residuals(inst[k].v_vm, inst[k].beta_vm, inst[k].ym, inst[k].vi_neigh, inst[k].dm, beta_copy_m, k, res["mp"], res["md"], res_norm["mp"], res_norm["md"])
            calculate_residuals(inst[k].v_va, inst[k].beta_va, inst[k].ya, inst[k].vi_neigh, inst[k].da, beta_copy_a, k, res["ap"], res["ad"], res_norm["ap"], res_norm["ad"])
            calculate_residuals(inst[k].v_pt, inst[k].beta_pt, inst[k].ypt, inst[k].ei_neigh, inst[k].dpq, beta_copy_pt, k, res["ptp"], res["ptd"], res_norm["ptp"], res_norm["ptd"])
            calculate_residuals(inst[k].v_pf, inst[k].beta_pf, inst[k].ypf, inst[k].ei_neigh, inst[k].dpq, beta_copy_pf, k, res["pfp"], res["pfd"], res_norm["pfp"], res_norm["pfd"])
            calculate_residuals(inst[k].v_qt, inst[k].beta_qt, inst[k].yqt, inst[k].ei_neigh, inst[k].dpq, beta_copy_qt, k, res["qtp"], res["qtd"], res_norm["qtp"], res_norm["qtd"])
            calculate_residuals(inst[k].v_qf, inst[k].beta_qf, inst[k].yqf, inst[k].ei_neigh, inst[k].dpq, beta_copy_qf, k, res["qfp"], res["qfd"], res_norm["qfp"], res_norm["qfd"])

            
            temp_d[k] = (temp_d_m, temp_d_a, temp_d_pq)
            

        for j in reg_bel:
            for k in reg_bel[j]:
                inst[k].dm[j] = temp_d[k][0]
                inst[k].da[j] = temp_d[k][1]
        for j in reg_bel_eg:
            for k in reg_bel_eg[j]:
                inst[k].dpq[j] = temp_d[k][2]
        


        P_da = sum([value(inst[k].obj) for k in range(1,num_R+1,1)])
        print("Problem {} : iteration is {}, current obj value is {} ".format(problem,  i, P_da))

        res_com = 0
        
        parts_p = ["mp", "ap", "ptp", "pfp", "qtp", "qfp"]
        parts_d = ["md", "ad", "ptd", "pfd", "qtd", "qfd"]
        for k in range(1,num_R+1,1):
            temp_r_p = sum([res[p][k-1] for p in parts_p])
            temp_r_d = sum([res[p][k-1] for p in parts_d])
            n1 = sum([res_norm[p][k-1][0] for p in parts_p])
            n2 = sum([res_norm[p][k-1][1] for p in parts_p])
            n3 = sum([res_norm[p][k-1] for p in parts_d])
            temp_res = max(sqrt(temp_r_p)/max([sqrt(n1), sqrt(n2)]), sqrt(temp_r_d)/sqrt(n3))
            res_com = max([res_com, temp_res])
            res_progress[k].append(temp_res)
        min_res_progress.append(res_com)

        if res_com < tol:
            break
    
    result["itr"] = i
    result["res_prog"] = res_progress
    result["min_res_prog"] = min_res_progress

    inst_ipm = mo.create_instance(folder+problem+".dat")

    SolverFactory('ipopt').solve(inst_ipm)
    # print(ipopt_res.solver.status)
    P_ipm = value(inst_ipm.obj)
    print("Optimal value by IPM = ", P_ipm)

    GAP = abs(P_ipm - P_da) / P_ipm
    
    
    result["P_ipm"] = P_ipm
    result["P_da"] = P_da
    result["GAP"] = GAP 
    
    plt.clf()
    plt.plot(numpy.log10(min_res_progress[:]))

    plt.xlabel('# of iterations')
    plt.ylabel('Log of the minimum residual')
    # plt.title("A simple line graph")
    plt.savefig('Figures/fig-'+problem+'-'+str(rho)+'.png')
    # plt.show()

    return result

