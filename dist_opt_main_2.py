# The main method that implements the DiCARP algorithm and returns the results
# and statistics of the algorithm. 
#
# Copyright (c) 2024, by 
# Mehdi Karimi


from opt_agent_2 import *
from OPF_opt import *
from collections import defaultdict
from radial_decom import *
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
import pickle
from statistics import mean


folder = "Data/"
tol = .0001
rho = 50


def DCA_algorithm(problem, tol = .0001, rhov = 1000, rhop = 100,  iter = 3000, adaptive = False, component = False):
    if component:
        regions = create_regions_2_component(folder+problem)
    else:
        regions = create_regions_2(folder+problem)
    # print(regions)
    num_R = len(regions)
    # print(num_R)

    result = {}
    result["num_R"] = num_R
    result["regions"] = regions
    
    file_path = folder+problem+'.pkl'
    if False and os.path.exists(file_path):
        with open(file_path, 'rb') as file:
            inst, reg_bel, reg_bel_eg, reg_bel_g = pickle.load(file)
            print("File loaded successfully.")
    else:
        inst = {}
        reg_bel = defaultdict(list)
        reg_bel_eg = defaultdict(list)
        reg_bel_g = defaultdict(list)
        for i in range(num_R):
            data = DataPortal()
            data["vi_region_all"] = regions[i+1]["vi_region_all"]
            data["vi_region"] = regions[i+1]["vi_region"]
            data["vi_neigh"] = regions[i+1]["vi_neigh"]
            data["ei_region"] = regions[i+1]["ei_region"]
            data["ei_neigh"] = regions[i+1]["ei_neigh"]
            data["gi_region"] = regions[i+1]["gi_region"]
            data["gi_neigh"] = regions[i+1]["gi_neigh"]

            # Paramter rho for the consensus algorithm
            data["rhov"] = {None: rhov}
            data["rhop"] = {None: rhop}
            
            data.load(filename= folder+problem+".dat", model=mo_agent_2)
            inst[i+1] = mo_agent_2.create_instance(data, name= "instance--"+str(i+1))
            for nod in list(inst[i+1].vi_neigh):
                reg_bel[nod].append(i+1)
            for eg in list(inst[i+1].ei_neigh):
                reg_bel_eg[eg].append(i+1)
            for g in list(inst[i+1].gi_neigh):
                reg_bel_g[g].append(i+1)

        with open(file_path, 'wb') as file:
            pickle.dump((inst, reg_bel, reg_bel_eg, reg_bel_g), file)
            

    # type_shared_var = ["m", "a", "pt", "pf", "qt", "qf"]
    # beta_copy={i: {} for i in type_shared_var}
    # y_copy={i: {} for i in type_shared_var}
    
    res_progress ={i+1:[] for i in range(num_R)}
    min_res_progress = []
    
    print(reg_bel)
    
    result_prog = {}
    adap_w = {"dm": {}, "da": {}, "dp":  {}, "dq" : {}, "dpg" : {}, "dqg" : {}}
    for i in range(iter):
        inst_copy = copy.deepcopy(inst)
        # if i==1:
        #     inst_copy_5 = copy.deepcopy(inst)
        solv_time = []
        for m in inst:
            ipopt_res = SolverFactory('ipopt').solve(inst[m], tee=False)
            solv_time.append(ipopt_res.Solver.Time)
        # print("avg:", np.average(np.array(solv_time)))
        # print("max:", np.max(np.array(solv_time)))
        

        for j in reg_bel:
            # temp_vm = sum([inst[k].v_vm[j]+(1/inst[k].dm[j])*inst[k].ym[j] for k in reg_bel[j]]) / len(reg_bel[j])
            # temp_va = sum([inst[k].v_va[j]+(1/inst[k].da[j])*inst[k].ya[j] for k in reg_bel[j]]) / len(reg_bel[j])
            temp_vm = sum([inst[k].dm[j]*inst[k].v_vm[j]+inst[k].ym[j] for k in reg_bel[j]]) / sum([inst[k].dm[j] for k in reg_bel[j]])
            temp_va = sum([inst[k].da[j]*inst[k].v_va[j]+inst[k].ya[j] for k in reg_bel[j]]) / sum([inst[k].da[j] for k in reg_bel[j]])
            for k in reg_bel[j]:
                inst[k].beta_vm[j] = temp_vm
                inst[k].beta_va[j] = temp_va
        for j in reg_bel:
            for k in reg_bel[j]:
                inst[k].yhm[j]  = inst[k].ym[j]+inst[k].dm[j]*(inst[k].v_vm[j]-inst_copy[k].beta_vm[j])
                inst[k].yha[j] = inst[k].ya[j]+inst[k].da[j]*(inst[k].v_va[j]-inst_copy[k].beta_va[j])
                inst[k].ym[j] = inst[k].ym[j]+inst[k].dm[j]*(inst[k].v_vm[j]-inst[k].beta_vm[j])
                inst[k].ya[j] = inst[k].ya[j]+inst[k].da[j]*(inst[k].v_va[j]-inst[k].beta_va[j])
        
        for j in reg_bel_eg:
            temp_pt = sum([inst[k].dp[j]*inst[k].v_pt[j]+inst[k].ypt[j] for k in reg_bel_eg[j]]) / sum([inst[k].dp[j] for k in reg_bel_eg[j]])
            temp_pf = sum([inst[k].dp[j]*inst[k].v_pf[j]+inst[k].ypf[j] for k in reg_bel_eg[j]]) / sum([inst[k].dp[j] for k in reg_bel_eg[j]])
            temp_qt = sum([inst[k].dq[j]*inst[k].v_qt[j]+inst[k].yqt[j] for k in reg_bel_eg[j]]) / sum([inst[k].dq[j] for k in reg_bel_eg[j]])
            temp_qf = sum([inst[k].dq[j]*inst[k].v_qf[j]+inst[k].yqf[j] for k in reg_bel_eg[j]]) / sum([inst[k].dq[j] for k in reg_bel_eg[j]])
            for k in reg_bel_eg[j]:
                inst[k].beta_pt[j] = temp_pt
                inst[k].beta_pf[j] = temp_pf
                inst[k].beta_qt[j] = temp_qt
                inst[k].beta_qf[j] = temp_qf
        for j in reg_bel_eg:
            for k in reg_bel_eg[j]:
                inst[k].yhpt[j] = inst[k].ypt[j]+inst[k].dp[j]*(inst[k].v_pt[j]-inst_copy[k].beta_pt[j])
                inst[k].yhpf[j] = inst[k].ypf[j]+inst[k].dp[j]*(inst[k].v_pf[j]-inst_copy[k].beta_pf[j])
                inst[k].yhqt[j] = inst[k].yqt[j]+inst[k].dq[j]*(inst[k].v_qt[j]-inst_copy[k].beta_qt[j])
                inst[k].yhqf[j] = inst[k].yqf[j]+inst[k].dq[j]*(inst[k].v_qf[j]-inst_copy[k].beta_qf[j])
                inst[k].ypt[j] = inst[k].ypt[j]+inst[k].dp[j]*(inst[k].v_pt[j]-inst[k].beta_pt[j])
                inst[k].ypf[j] = inst[k].ypf[j]+inst[k].dp[j]*(inst[k].v_pf[j]-inst[k].beta_pf[j])
                inst[k].yqt[j] = inst[k].yqt[j]+inst[k].dq[j]*(inst[k].v_qt[j]-inst[k].beta_qt[j])
                inst[k].yqf[j] = inst[k].yqf[j]+inst[k].dq[j]*(inst[k].v_qf[j]-inst[k].beta_qf[j])
                
        for j in reg_bel_g:
            temp_pg = sum([inst[k].dpg[j]*inst[k].v_pg[j]+inst[k].ypg[j] for k in reg_bel_g[j]]) / sum([inst[k].dpg[j] for k in reg_bel_g[j]])
            temp_qg = sum([inst[k].dqg[j]*inst[k].v_qg[j]+inst[k].yqg[j] for k in reg_bel_g[j]]) / sum([inst[k].dqg[j] for k in reg_bel_g[j]])
            for k in reg_bel_g[j]:
                inst[k].beta_pg[j] = temp_pg
                inst[k].beta_qg[j] = temp_qg
        for j in reg_bel_g:
            for k in reg_bel_g[j]:
                inst[k].yhpg[j] = inst[k].ypg[j]+inst[k].dpg[j]*(inst[k].v_pg[j]-inst_copy[k].beta_pg[j])
                inst[k].yhqg[j] = inst[k].yqg[j]+inst[k].dqg[j]*(inst[k].v_qg[j]-inst_copy[k].beta_qg[j])
                inst[k].ypg[j] = inst[k].ypg[j]+inst[k].dpg[j]*(inst[k].v_pg[j]-inst[k].beta_pg[j])
                inst[k].yqg[j] = inst[k].yqg[j]+inst[k].dqg[j]*(inst[k].v_qg[j]-inst[k].beta_qg[j])
                
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
            temp_d_p = value(list(inst[k].dp.values())[0])
            temp_d_q = value(list(inst[k].dq.values())[0])
                
            def  calculate_residuals(inst_v, inst_beta, inst_y, lis, inst_d, i_copy, k, res_p, res_d, res_n_p, res_n_d):
                f2 = sum([(inst_v[i]-inst_beta[i])**2 for i in lis])
                f3 = sum([(inst_d[i]*(i_copy[i]-inst_beta[i]))**2 for i in lis])
                n1 = sum([(inst_beta[i])**2 for i in lis])
                n2 = sum([(inst_v[i])**2 for i in lis])
                n3 = sum([(inst_y[i])**2 for i in lis])
                res_p.append(value(f2))
                res_d.append(value(f3))
                res_n_p.append((value(n1),value(n2)))
                res_n_d.append(value(n3))
        
            
            calculate_residuals(inst[k].v_vm, inst[k].beta_vm, inst[k].ym, inst[k].vi_neigh, inst[k].dm, inst_copy[k].beta_vm, k, res["mp"], res["md"], res_norm["mp"], res_norm["md"])
            calculate_residuals(inst[k].v_va, inst[k].beta_va, inst[k].ya, inst[k].vi_neigh, inst[k].da, inst_copy[k].beta_va, k, res["ap"], res["ad"], res_norm["ap"], res_norm["ad"])
            calculate_residuals(inst[k].v_pt, inst[k].beta_pt, inst[k].ypt, inst[k].ei_neigh, inst[k].dp, inst_copy[k].beta_pt, k, res["ptp"], res["ptd"], res_norm["ptp"], res_norm["ptd"])
            calculate_residuals(inst[k].v_pf, inst[k].beta_pf, inst[k].ypf, inst[k].ei_neigh, inst[k].dp, inst_copy[k].beta_pf, k, res["pfp"], res["pfd"], res_norm["pfp"], res_norm["pfd"])
            calculate_residuals(inst[k].v_qt, inst[k].beta_qt, inst[k].yqt, inst[k].ei_neigh, inst[k].dq, inst_copy[k].beta_qt, k, res["qtp"], res["qtd"], res_norm["qtp"], res_norm["qtd"])
            calculate_residuals(inst[k].v_qf, inst[k].beta_qf, inst[k].yqf, inst[k].ei_neigh, inst[k].dq, inst_copy[k].beta_qf, k, res["qfp"], res["qfd"], res_norm["qfp"], res_norm["qfd"])

            
            temp_d[k] = (temp_d_m, temp_d_a, temp_d_p, temp_d_q)
        
        if i > 1:
            update_parameter_adaptive_2(inst,inst_copy,num_R,reg_bel,reg_bel_eg,reg_bel_g,adap_w,i%10)
        if adaptive: 
            if i > 10 and i%10 == 9:
                # update_parameter_adaptive(inst,inst_copy,num_R)
                for j in reg_bel:
                    for k in reg_bel[j]:
                        inst[k].dm[j] = mean([adap_w["dm"][(x,k,j)] for x in range(9)])
                        inst[k].da[j] = mean([adap_w["da"][(x,k,j)] for x in range(9)])
                for j in reg_bel_eg:
                    for k in reg_bel_eg[j]:
                        inst[k].dp[j] = mean([adap_w["dp"][(x,k,j)] for x in range(9)]) 
                        inst[k].dq[j] = mean([adap_w["dq"][(x,k,j)] for x in range(9)]) 
                for j in reg_bel_g:
                    for k in reg_bel_g[j]:
                        inst[k].dpg[j] = mean([adap_w["dpg"][(x,k,j)] for x in range(9)]) 
                        inst[k].dqg[j] = mean([adap_w["dqg"][(x,k,j)] for x in range(9)])
                # print(value(inst[1].dm[8]))
        # else:
            # update_parameter(inst, temp_d, reg_bel, reg_bel_eg)
        


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

        if i > 1 and i%100==0:
            result_prog[i] = min_res_progress
            with open(folder+problem+'iter'+str(i)+'_3.pkl', 'wb') as file:
            
                pickle.dump(result_prog, file)
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
    plt.plot(np.log10(min_res_progress[:]))

    plt.xlabel('# of iterations')
    plt.ylabel('Log of the minimum residual')
    # plt.title("A simple line graph")
    path_name = 'Figures/fig-'+problem+'-'+str(rho)+'.png'
    if adaptive:
        path_name = 'Figures/fig-'+problem+'-'+str(rho)+'adapt.png'
    plt.savefig(path_name)
    # plt.show()

    return result


def update_parameter(inst, temp_d, reg_bel, reg_bel_eg):
    for j in reg_bel:
        for k in reg_bel[j]:
            inst[k].dm[j] = temp_d[k][0]
            inst[k].da[j] = temp_d[k][1]
    for j in reg_bel_eg:
        for k in reg_bel_eg[j]:
            inst[k].dp[j] = temp_d[k][2]
            inst[k].dq[j] = temp_d[k][3]


def update_parameter_adaptive(inst,inst_copy,num_R):
    
    def ratio(x,y):
        x = [value(i) for i in x]
        y = [value(i) for i in y]
        x = np.array(x)
        y = np.array(y)
        # temp1 = sum([x[i]**2 for i in range(len(x))])
        # temp2 = sum([x[i]*y[i] for  i in range(len(x)) ])
        nx = np.linalg.norm(x)
        ny = np.linalg.norm(y)
        inner = np.inner(x,y)
        temp_SD = value((nx**2)/inner)
        temp_MG = value(inner/(ny**2))
        corr = inner/(nx*ny)
        if 2*temp_MG > temp_SD:
            return temp_MG, corr
        else:
            return temp_SD-temp_MG/2, corr
    
    def safe_guard(al,be,cor_al,cor_be,d,eps=0.5):
        temp = d
        if cor_al > eps and cor_be > eps:
            temp=sqrt(al*be)
        elif cor_al > eps and cor_be <= eps:
            temp=al
        elif cor_al <= eps and cor_be > eps:
            temp=be
        cost = 2
        if temp > cost*d:
            temp= cost*d
        elif temp < d/cost:
            temp= d/cost
        
        # print(temp)
        return max(min(temp,5000),10)
        
    for k in range(1,num_R+1,1):
        x1 = [(inst[k].yhm[j]-inst_copy[k].yhm[j]) for j in list(inst[k].vi_neigh)]
        y1 = [(inst[k].v_vm[j]-inst_copy[k].v_vm[j]) for j in list(inst[k].vi_neigh)]
        al_m, cor_al_m = ratio(x1,y1)
        x2 = [(inst[k].ym[j]-inst_copy[k].ym[j]) for j in list(inst[k].vi_neigh)]
        y2 = [-(inst[k].beta_vm[j]-inst_copy[k].beta_vm[j]) for j in list(inst[k].vi_neigh)]
        be_m, cor_be_m = ratio(x2,y2)

        x1 = [(inst[k].yha[j]-inst_copy[k].yha[j]) for j in list(inst[k].vi_neigh)]
        y1 = [(inst[k].v_va[j]-inst_copy[k].v_va[j]) for j in list(inst[k].vi_neigh)]
        al_a, cor_al_a = ratio(x1,y1)
        x2 = [(inst[k].ya[j]-inst_copy[k].ya[j]) for j in list(inst[k].vi_neigh)]
        y2 = [-(inst[k].beta_va[j]-inst_copy[k].beta_va[j]) for j in list(inst[k].vi_neigh)]
        be_a, cor_be_a = ratio(x2,y2)
        
        for j in list(inst[k].vi_neigh):
            inst[k].dm[j] = safe_guard(al_m,be_m,cor_al_m,cor_be_m,value(inst[k].dm[j]))
            inst[k].da[j] = safe_guard(al_a,be_a,cor_al_a,cor_be_a,value(inst[k].da[j]))

        x1 = [(inst[k].yhpt[j]-inst_copy[k].yhpt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yhpf[j]-inst_copy[k].yhpf[j]) for j in list(inst[k].ei_neigh)]
        y1 = [(inst[k].v_pt[j]-inst_copy[k].v_pt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].v_pf[j]-inst_copy[k].v_pf[j]) for j in list(inst[k].ei_neigh)]
        al_p, cor_al_p = ratio(x1,y1)
        x2 = [(inst[k].ypt[j]-inst_copy[k].ypt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].ypf[j]-inst_copy[k].ypf[j]) for j in list(inst[k].ei_neigh)]
        y2 = [-(inst[k].beta_pt[j]-inst_copy[k].beta_pt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].beta_pf[j]-inst_copy[k].beta_pf[j]) for j in list(inst[k].ei_neigh)]
        be_p, cor_be_p= ratio(x2,y2)
        
        x1 = [(inst[k].yhqt[j]-inst_copy[k].yhqt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yhqf[j]-inst_copy[k].yhqf[j]) for j in list(inst[k].ei_neigh)]
        y1 = [(inst[k].v_qt[j]-inst_copy[k].v_qt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].v_qf[j]-inst_copy[k].v_qf[j]) for j in list(inst[k].ei_neigh)]
        al_q, cor_al_q = ratio(x1,y1)
        x2 = [(inst[k].yqt[j]-inst_copy[k].yqt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yqf[j]-inst_copy[k].yqf[j]) for j in list(inst[k].ei_neigh)]
        y2 = [-(inst[k].beta_qt[j]-inst_copy[k].beta_qt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].beta_qf[j]-inst_copy[k].beta_qf[j]) for j in list(inst[k].ei_neigh)]
        be_q, cor_be_q= ratio(x2,y2)

        for j in list(inst[k].ei_neigh):
            inst[k].dp[j] = safe_guard(al_p,be_p,cor_al_p,cor_be_p,value(inst[k].dp[j]))
            inst[k].dq[j] = safe_guard(al_q,be_q,cor_al_q,cor_be_q,value(inst[k].dq[j]))


def update_parameter_adaptive_2(inst,inst_copy,num_R,reg_bel,reg_bel_eg,reg_bel_g,adap_w,iter):
    
    def ratio(x,y):
        x = [value(i) for i in x]
        y = [value(i) for i in y]
        x = np.array(x)
        y = np.array(y)
        # temp1 = sum([x[i]**2 for i in range(len(x))])
        # temp2 = sum([x[i]*y[i] for  i in range(len(x)) ])
        nx = np.linalg.norm(x)
        ny = np.linalg.norm(y)
        inner = np.inner(x,y)
        temp_SD = value((nx**2)/inner) if inner != 0 else 0
        temp_MG = value(inner/(ny**2)) if ny != 0 else 0
        corr = inner/(nx*ny) if ny !=0 else 0
        if 2*temp_MG > temp_SD:
            return temp_MG, corr
        else:
            return temp_SD-temp_MG/2, corr
    
    def safe_guard(al,be,cor_al,cor_be,d,eps=0.2):
        temp = d
        if cor_al > eps and cor_be > eps:
            temp=sqrt(al*be)
        elif cor_al > eps and cor_be <= eps:
            temp=al
        elif cor_al <= eps and cor_be > eps:
            temp=be
        cost = 1.2
        if temp > cost*d:
            temp= cost*d
        elif temp < d/cost:
            temp= d/cost
        
        # print(temp)
        return max(min(temp,20000),10)

    for j in reg_bel:
        x1 = [(inst[k].yhm[j]-inst_copy[k].yhm[j]) for k in reg_bel[j]]
        y1 = [(inst[k].v_vm[j]-inst_copy[k].v_vm[j]) for k in reg_bel[j]]
        al_m, cor_al_m = ratio(x1,y1)
        x2 = [(inst[k].ym[j]-inst_copy[k].ym[j]) for k in reg_bel[j]]
        y2 = [-(inst[k].beta_vm[j]-inst_copy[k].beta_vm[j]) for k in reg_bel[j]]
        be_m, cor_be_m = ratio(x2,y2)

        x1 = [(inst[k].yha[j]-inst_copy[k].yha[j]) for k in reg_bel[j]]
        y1 = [(inst[k].v_va[j]-inst_copy[k].v_va[j]) for k in reg_bel[j]]
        al_a, cor_al_a = ratio(x1,y1)
        x2 = [(inst[k].ya[j]-inst_copy[k].ya[j]) for k in reg_bel[j]]
        y2 = [-(inst[k].beta_va[j]-inst_copy[k].beta_va[j]) for k in reg_bel[j]]
        be_a, cor_be_a = ratio(x2,y2)
        
        for k in reg_bel[j]:
            adap_w["dm"][(iter,k,j)] = safe_guard(al_m,be_m,cor_al_m,cor_be_m,value(inst[k].dm[j]))
            adap_w["da"][(iter,k,j)] = safe_guard(al_a,be_a,cor_al_a,cor_be_a,value(inst[k].da[j]))
            # inst[k].dm[j] = safe_guard(al_m,be_m,cor_al_m,cor_be_m,value(inst[k].dm[j]))
            # inst[k].da[j] = safe_guard(al_a,be_a,cor_al_a,cor_be_a,value(inst[k].da[j]))

    for j in reg_bel_eg:
        x1 = [(inst[k].yhpt[j]-inst_copy[k].yhpt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yhpf[j]-inst_copy[k].yhpf[j]) for k in reg_bel_eg[j]]
        y1 = [(inst[k].v_pt[j]-inst_copy[k].v_pt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].v_pf[j]-inst_copy[k].v_pf[j]) for k in reg_bel_eg[j]]
        al_p, cor_al_p = ratio(x1,y1)
        x2 = [(inst[k].ypt[j]-inst_copy[k].ypt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].ypf[j]-inst_copy[k].ypf[j]) for k in reg_bel_eg[j]]
        y2 = [-(inst[k].beta_pt[j]-inst_copy[k].beta_pt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].beta_pf[j]-inst_copy[k].beta_pf[j]) for k in reg_bel_eg[j]]
        be_p, cor_be_p= ratio(x2,y2)
        
        x1 = [(inst[k].yhqt[j]-inst_copy[k].yhqt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yhqf[j]-inst_copy[k].yhqf[j]) for k in reg_bel_eg[j]]
        y1 = [(inst[k].v_qt[j]-inst_copy[k].v_qt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].v_qf[j]-inst_copy[k].v_qf[j]) for k in reg_bel_eg[j]]
        al_q, cor_al_q = ratio(x1,y1)
        x2 = [(inst[k].yqt[j]-inst_copy[k].yqt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].yqf[j]-inst_copy[k].yqf[j]) for k in reg_bel_eg[j]]
        y2 = [-(inst[k].beta_qt[j]-inst_copy[k].beta_qt[j]) for j in list(inst[k].ei_neigh)]+[(inst[k].beta_qf[j]-inst_copy[k].beta_qf[j]) for k in reg_bel_eg[j]]
        be_q, cor_be_q= ratio(x2,y2)

        for k in reg_bel_eg[j]:
            adap_w["dp"][(iter,k,j)] = safe_guard(al_p,be_p,cor_al_p,cor_be_p,value(inst[k].dp[j]))
            adap_w["dq"][(iter,k,j)] = safe_guard(al_q,be_q,cor_al_q,cor_be_q,value(inst[k].dq[j]))
            # inst[k].dp[j] = safe_guard(al_p,be_p,cor_al_p,cor_be_p,value(inst[k].dp[j]))
            # inst[k].dq[j] = safe_guard(al_q,be_q,cor_al_q,cor_be_q,value(inst[k].dq[j]))
    
    for j in reg_bel_g:
        x1 = [(inst[k].yhpg[j]-inst_copy[k].yhpg[j]) for k in reg_bel_g[j]]
        y1 = [(inst[k].v_pg[j]-inst_copy[k].v_pg[j]) for k in reg_bel_g[j]]
        al_pg, cor_al_pg = ratio(x1,y1)
        x2 = [(inst[k].ypg[j]-inst_copy[k].ypg[j]) for k in reg_bel_g[j]]
        y2 = [-(inst[k].beta_pg[j]-inst_copy[k].beta_pg[j]) for k in reg_bel_g[j]]
        be_pg, cor_be_pg = ratio(x2,y2)

        x1 = [(inst[k].yhqg[j]-inst_copy[k].yhqg[j]) for k in reg_bel_g[j]]
        y1 = [(inst[k].v_qg[j]-inst_copy[k].v_qg[j]) for k in reg_bel_g[j]]
        al_qg, cor_al_qg = ratio(x1,y1)
        x2 = [(inst[k].yqg[j]-inst_copy[k].yqg[j]) for k in reg_bel_g[j]]
        y2 = [-(inst[k].beta_qg[j]-inst_copy[k].beta_qg[j]) for k in reg_bel_g[j]]
        be_qg, cor_be_qg = ratio(x2,y2)
        
        for k in reg_bel_g[j]:
            adap_w["dpg"][(iter,k,j)] = safe_guard(al_pg,be_pg,cor_al_pg,cor_be_pg,value(inst[k].dpg[j]))
            adap_w["dqg"][(iter,k,j)] = safe_guard(al_qg,be_qg,cor_al_qg,cor_be_qg,value(inst[k].dqg[j]))
