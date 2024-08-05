# This file creates a Pyomo Abstract model for each of sub-problems aganets for the OPF problem. 
#
# Copyright (c) 2023, by 
# Mehdi Karimi


from pyomo.environ import *
mo_agent_2 = AbstractModel()

mo_agent_2.baseMVA = Param()
mo_agent_2.nv = Param(within=NonNegativeIntegers)
mo_agent_2.ne = Param(within=NonNegativeIntegers)
mo_agent_2.ng = Param(within=NonNegativeIntegers)


mo_agent_2.vi = RangeSet(1,mo_agent_2.nv)
mo_agent_2.ei = RangeSet(1,mo_agent_2.ne)

mo_agent_2.gi = Set(within=mo_agent_2.vi * NonNegativeIntegers)

mo_agent_2.vi_region_all = Set(ordered=False, within=mo_agent_2.vi)
mo_agent_2.vi_region = Set(ordered=False, within=mo_agent_2.vi_region_all)
mo_agent_2.vi_neigh = Set(ordered=False, within=mo_agent_2.vi_region_all)
mo_agent_2.ei_region = Set(ordered=False, within=mo_agent_2.ei)
mo_agent_2.ei_neigh = Set(ordered=False, within=mo_agent_2.ei)
mo_agent_2.gi_region = Set(ordered=False, within=mo_agent_2.gi)
mo_agent_2.gi_neigh = Set(ordered=False, within=mo_agent_2.gi)

mo_agent_2.slack = Param(within=mo_agent_2.vi)
mo_agent_2.angslack = Param()

mo_agent_2.pL = Param(mo_agent_2.vi, within=Reals)
mo_agent_2.qL = Param(mo_agent_2.vi, within=Reals)

mo_agent_2.vL = Param(mo_agent_2.vi, within=Reals)
mo_agent_2.vU = Param(mo_agent_2.vi, within=Reals)

mo_agent_2.gs = Param(mo_agent_2.vi, within=Reals)
mo_agent_2.bs = Param(mo_agent_2.vi, within=Reals)

mo_agent_2.fr = Param(mo_agent_2.ei, within = mo_agent_2.vi)
mo_agent_2.to = Param(mo_agent_2.ei, within = mo_agent_2.vi)

mo_agent_2.Y_ind = Set(initialize=['Y11r', 'Y11i', 'Y12r', 'Y12i', 'Y21r', 'Y21i', 'Y22r', 'Y22i'])
mo_agent_2.Y = Param(mo_agent_2.ei, mo_agent_2.Y_ind, within = Reals)
mo_agent_2.e_con = Param(mo_agent_2.ei, within = Reals)
mo_agent_2.e_sus = Param(mo_agent_2.ei, within = Reals)
mo_agent_2.e_bsh = Param(mo_agent_2.ei, within = Reals)

mo_agent_2.e_su = Param(mo_agent_2.ei, within= Reals)
mo_agent_2.e_dl = Param(mo_agent_2.ei, within= Reals)
mo_agent_2.e_du = Param(mo_agent_2.ei, within= Reals)


mo_agent_2.pGl = Param(mo_agent_2.gi)
mo_agent_2.pGu = Param(mo_agent_2.gi)
mo_agent_2.qGl = Param(mo_agent_2.gi)
mo_agent_2.qGu = Param(mo_agent_2.gi)

mo_agent_2.c2 = Param(mo_agent_2.gi, within = NonNegativeReals)
mo_agent_2.c1 = Param(mo_agent_2.gi, within = NonNegativeReals)
mo_agent_2.c0 = Param(mo_agent_2.gi, within = NonNegativeReals)

def init_rule_pg(m,i,j):
    i = (i,j)
    if m.pGl[i] == '-Inf' and m.pGu[i] == 'Inf':
      return 0
    elif m.pGl[i] == '-Inf':
      return 0
    elif m.pGu[i] == 'Inf':
      return 0
    else: 
      return  (m.pGl[i]+m.pGu[i])/2
mo_agent_2.v_pg = Var(mo_agent_2.gi_region | mo_agent_2.gi_neigh, initialize=init_rule_pg)
def init_rule_qg(m,i,j):
    i = (i,j)
    if m.qGl[i] == '-Inf' and m.qGu[i] == 'Inf':
      return 0
    elif m.qGl[i] == '-Inf':
      return 0
    elif m.qGu[i] == 'Inf':
      return 0
    else: 
      return  (m.qGl[i]+m.qGu[i])/2
mo_agent_2.v_qg = Var(mo_agent_2.gi_region | mo_agent_2.gi_neigh, initialize=init_rule_qg)
def init_rule_vm(model, i):
    return 1
mo_agent_2.v_vm = Var(mo_agent_2.vi_region_all, initialize=init_rule_vm)
def init_rule_va(model, i):
    return 0
mo_agent_2.v_va = Var(mo_agent_2.vi_region_all, initialize=init_rule_va)
mo_agent_2.v_pt = Var(mo_agent_2.ei_region)
mo_agent_2.v_pf = Var(mo_agent_2.ei_region)
mo_agent_2.v_qt = Var(mo_agent_2.ei_region)
mo_agent_2.v_qf = Var(mo_agent_2.ei_region)

mo_agent_2.rhov = Param(within = NonNegativeReals)
mo_agent_2.rhop = Param(within = NonNegativeReals)

def beta_vm_ini(m,i):
    return 1
mo_agent_2.beta_vm = Param(mo_agent_2.vi_neigh, mutable =True, initialize=beta_vm_ini)

def beta_va_ini(m,i):
    return 0
mo_agent_2.beta_va = Param(mo_agent_2.vi_neigh, mutable =True, initialize=beta_va_ini)

def beta_pt_ini(m,i):
    return 1
mo_agent_2.beta_pt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=beta_pt_ini)

def beta_pf_ini(m,i):
    return 1
mo_agent_2.beta_pf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=beta_pf_ini)

def beta_qt_ini(m,i):
    return 1
mo_agent_2.beta_qt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=beta_qt_ini)

def beta_qf_ini(m,i):
    return 1
mo_agent_2.beta_qf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=beta_qf_ini)

def beta_pg_ini(m,i,j):
    return 1
mo_agent_2.beta_pg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=beta_pg_ini)

def beta_qg_ini(m,i,j):
    return 1
mo_agent_2.beta_qg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=beta_qg_ini)

def dm_ini(m,i):
    return m.rhov #/(len(list(m.vi_neigh))/m.nv)
mo_agent_2.dm = Param(mo_agent_2.vi_neigh, mutable =True, initialize=dm_ini)

def da_ini(m,i):
    return m.rhov # /(len(list(m.vi_neigh))/m.nv)
mo_agent_2.da = Param(mo_agent_2.vi_neigh, mutable =True, initialize=da_ini)

def dp_ini(m,i):
    return m.rhop # /(len(list(m.vi_neigh))/m.nv)
mo_agent_2.dp = Param(mo_agent_2.ei_neigh, mutable =True, initialize=dp_ini)

def dq_ini(m,i):
    return m.rhop # /(len(list(m.vi_neigh))/m.nv)
mo_agent_2.dq = Param(mo_agent_2.ei_neigh, mutable =True, initialize=dq_ini)

def dpg_ini(m,i,j):
    return m.rhop # /(len(list(m.vi_neigh))/m.nv)
mo_agent_2.dpg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=dpg_ini)

def dqg_ini(m,i,j):
    return m.rhop # /(len(list(m.vi_neigh))/m.nv)
mo_agent_2.dqg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=dqg_ini)

def ym_ini(m,i):
    return 0
mo_agent_2.ym = Param(mo_agent_2.vi_neigh, mutable =True, initialize=ym_ini)

def ya_ini(m,i):
    return 0
mo_agent_2.ya = Param(mo_agent_2.vi_neigh, mutable =True, initialize=ya_ini)

def ypt_ini(m,i):
    return 0
mo_agent_2.ypt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=ypt_ini)

def ypf_ini(m,i):
    return 0
mo_agent_2.ypf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=ypf_ini)

def yqt_ini(m,i):
    return 0
mo_agent_2.yqt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yqt_ini)

def yqf_ini(m,i):
    return 0
mo_agent_2.yqf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yqf_ini)

def ypg_ini(m,i,j):
    return 0
mo_agent_2.ypg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=ypg_ini)

def yqg_ini(m,i,j):
    return 0
mo_agent_2.yqg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=yqg_ini)

def yhm_ini(m,i):
    return 0
mo_agent_2.yhm = Param(mo_agent_2.vi_neigh, mutable =True, initialize=yhm_ini)

def yha_ini(m,i):
    return 0
mo_agent_2.yha = Param(mo_agent_2.vi_neigh, mutable =True, initialize=yha_ini)

def yhpt_ini(m,i):
    return 0
mo_agent_2.yhpt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yhpt_ini)

def yhpf_ini(m,i):
    return 0
mo_agent_2.yhpf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yhpf_ini)

def yhqt_ini(m,i):
    return 0
mo_agent_2.yhqt = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yhqt_ini)

def yhqf_ini(m,i):
    return 0
mo_agent_2.yhqf = Param(mo_agent_2.ei_neigh, mutable =True, initialize=yhqf_ini)

def yhpg_ini(m,i,j):
    return 0
mo_agent_2.yhpg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=yhpg_ini)

def yhqg_ini(m,i,j):
    return 0
mo_agent_2.yhqg = Param(mo_agent_2.gi_neigh, mutable =True, initialize=yhqg_ini)

def active_power_rule(m, i,j):
  i = (i,j)
  if m.pGl[i] == '-Inf' and m.pGu[i] == 'Inf':
    return Constraint.Feasible
  elif m.pGl[i] == '-Inf':
    return m.v_pg[i] <= m.pGu[i]
  elif m.pGu[i] == 'Inf':
    return m.pGl[i] <= m.v_pg[i]
  else: 
    return  (m.pGl[i], m.v_pg[i], m.pGu[i]) 
mo_agent_2.active_power_const = Constraint(mo_agent_2.gi_region | mo_agent_2.gi_neigh, rule=active_power_rule)

def reactive_power_rule(m, i, j):
  i = (i,j)
  if m.qGl[i] == '-Inf' and m.qGu[i] == 'Inf':
    return Constraint.Feasible
  elif m.qGl[i] == '-Inf':
    return m.v_qg[i] <= m.qGu[i]
  elif m.qGu[i] == 'Inf':
    return m.qGl[i] <= m.v_qg[i]
  else: 
    return  (m.qGl[i], m.v_qg[i], m.qGu[i]) 
mo_agent_2.reactive_power_const = Constraint(mo_agent_2.gi_region | mo_agent_2.gi_neigh, rule=reactive_power_rule)

def active_power_eqn(m,i):
  f1 = sum(m.v_pf[j] for j in m.ei if m.fr[j] == i)
  f2 = sum(m.v_pt[j] for j in m.ei if m.to[j] == i)
  if (i,1) in m.gi:
    temp = sum(m.v_pg[j] for j in m.gi if j[0] == i)
    return temp-m.pL[i]-m.gs[i]*m.v_vm[i]**2 == f1+f2
  else:
    return -m.pL[i]-m.gs[i]*m.v_vm[i]**2 == f1+f2
mo_agent_2.active_power_eqn_const = Constraint(mo_agent_2.vi_region, rule=active_power_eqn)

def reactive_power_eqn(m,i):
  f1 = sum(m.v_qf[j] for j in m.ei if m.fr[j] == i)
  f2 = sum(m.v_qt[j] for j in m.ei if m.to[j] == i)
  if (i,1) in m.gi:
    temp = sum(m.v_qg[j] for j in m.gi if j[0] == i)
    return temp-m.qL[i]+m.bs[i]*m.v_vm[i]**2 == f1+f2
  else:
    return -m.qL[i]+m.bs[i]*m.v_vm[i]**2 == f1+f2
mo_agent_2.reactive_power_eqn_const = Constraint(mo_agent_2.vi_region, rule=reactive_power_eqn)

def active_power_fr(m,i):
  gijc = m.Y[i,'Y11r']
  gij = m.Y[i,'Y12r']
  bij = m.Y[i,'Y12i']
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_pf[i] == gijc*m.v_vm[fi]**2 - gij*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[fi]-m.v_va[ti]) + bij*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[fi]-m.v_va[ti])
mo_agent_2.active_power_fr_const = Constraint(mo_agent_2.ei_region, rule=active_power_fr)

def active_power_to(m,i):
  gjic = m.Y[i,'Y22r']
  gji = m.Y[i,'Y21r']
  bji = m.Y[i,'Y21i']
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_pt[i] == gjic*m.v_vm[ti]**2 - gji*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[ti]-m.v_va[fi]) + bji*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[ti]-m.v_va[fi])
mo_agent_2.active_power_to_const = Constraint(mo_agent_2.ei_region, rule=active_power_to)

def reactive_power_fr(m,i):
  bijc = m.Y[i,'Y11i']
  gij = m.Y[i,'Y12r']
  bij = m.Y[i,'Y12i'] 
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_qf[i] == bijc*m.v_vm[fi]**2 -bij*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[fi]-m.v_va[ti]) - gij*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[fi]-m.v_va[ti])
mo_agent_2.reactive_power_fr_const = Constraint(mo_agent_2.ei_region, rule=reactive_power_fr)

def reactive_power_to(m,i):
  bjic = m.Y[i,'Y22i']
  gji = m.Y[i,'Y21r']
  bji = m.Y[i,'Y21i'] 
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_qt[i] == bjic*m.v_vm[ti]**2 -bji*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[ti]-m.v_va[fi]) - gji*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[ti]-m.v_va[fi])
mo_agent_2.reactive_power_to_const = Constraint(mo_agent_2.ei_region, rule=reactive_power_to)

def flow_limit_from_rule(m, i):
  if m.e_su[i] > 0:
    return m.v_pf[i]**2+m.v_qf[i]**2 <= m.e_su[i]**2
  else:
    return Constraint.Feasible
mo_agent_2.flow_limit_from_const = Constraint(mo_agent_2.ei_region, rule=flow_limit_from_rule)

def flow_limit_to_rule(m, i):
  if m.e_su[i] > 0:
    return m.v_pt[i]**2+m.v_qt[i]**2 <= m.e_su[i]**2
  else:
    return Constraint.Feasible
mo_agent_2.flow_limit_to_const = Constraint(mo_agent_2.ei_region, rule=flow_limit_to_rule)

def vol_rule(m,i):
  return (m.vL[i] , m.v_vm[i] , m.vU[i]) 
mo_agent_2.vol_const = Constraint(mo_agent_2.vi_region_all, rule=vol_rule)

def angle_rule(m,i):
  return (m.e_dl[i], m.v_va[m.fr[i]]-m.v_va[m.to[i]],m.e_du[i])
mo_agent_2.angle_const = Constraint(mo_agent_2.ei_region, rule=angle_rule)

def ObjRule(m):
    f1 = 1*sum([m.c2[i]*m.v_pg[i]**2+m.c1[i]*m.v_pg[i]+m.c0[i] for i in m.gi_region])
    f2m = sum([(m.dm[i]/2)*(m.v_vm[i]-m.beta_vm[i])**2 for i in m.vi_neigh])
    f3m = sum([m.ym[i]*(m.v_vm[i]-m.beta_vm[i]) for i in m.vi_neigh])
    f2a = sum([(m.da[i]/2)*(m.v_va[i]-m.beta_va[i])**2 for i in m.vi_neigh])
    f3a = sum([m.ya[i]*(m.v_va[i]-m.beta_va[i]) for i in m.vi_neigh])
    f2pt = sum([(m.dp[i]/2)*(m.v_pt[i]-m.beta_pt[i])**2 for i in m.ei_neigh])
    f2pf = sum([(m.dp[i]/2)*(m.v_pf[i]-m.beta_pf[i])**2 for i in m.ei_neigh])
    f2qt = sum([(m.dq[i]/2)*(m.v_qt[i]-m.beta_qt[i])**2 for i in m.ei_neigh])
    f2qf = sum([(m.dq[i]/2)*(m.v_qf[i]-m.beta_qf[i])**2 for i in m.ei_neigh])
    f2pg = sum([(m.dpg[i]/2)*(m.v_pg[i]-m.beta_pg[i])**2 for i in m.gi_neigh])
    f2qg = sum([(m.dqg[i]/2)*(m.v_qg[i]-m.beta_qg[i])**2 for i in m.gi_neigh])
    f3pt = sum([m.ypt[i]*(m.v_pt[i]-m.beta_pt[i]) for i in m.ei_neigh])
    f3pf = sum([m.ypf[i]*(m.v_pf[i]-m.beta_pf[i]) for i in m.ei_neigh])
    f3qt = sum([m.yqt[i]*(m.v_qt[i]-m.beta_qt[i]) for i in m.ei_neigh])
    f3qf = sum([m.yqf[i]*(m.v_qf[i]-m.beta_qf[i]) for i in m.ei_neigh])
    f3pg = sum([m.ypg[i]*(m.v_pg[i]-m.beta_pg[i]) for i in m.gi_neigh])
    f3qg = sum([m.yqg[i]*(m.v_qg[i]-m.beta_qg[i]) for i in m.gi_neigh])
    return f1+f2m+f3m+f2a+f3a+f2pt+f2pf+f2qt+f2qf+f3pt+f3pf+f3qt+f3qf+f2pg+f2qg+f3pg+f3qg
mo_agent_2.obj = Objective(rule=ObjRule)
