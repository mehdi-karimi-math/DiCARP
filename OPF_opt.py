# This file creates a Pyomo Abstract model for solving the OPF for a given problem instance. 
#
# Copyright (c) 2023, by 
# Mehdi Karimi

from pyomo.environ import *
mo = AbstractModel()

mo.baseMVA = Param()
mo.nv = Param(within=NonNegativeIntegers)
mo.ne = Param(within=NonNegativeIntegers)
mo.ng = Param(within=NonNegativeIntegers)

mo.vi = RangeSet(1,mo.nv)
mo.ei = RangeSet(1,mo.ne)
mo.gi = Set(within=mo.vi * NonNegativeIntegers)

mo.slack = Param(within=mo.vi)
mo.angslack = Param()

mo.pL = Param(mo.vi, within=Reals)
mo.qL = Param(mo.vi, within=Reals)

mo.vL = Param(mo.vi, within=Reals)
mo.vU = Param(mo.vi, within=Reals)

mo.gs = Param(mo.vi, within=Reals)
mo.bs = Param(mo.vi, within=Reals)

mo.fr = Param(mo.ei, within = mo.vi)
mo.to = Param(mo.ei, within = mo.vi)

mo.Y_ind = Set(initialize=['Y11r', 'Y11i', 'Y12r', 'Y12i', 'Y21r', 'Y21i', 'Y22r', 'Y22i'])
mo.Y = Param(mo.ei, mo.Y_ind, within = Reals)
mo.e_con = Param(mo.ei, within = Reals)
mo.e_sus = Param(mo.ei, within = Reals)
mo.e_bsh = Param(mo.ei, within = Reals)

mo.e_su = Param(mo.ei, within= Reals)
mo.e_dl = Param(mo.ei, within= Reals)
mo.e_du = Param(mo.ei, within= Reals)


mo.pGl = Param(mo.gi)
mo.pGu = Param(mo.gi)
mo.qGl = Param(mo.gi)
mo.qGu = Param(mo.gi)

mo.c2 = Param(mo.gi, within = NonNegativeReals)
mo.c1 = Param(mo.gi, within = NonNegativeReals)
mo.c0 = Param(mo.gi, within = NonNegativeReals)

mo.v_pg = Var(mo.gi)
mo.v_qg = Var(mo.gi)
mo.v_vm = Var(mo.vi)
mo.v_va = Var(mo.vi)
mo.v_pt = Var(mo.ei)
mo.v_pf = Var(mo.ei)
mo.v_qt = Var(mo.ei)
mo.v_qf = Var(mo.ei)


def active_power_rule(m, i, j):
  i= (i,j)
  if m.pGl[i] == '-Inf' and m.pGu[i] == 'Inf':
    return Constraint.Feasible
  elif m.pGl[i] == '-Inf':
    return m.v_pg[i] <= m.pGu[i]
  elif m.pGu[i] == 'Inf':
    return m.pGl[i] <= m.v_pg[i]
  else: 
    return  (m.pGl[i], m.v_pg[i], m.pGu[i]) 
mo.active_power_const = Constraint(mo.gi, rule=active_power_rule)

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
mo.reactive_power_const = Constraint(mo.gi, rule=reactive_power_rule)

def active_power_eqn(m,i):
  f1 = sum(m.v_pf[j] for j in m.ei if m.fr[j] == i)
  f2 = sum(m.v_pt[j] for j in m.ei if m.to[j] == i)
  if (i,1) in m.gi:
     temp = sum(m.v_pg[j] for j in m.gi if j[0] == i)
     return temp-m.pL[i]-m.gs[i]*m.v_vm[i]**2  == f1+f2
  else:
    return -m.pL[i]-m.gs[i]*m.v_vm[i]**2 == f1+f2
mo.active_power_eqn_const = Constraint(mo.vi, rule=active_power_eqn)

def reactive_power_eqn(m,i):
  f1 = sum(m.v_qf[j] for j in m.ei if m.fr[j] == i)
  f2 = sum(m.v_qt[j] for j in m.ei if m.to[j] == i)
  if (i,1) in m.gi:
    temp = sum(m.v_qg[j] for j in m.gi if j[0] == i)
    return temp-m.qL[i]+m.bs[i]*m.v_vm[i]**2 == f1+f2
  else:
    return -m.qL[i]+m.bs[i]*m.v_vm[i]**2 == f1+f2
mo.reactive_power_eqn_const = Constraint(mo.vi, rule=reactive_power_eqn)

def active_power_fr(m,i):
  gijc = m.Y[i,'Y11r']
  gij = m.Y[i,'Y12r']
  bij = m.Y[i,'Y12i']
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_pf[i] == gijc*m.v_vm[fi]**2 - gij*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[fi]-m.v_va[ti]) + bij*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[fi]-m.v_va[ti])
mo.active_power_fr_const = Constraint(mo.ei, rule=active_power_fr)

def active_power_to(m,i):
  gjic = m.Y[i,'Y22r']
  gji = m.Y[i,'Y21r']
  bji = m.Y[i,'Y21i']
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_pt[i] == gjic*m.v_vm[ti]**2 - gji*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[ti]-m.v_va[fi]) + bji*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[ti]-m.v_va[fi])
mo.active_power_to_const = Constraint(mo.ei, rule=active_power_to)

def reactive_power_fr(m,i):
  bijc = m.Y[i,'Y11i']
  gij = m.Y[i,'Y12r']
  bij = m.Y[i,'Y12i'] 
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_qf[i] == bijc*m.v_vm[fi]**2 -bij*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[fi]-m.v_va[ti]) - gij*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[fi]-m.v_va[ti])
mo.reactive_power_fr_const = Constraint(mo.ei, rule=reactive_power_fr)

def reactive_power_to(m,i):
  bjic = m.Y[i,'Y22i']
  gji = m.Y[i,'Y21r']
  bji = m.Y[i,'Y21i'] 
  fi = m.fr[i]
  ti = m.to[i]
  return m.v_qt[i] == bjic*m.v_vm[ti]**2 -bji*m.v_vm[fi]*m.v_vm[ti]*cos(m.v_va[ti]-m.v_va[fi]) - gji*m.v_vm[fi]*m.v_vm[ti]*sin(m.v_va[ti]-m.v_va[fi])
mo.reactive_power_to_const = Constraint(mo.ei, rule=reactive_power_to)

def flow_limit_from_rule(m, i):
  if m.e_su[i] > 0:
    return m.v_pf[i]**2+m.v_qf[i]**2 <= m.e_su[i]**2
  else:
    return Constraint.Feasible
mo.flow_limit_from_const = Constraint(mo.ei, rule=flow_limit_from_rule)

def flow_limit_to_rule(m, i):
  if m.e_su[i] > 0:
    return m.v_pt[i]**2+m.v_qt[i]**2 <= m.e_su[i]**2
  else:
    return Constraint.Feasible
mo.flow_limit_to_const = Constraint(mo.ei, rule=flow_limit_to_rule)

def vol_rule(m,i):
  return (m.vL[i] , m.v_vm[i] , m.vU[i]) 
mo.vol_const = Constraint(mo.vi, rule=vol_rule)

def angle_rule(m,i):
  return (m.e_dl[i], m.v_va[m.fr[i]]-m.v_va[m.to[i]],m.e_du[i])
mo.angle_const = Constraint(mo.ei, rule=angle_rule)

def ObjRule(m):
  return sum([m.c2[i]*m.v_pg[i]**2+m.c1[i]*m.v_pg[i]+m.c0[i] for i in m.gi])
mo.obj = Objective(rule=ObjRule)

# instance=mo.create_instance("case9.dat")

# ipopt_res = SolverFactory('ipopt').solve(instance, tee=True)

# print(ipopt_res.Solver)

# print(value(instance.obj))