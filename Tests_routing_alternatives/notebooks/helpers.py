import numpy as np


def BPR_int(alpha,beta,phi,x,kappa):
    l=x.shape[0]
    INT_COST=0
    for i in range(l):
        INT_COST+=phi[i]*(x[i]+alpha/(beta+1)*np.power(x[i],(beta+1))/np.power(kappa[i],beta))
    return INT_COST

def BPR(alpha,beta,phi,x,kappa):
    return phi*(1+alpha*(np.divide(x,kappa))**beta)

def phi(l,t):
    return 36*l/t

def cost_per_edge(alpha,beta,phi_vec,flow_vec,kappa_vec,K_vec):
    c=BPR(alpha,beta,phi_vec,flow_vec,kappa_vec)-K_vec
    return c

def TotalCost(alpha,beta,phi_vec,flow_vec,kappa_vec,K_vec):
    return cp.sum(BPR_int(alpha,beta,phi_vec,flow_vec,kappa_vec)-K_vec*flow_vec)