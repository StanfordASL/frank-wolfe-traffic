import numpy as np
import cvxpy as cp

def BPR_int(phi,x,kappa,alpha=0.15,beta=4):
    return phi*(x+alpha/(beta+1)*np.power(x,(beta+1))/np.power(kappa,beta))


# def BPR_int(phi,x,kappa,alpha=0.15,beta=4):
#     l=x.shape[0]
#     INT_COST=0
#     for i in range(l):
#         INT_COST+=phi[i]*(x[i]+alpha/(beta+1)*np.power(x[i],(beta+1))/np.power(kappa[i],beta))
#     return INT_COST

def BPR(phi,x,kappa,alpha=0.15,beta=4):
    return phi*(1+alpha*(np.divide(x,kappa))**beta)

def phi(l,t):
    return 36*l/t

def cost_per_edge(alpha,beta,phi_vec,flow_vec,kappa_vec,K_vec):
    c=BPR(alpha,beta,phi_vec,flow_vec,kappa_vec)-K_vec
    return c