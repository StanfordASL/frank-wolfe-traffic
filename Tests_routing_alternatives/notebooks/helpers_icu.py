import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import networkx as nx

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

def plot_edge_attrs(G_list,attrs):
    G_=G_list[0]
    f,axes=plt.subplots(len(G_.edges()),len(attrs),figsize=(18,5*len(G_.edges())))
    i=0
    for e in G_.edges():
        for j in range(len(attrs)):
            att=[]
            for G in G_list: 
                att.append(G[e[0]][e[1]][attrs[j]])
            axes[i,j].plot(att,'--o',label=attrs[j])
            axes[i,j].grid(True)
            axes[i,j].set_xlabel('Iteration #')
            axes[i,j].set_title(' Edge : ' + str(e))
            axes[i,j].legend()
        i+=1

def plot_node_attrs(G_list,attrs):
    G_=G_list[0]
    f,axes=plt.subplots(len(G_.nodes()),len(attrs),figsize=(18,5*len(G_.nodes())))
    i=0
    for n in G_.nodes():
        for j in range(len(attrs)):
            att=[]
            for G in G_list: 
                att.append(G.nodes[n][attrs[j]])
            
            if len(attrs)==1:
                axes[i].plot(att,'--o',label=attrs[j])
                axes[i].grid(True)
                axes[i].set_xlabel('Iteration #')
                axes[i].set_title(' node : ' + str(n))
            
            else:
                axes[i,j].plot(att,'--o',label=attrs[j])
                axes[i,j].grid(True)
                axes[i,j].set_xlabel('Iteration #')
                axes[i,j].set_title(' node : ' + str(n))
                axes[i,j].legend()
        i+=1