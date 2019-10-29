import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import networkx as nx

def BPR_int(phi,x,kappa,alpha=0.15,beta=4):
    return phi*(x+alpha/(beta+1)*cp.power(x,(beta+1))/np.power(kappa,beta))

#returns the value of BPR int, not just an expression as is the case above
def BPR_int_val(phi,x,kappa,alpha=0.15,beta=4):
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

def plot_edge_attrs(G_list,y_list,attrs):
    G_=G_list[0]
    _,axes=plt.subplots(len(G_.edges()),len(attrs),figsize=(20,5*len(G_.edges())))
    i=0
    for e in G_.edges():
        for j in range(len(attrs)):
            att=[]
            if attrs[j]=="y_m":
                for y_k in y_list:
                    att.append(y_k[e,'f_m'])
            elif attrs[j]=="y_r":
                for y_k in y_list:
                    att.append(y_k[e,'f_r']) 
            else:
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

####################################################################
########### DEBUG HELPERS 
####################################################################

def disp_both_costs(G):

    print("Cost for dummy edge: ", G['1']['2_p']['cost'])
    print("Cost for normal edge (1,2),(2,2_p): ", G['1']['2']['cost'], " --- ",G['2']['2_p']['cost'])
    print("Cost for normal edge (1,2,2_p): ", G['1']['2']['cost']+G['2']['2_p']['cost'])

def print_final_flows(G_k):
    G_end=G_k[-1]
    for e in G_end.edges():
        print(e," : ",G_end[e[0]][e[1]]['f_m']+G_end[e[0]][e[1]]['f_r'])