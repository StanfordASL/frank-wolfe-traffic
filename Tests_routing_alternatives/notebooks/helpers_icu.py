import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import networkx as nx

def BPR_int(phi,x,kappa,alpha=0.15,beta=4):
    return phi*(x+alpha/(beta+1)*cp.power(x,(beta+1))/np.power(kappa,beta))

#returns the value of BPR int, not just an expression as is the case above
def BPR_int_val(phi,x,kappa,alpha=0.15,beta=4):
    return phi*(x+alpha/(beta+1)*np.power(x,(beta+1))/np.power(kappa,beta))

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

def disp_costs(G):

    print("Cost for dummy edge: ", G['1']['2_p']['cost'])
    print("Cost for normal edge (1,2),(2,2_p): ", G['1']['2']['cost'], " --- ",G['2']['2_p']['cost'])
    print("Cost for normal edge (1,2,2_p): ", G['1']['2']['cost']+G['2']['2_p']['cost'])
    print("Cost for 1-R: ", G['1']['R']['cost'])
    print("Cost for 2-R: ", G['2']['R']['cost'])

def print_final_flows(G_k):
    G_end=G_k[-1]
    for e in G_end.edges():
        print(e," : ",G_end[e[0]][e[1]]['f_m']+G_end[e[0]][e[1]]['f_r'])

def print_final_cost(G_k):
    G_end=G_k[-1]
    for e in G_end.edges():
        print(e," : ",G_end[e[0]][e[1]]['cost'])

def sanity_check_N(G_k):
    G_end=G_k[-1]
    print("1 to 2: ", G_end['1']['2_p']['f_m']+G_end['1']['2']['f_m'])
    print("2 to 1: ", G_end['2']['1_p']['f_m']+G_end['2']['1']['f_m'])
            
            
def sanity_check_cost(G_k):
    G_end=G_k[-1]
    print("1 to 2: ", G_end['1']['2_p']['cost']," ===== ", G_end['1']['2']['cost']+G_end['2']['2_p']['cost'])
    print("2 to 1: ", G_end['2']['1_p']['cost']," ===== ", G_end['2']['1']['cost']+G_end['1']['1_p']['cost'])
        
    
def plot_sanity_check_cost(G_k,scale='log'):
    c_12=[]
    c_21=[]
    for G in G_k:
        c_12.append(np.abs((G['1']['2']['cost']+G['2']['2_p']['cost']-G['1']['2_p']['cost'])/G['1']['2_p']['cost']))
        c_21.append(np.abs((G['2']['1']['cost']+G['1']['1_p']['cost']-G['2']['1_p']['cost'])/G['2']['1_p']['cost']))

    plt.figure(figsize=(13,5))
    plt.plot(c_12,'--o',label="(1,2)")
    plt.plot(c_21,'--o',label="(2,1)")
    plt.legend()
    plt.grid(True)
    plt.yscale(scale)
    plt.xlabel("Iteration #")
    plt.ylabel("error")


def analyze_cost_oscillations(G_k,o,d,scale='log'):
    c=[]
    for G in G_k:
        c.append(((G[o][d]['cost']+G[d][d+'_p']['cost']-G[o][d+'_p']['cost'])/G[o][d+'_p']['cost']))
    c=np.array(c)
    plt.figure(figsize=(13,5))
    plt.plot(c,'o',label="("+o+","+d+")")
    plt.plot(-c,'ro')
    plt.legend()
    plt.grid(True)
    plt.yscale(scale)
    plt.xlabel("Iteration #")
    plt.ylabel("error")