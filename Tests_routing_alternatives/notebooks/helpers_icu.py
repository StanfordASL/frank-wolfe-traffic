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

def plot_edge_attrs(G_list,y_list,attrs,dots=True,lims=None):
    G_=G_list[0]
    _,axes=plt.subplots(len(G_.edges()),len(attrs),figsize=(20,5*len(G_.edges())))
    i=0


    if dots:
        options='--o'
    else:
        options='--'

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
            axes[i,j].plot(att,options,label=attrs[j])
            axes[i,j].grid(True)
            axes[i,j].set_xlabel('Iteration #')
            axes[i,j].set_title(' Edge : ' + str(e))
            axes[i,j].legend()
            axes[i,j].set_xlim(lims)
        i+=1

def plot_node_attrs(G_list,attrs):
    G_=G_list[0]
    _,axes=plt.subplots(len(G_.nodes()),len(attrs),figsize=(18,5*len(G_.nodes())))
    i=0
    for n in G_.nodes():
        for j in range(len(attrs)):
            att=[]
            for G in G_list: 
                att.append(G.nodes[n][attrs[j]])
            
            if len(attrs)==1:
                axes[i].plot(att,'--',label=attrs[j])
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
        
    
def plot_errors(G_k,scale='log',dots=True,lims=None):
    c_12=[]
    c_21=[]
    for G in G_k:
        c_12.append(np.abs((G['1']['2']['cost']+G['2']['2_p']['cost']-G['1']['2_p']['cost'])/G['1']['2_p']['cost']))
        c_21.append(np.abs((G['2']['1']['cost']+G['1']['1_p']['cost']-G['2']['1_p']['cost'])/G['2']['1_p']['cost']))

    if dots:
        options='--o'
    else:
        options='--'

    _,axes=plt.subplots(2,1,figsize=(13,10))
    axes[0].plot(c_12,options,label="(1,2)")
    axes[1].plot(c_21,options,label="(2,1)")
    for i in [0,1]:
        axes[i].legend()
        axes[i].grid(True)
        axes[i].set_yscale(scale)
        axes[i].set_ylabel("relative error")
        axes[i].set_xlabel("Iteration #")
        axes[i].set_xlim(lims)


def analyze_cost_oscillations(G_k,o,d,lims=None,scale='log'):
    c=[]
    f_m=[]
    f_r=[]
    for G in G_k:
        c.append(((G[o][d]['cost']+G[d][d+'_p']['cost']-G[o][d+'_p']['cost'])/G[o][d+'_p']['cost']))
        f_m.append(G[o][d]['f_m'])
        f_r.append(G[o][d]['f_r'])
    c=np.abs(np.array(c))
    f_m=np.array(f_m)
    f_r=np.array(f_r)
    fig, ax1 = plt.subplots(figsize=(18,5))
    ax1.plot(c,'ro',label="("+o+","+d+")")
    # ax1.plot(-c,'ro')
    ax2=ax1.twinx()
    ax2.plot(f_m,"--o",label="f_m")
    ax2.plot(f_r,"--o",label="f_r")
    ax2.plot(f_m+f_r,"--o",label="sum")
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("error")
    ax2.set_ylabel("Flow")
    if not lims == None:
        ax1.set_xlim(lims)


def analyze_cost_oscillations_2(G_k,o,d,lims=None,scale='linear'):
    tgt_cost=90
    tgt_flow=8.65
    c=[]
    f_m=[]
    f_r=[]
    for G in G_k:
        c.append(G[o][d]['cost']+G[d][d+'_p']['cost'])
        f_m.append(G[o][d]['f_m'])
        f_r.append(G[o][d]['f_r'])
    c=np.abs(np.array(c))

    #drop the first component
    c=c[1:]

    f_m=np.array(f_m)
    f_r=np.array(f_r)
    fig, ax1 = plt.subplots(figsize=(18,5))
    ax1.plot(c,'ro--',label="cost")
    ax1.plot(np.linspace(1,c.shape[0],50),tgt_cost*np.ones(50),'g--',label='Target value')
    # ax1.plot(-c,'ro')
    ax2=ax1.twinx()
    # ax2.plot(f_m,"--o",label="f_m")
    # ax2.plot(f_r,"--o",label="f_r")
    ax2.plot(f_m+f_r,"--o",label="total flow (m + r)")
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("Total cost")
    ax2.set_ylabel("Flow")
    ax1.set_ylim([88,92])
    ax2.set_ylim([tgt_flow*0.98,tgt_flow*1.02])

    if not lims == None:
        ax1.set_xlim(lims)


def analyze_cost_oscillations_3(G_k,y_k,o,d,lims=None,scale='linear'):
    #same analysis as above but focusing on the dummy edges
    
    tgt_cost=90
    tgt_flow=8.65
    c=[]
    f_m=[]
    f_r=[]
    y_dummy=[]
    for G in G_k:
        c.append(G[o][d]['cost']+G[d][d+'_p']['cost'])
        f_m.append(G[o][d+"_p"]['f_m'])
        f_r.append(G[o][d+"_p"]['f_r'])
    for y in y_k:
        y_dummy.append(y[(o,d+"_p"),'f_m'])

    #drop the first component
    c=c[1:]

    c=np.abs(np.array(c))
    f_m=np.array(f_m)
    f_r=np.array(f_r)
    y_dummy=np.array(y_dummy)
    fig, ax1 = plt.subplots(figsize=(18,5))
    # ax1.plot(-c,'ro')
    ax2=ax1.twinx()
    for k in range(y_dummy.shape[0]):
        if y_dummy[k]!=0:
            ax2.plot(k*np.ones(50),np.linspace(1,10,50),'k-',linewidth=1.3)
    ax1.plot(c,'ro',label="Cost of the normal edge")
    ax1.plot(np.linspace(1,c.shape[0],50),tgt_cost*np.ones(50),'g--',label='Target value')
    ax2.plot(f_m,"--o",label="f_m dummy edge")
    ax2.plot(f_r,"--o",label="f_r dummy edge")
    # ax2.plot(y_dummy,label='y_m')
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("Total cost")
    ax2.set_ylabel("flow")
    ax1.set_ylim([88,92])
    avg_=np.mean(f_m[-7:-1]+f_r[-7:-1])
    ax2.set_ylim([(avg_)*0.95,(avg_)*1.05])
    if not lims ==None:
        ax1.set_xlim(lims)
        
    print("Black lines indicate that the assignment variable y_m is activated for the dummy edge")