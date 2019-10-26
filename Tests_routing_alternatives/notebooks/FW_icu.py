import numpy as np 
import cvxpy as cp
import networkx as nx
from routines_icu import *
from helpers_icu import *



def init_flows(G):
    for e in G.edges:
        if e[1]=='R':
            G[e[0]][e[1]]['f_m']=0
            G[e[0]][e[1]]['f_r']=G[e[0]][e[1]]['k']
        else:
            G[e[0]][e[1]]['f_m']=0
            G[e[0]][e[1]]['f_r']=0 
    return G

def AoN(G,OD):
    
    y_k=init_y(G)
    eps=10**-6
    for (o,d) in OD.keys():
        N=OD[o,d]
        if N>eps:
            if d=='R':
                flag='f_r'
            else:
                flag='f_m'

            print("AON, (o,d):", o,d)
            path=nx.shortest_path(G,source=o,target=d,weight='cost')
            
            
            for i in range(len(path)-1):
                y_k[(path[i],path[i+1]),flag]+=N
            
    return y_k

def FW(G_0,OD,edge_list,dummy_nodes):
    
    G_list=[]
    G_list.append(G_0)
    i=1
    G_k=G_0
    while i<10:  
        print("###############################") 
        print("iteration # ", i)
        G_crt=G_k.copy()
        ri_k,G_crt=estimate_ri_k(G_crt,dummy_nodes)
        OD=update_OD(OD,ri_k,dummy_nodes)
        G_crt=update_capacities(G_crt,ri_k)
        print("oD:", OD)
        y_k=AoN(G_crt,OD)
        # print("Y_K:", y_k)
        a_k=line_search(G_crt,y_k,edge_list)
        print('ALPHA:', a_k)
        G_crt=update_flows(G_crt,y_k,a_k,edge_list)
        G_crt=update_costs(G_crt, 80)
        G_k=G_crt
        G_list.append(G_k)
        i+=1
    return G_list,OD

def estimate_ri_k(G,dummy_nodes):
    ri_k=dict()
    for n in G.nodes():
        ri_k[n]=0
    for e in G.edges():
        if e[1] not in dummy_nodes.keys():
            ri_k[e[0]]+=G[e[0]][e[1]]['f_m']
            ri_k[e[1]]-=G[e[0]][e[1]]['f_m']
    for n in G.nodes():
        G.nodes[n]["ri"]=ri_k[n]
    return ri_k,G

def update_OD(OD,ri_k, dummy_nodes):
    eps=10**-6
    for n in ri_k.keys():
        if not n=='R':
            if ri_k[n]<-eps: #you are in excess
                # n=dummy_nodes[n]
                OD[(n,'R')]=-ri_k[n]
            else:
                OD[(n,'R')]=0
    return OD

def update_capacities(G,ri_k): #not so sure whether we need to update costs as well
    eps=10**-6
    for n in G.nodes():
        if ri_k[n]>eps: #this guy needs rebalancers to go through it
            G[n]['R']['k']=ri_k[n]
            
    return G

def init_y(G):
    y_k=dict()
    for e in G.edges():
        y_k[e,'f_r']=0
        y_k[e,'f_m']=0
    return y_k

#computes the total cost under the current model
def Total_Cost(G,y_k,a_k,edge_list):
    F_E=0
    for i in range(len(edge_list)):#you know for sure exactly what edge it is for
        e=edge_list[i]
        x_k_e=G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_r'] # retrieve the flow, total
        y_k_e=y_k[(e[0],e[1]),'f_m']+y_k[(e[0],e[1]),'f_r'] #retrieve the flow from the manual assignment, total
        
        flow_tmp=x_k_e+a_k*(y_k_e-x_k_e)
        
        #retrieve parameters to compute the BPR
        phi=G[e[0]][e[1]]['phi']
        k=G[e[0]][e[1]]['k']
        sign=G[e[0]][e[1]]['sign']
        
        F_E+=BPR_int(phi,flow_tmp,k,alpha=0.15,beta=4)#I am assuming there will be syntaxic problems there
        
        
        #this has to be included because it is directly included in the definition of the cost function
        if G[e[0]][e[1]]['sign']==(-1): #we have a negative edge
            F_E-=flow_tmp*80#INVERSE_DEMAND_SHIFT
        

        # if 'pot' in G.nodes[e[1]]:
        #     F_E+=G.nodes[e[1]]['pot']*flow_tmp
            
    return F_E

def line_search(G,y_k,edge_list):
    a_k=cp.Variable()
    constraints=[a_k>=0, a_k<=1]
    obj=Total_Cost(G,y_k,a_k,edge_list)
    # print(obj)
    prob=cp.Problem(cp.Minimize(obj),constraints)
    prob.solve(verbose=False)
    print(prob.status)
    return a_k.value

def update_flows(G,y_k,a_k,edge_list):
    
    for i in range(len(edge_list)):
        e=edge_list[i]
        for flag in ['f_m', 'f_r']:
            x_k_e=G[e[0]][e[1]][flag] # retrieve the flow
            y_k_e=y_k[(e[0],e[1]),flag] #retrieve the flow from the manual assignment
            G[e[0]][e[1]][flag]+=a_k*(y_k_e-x_k_e)
        
    G=update_costs(G,80)#inverse demand shift replaced by nominal value currently
    return G   