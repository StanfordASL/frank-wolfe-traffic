#routines to support in the numerical tests for the routing algorithm
from helpers_icu import BPR, BPR_int
import numpy as np

 


#sign not necessary (think about what is actually D-1 and about the fact that it has a negative sign in front)
def update_costs(G,INVERSE_DEMAND_SHIFT):
    for e in G.edges:
        x=G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_m']
        phi=G[e[0]][e[1]]['phi']
        k=G[e[0]][e[1]]['k']
        G[e[0]][e[1]]['cost']=BPR(phi,x,k)

        if G[e[0]][e[1]]['sign']==(-1): #we have a negative edge
            G[e[0]][e[1]]['cost']-=INVERSE_DEMAND_SHIFT
            
        if 'pot' in G.nodes[e[1]]:
            G[e[0]][e[1]]['cost']+=G.nodes[e[1]]['pot']

    return G



def initEdgeAttr(G,edge_list,k_list,phi_list,is_negative):

    for i in range(len(edge_list)):
        e=edge_list[i]
        G[e[0]][e[1]]['k']=k_list[i]
        G[e[0]][e[1]]['phi']=phi_list[i]
        G[e[0]][e[1]]['sign']=(-1)**is_negative[i]
    
    return G

def initNodeAttr(G,nodes_pots):

    #Node potential
    for (n,p) in nodes_pots:
        G.nodes[n]['pot']=p
    
    #ri_k
    for n in G.nodes():
        G.nodes[n]['ri']=0

    return G

def vect_attribute(G,att):
    # x=[]
    x=dict()
    for e in G.edges():
        # x.append(G[e[0]][e[1]][att])
        x[e]=G[e[0]][e[1]][att]
    return x

