#routines to support in the numerical tests for the routing algorithm
from helpers import BPR, BPR_int
import numpy as np

def init_flows(G):
    for e in G.edges:
        G[e[0]][e[1]]['f']=0
    return G


#sign not necessary (think about what is actually D-1 and about the fact that it has a negative sign in front)
def update_costs(G):
    for e in G.edges:
        x=G[e[0]][e[1]]['f']
        phi=G[e[0]][e[1]]['phi']
        k=G[e[0]][e[1]]['k']
        if 'pot' in G.nodes[e[1]]:
            G[e[0]][e[1]]['cost']=BPR(phi,x,k)-G.nodes[e[1]]['pot']
        else:
            G[e[0]][e[1]]['cost']=BPR(phi,x,k)
    return G



def initEdgeAttr(G,edge_list,k_list,phi_list,is_negative):

    for i in range(len(edge_list)):
        e=edge_list[i]
        G[e[0]][e[1]]['k']=k_list[i]
        G[e[0]][e[1]]['phi']=phi_list[i]
        G[e[0]][e[1]]['sign']=(-1)**is_negative[i]
    
    return G

def initPotShifts(G,nodes_pots):
    for (n,p) in nodes_pots:
        G.nodes[n]['pot']=p
    return G

def vect_attribute(G,att):
    # x=[]
    x=dict()
    for e in G.edges():
        # x.append(G[e[0]][e[1]][att])
        x[e]=G[e[0]][e[1]][att]
    return x

