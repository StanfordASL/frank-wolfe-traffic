#This document contains the implementation of Kiril's algorithm to solve the assignment 
#problem with rebalancers with graph extension, with Elastic Demand too

import numpy as np 
import cvxpy as cp
import networkx as nx
from FW_icu import update_OD,update_capacities, AoN, estimate_ri_k
from helpers_icu import Value_Total_Cost

def solve(G_0,OD,edge_list,dummy_nodes,tol=10**-6):
    #while ri too different from d do: 
    #step 1: estimate ri
    #step 2: solve exactly the elastic demand problem as a function of ri
    #end while

    #Variables to store at each iterations
    i=1
    G_=[]
    ri_=[]


    #initialize certain values
    G_k=G_0
    ri_k=0

    #initialize flows and compute ri? 
    G_crt=init_flows(G_crt,OD)
    #estimate ri

    compute=True
    while compute:
        
        #solve system entirely for the given parameters
        G_list,y_list,opt_res,OD_list=FW_graph_extension(
            G_k,OD,edge_list,dummy_nodes,ri_k,tol=10**-6, 
            step='line_search', evolving_bounds=False)

        #solution at previous step
        #important because it basically contains the passenger flows
        #which directly informs the value of the rebalancing
        G_end=Glist[-1]

        #estimate #ri_k, update OD, assign, update costs
        ri_new,G_end=estimate_ri_k(G_end,dummy_nodes, ri_smoothing=False, a_k=0)
        
        if diff_ri(ri_k,ri_new)<tol:
            compute=False
        
        #The big question is what value of the flows and the graph you actually restart everything with
        #all the time? 

        #update the values for the new iteration

        ri_k = ri_new
        G_k = 

        #Save the different variables
        G_.append(G_end) #we keep the end version of the solved graph
        ri_.append(ri_new)
    return



def FW_graph_extension(
    G_0,OD,edge_list,dummy_nodes,ri_k,tol=10**-6, 
    step='line_search', evolving_bounds=True):
    #python implementation of Kiril's routine
    #ri_t are the estimate of ri at timestep k

    ###########################################
    # PARAMETERS
    ##########################################  
    y_list=[]
    opt_res=dict()
    opt_res['a_k']=[]
    opt_res['obj']=[]
    opt_res['dual_gap']=[]
    OD_list=[]
    G_list=[]
    G_list.append(G_0)
    G_k=G_0
    a_k=1
    i=1
    compute=True

    #################################
    # update the OD pairs and capacities
    #################################
    G_crt=G_k.copy()
    OD=update_OD(OD,ri_k,a_k,dummy_nodes, G_crt, evolving_bounds)
    #you update capacities because you have new values of ri_k
    G_crt=update_capacities(G_crt,ri_k, dummy_nodes)

    ###################################
    # Reinitialize
    ###################################
    #I think we need a new initialization at every step
    #because the problem is never the same 

    G_crt=init_flows(G_crt,OD)

    ###################################
    # Solve for the given ri_k
    ###################################

    while compute: #introduce stopping criterion for FW algorithm. See Jaggi.  
        #here we perform the assignment with OD remaining fixed
        #the structure of the graph is the same
        #we just want to reach the "perfect" assignment
        
        #perform AON assignment
        y_k=AoN(G_crt,OD,dummy_nodes) #check if there is a special feature of AON as previously setup that prevents its good implementation
        if step == 'line_search':
            # a_k,obj_k=line_search(G_crt,y_k,edge_list)#include the fixed step size,
            pass
        elif step == 'fixed':
            a_k,obj_k=fixed_step(G_crt,y_k,edge_list,i)
        else:
            print("wrong optim step chosen")
            return

        #compute the duality gap
        #based on the current version of y_k and x_k
        duality_gap=compute_duality_gap(G_crt, y_k)
        if duality_gap<tol:
            compute=False

        #update the flows
        G_crt=update_flows(G_crt,y_k,a_k,edge_list)

        #save for analyses
        opt_res['obj'].append(obj_k)
        opt_res['a_k'].append(a_k)
        opt_res['dual_gap'].append(duality_gap)
        G_k=G_crt
        G_list.append(G_k)
        y_list.append(y_k)
        OD_list.append(OD.copy())
        i+=1
        
    return G_list,y_list,opt_res,OD_list

def compute_duality_gap(G_k,y_k):
    #G_k : version of the graph (and therefore the flows) at iteration k
    #y_k : minimizer of linearized problem at iteration k 
    d_gap=0
    for e in G_k.edges():
        for flag in ['f_m','f_r']: 
            x_k_ij=G_k[e[0]][e[1]][flag]
            c_k_ij=G_k[e[0]][e[1]]['cost']
            y_k_ij=y_k[(e[0],e[1]),flag]
            d_gap+=(x_k_ij-y_k_ij)*c_k_ij

    return d_gap

###########################################################################
# Original version of functions in FW_icu
#
# Those functions are already present in the FW_icu file. We need another
# version because here we want to update all flows at the same time
###########################################################################

#TODO: consolidate them all in a single file, with flexible versions

def fixed_step(G,y_k,edge_list,k):
    #here we actually update both the rebalancing and the passenger flows
    #we need to take both into account because we are solving for a fixed value of the ri's

    gamma=2/(k+2)

    for i in range(len(edge_list)):
        e=edge_list[i]
        for flag in ['f_m','f_r']:
            x_k_e=G[e[0]][e[1]][flag] # retrieve the flow
            y_k_e=y_k[(e[0],e[1]),flag] #retrieve the flow from the manual assignment
            G[e[0]][e[1]][flag]=(1-gamma)*x_k_e+gamma*y_k_e
        
    return gamma, Value_Total_Cost(G)


def update_flows(G,y_k,a_k,edge_list):
    
    for i in range(len(edge_list)):
        e=edge_list[i]
        for flag in ['f_m', 'f_r']:
            x_k_e=G[e[0]][e[1]][flag] # retrieve the flow
            y_k_e=y_k[(e[0],e[1]),flag] #retrieve the flow from the manual assignment
            G[e[0]][e[1]][flag]=(1-a_k)*x_k_e + a_k * y_k_e
    return G   

def init_flows(G,OD):
    #Initiliaze the flows, with a feasible solution
    #still unclear whether we need to initialize the rebalancers appropriately too
    for (o,d) in OD.keys():
        path=nx.shortest_path(G,source=o,target=d,weight='cost')
        for i in range(len(path)-1):
            if d=='R':
                G[path[i]][path[i+1]]['f_r']+=OD[o,d]
            else:
                G[path[i]][path[i+1]]['f_m']+=OD[o,d]
    return G