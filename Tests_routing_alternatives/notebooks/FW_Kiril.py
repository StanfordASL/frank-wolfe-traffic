#This document contains the implementation of Kiril's algorithm to solve the assignment 
#problem with rebalancers with graph extension, with Elastic Demand too

import numpy as np 
import cvxpy as cp
import networkx as nx


def solve():
    #while ri too different from d do: 
    #step 1: estimate ri
    #step 2: solve exactly the elastic demand problem as a function of ri
    #end while
    return


def FW_graph_extention(G_0,OD,edge_list,dummy_nodes,ri_k,maxIter=50, step='line_search', evolving_bounds=True):
    #python implementation of Kiril's routine
    #ri_t are the estimate of ri at timestep k

    ###########################################
    # PARAMETERS
    ##########################################  
    y_list=[]
    opt_res=dict()
    opt_res['a_k']=[]
    opt_res['obj']=[]
    OD_list=[]
    G_list=[]
    G_list.append(G_0)
    G_k=G_0
    a_k=1
    i=1
    ############################################
    # FW loop
    ############################################

    while i<maxIter: #introduce stopping criterion for FW algorithm. See Jaggi.  
        G_crt=G_k.copy()

        OD=update_OD(OD,ri_k,a_k,dummy_nodes, G_crt, evolving_bounds)
        G_crt=update_capacities(G_crt,ri_k, dummy_nodes)

        G_crt=assign_rebalancers(G_crt,OD,rebalancer_smoothing,a_k)
        G_crt=update_costs(G_crt,80)
        #disp_costs(G_crt)#debug helper
        
        #perform AON assignment
        y_k=AoN(G_crt,OD,dummy_nodes)
        if step == 'line_search':
            a_k,obj_k=line_search(G_crt,y_k,edge_list)#include the fixed step size,
        elif step == 'fixed':
            a_k,obj_k=fixed_step(G_crt,y_k,edge_list,i)
        else:
            print("wrong optim step chosen")
            return

        #update the flows
        G_crt=update_flows(G_crt,y_k,a_k,edge_list)

        #save for analyses
        opt_res['obj'].append(obj_k)
        opt_res['a_k'].append(a_k)
        G_k=G_crt
        G_list.append(G_k)
        y_list.append(y_k)
        OD_list.append(OD.copy())
        i+=1
    return G_list,y_list,opt_res,OD_list