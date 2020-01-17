import numpy as np 
import cvxpy as cp
import networkx as nx
from routines_icu import *
from helpers_icu import *


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


#we currently assign a rebalancing flow to y_k, but this is useless currently
#as with the new cost function we do not take it into account
def AoN(G,OD,dummy_nodes):
    #perform the All Or Nothing assignment    
    y_k=init_y(G)
    eps=10**-6

    alpha=1.5 #TODO: what is this alpha for? 
    for (o,d) in OD.keys():
        U=OD[o,d]
        if U>eps:
            if d=='R':
                flag='f_r'
            else:
                flag='f_m'
            # print("AON, (o,d):", o,d)
            path=nx.shortest_path(G,source=o,target=d,weight='cost')

            #in this case, we are assigning to the dummy edge
            #again, not sure that the dummy edges are truly necessary in the end
            #TODO: figure out what the commented version here below actually meant
            #it seems to me like it was a fix to make sure we introduce the 
            #"complement" of the demand
            #we have to see whether or not it was helping in the ICU case... ? 
            if len(path)==2 and d!='R':
                # y_k[(o,d),'f_m']+=alpha*(U-G[dummy_nodes[d]][d]['f_m'])
                y_k[(o,d),'f_m']+=U
            else:
                for i in range(len(path)-1):
                    y_k[(path[i],path[i+1]),flag]+=U
    return y_k

def modified_FW(G_0,OD,edge_list,dummy_nodes,maxIter=50, step='line_search',
       rebalancer_smoothing=False, ri_smoothing=True,evolving_bounds=True):

    #might have to "reorganize" the order of operations, in order to make it more understandable/logical


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

    while i<maxIter: #introduce stopping criterion 
        G_crt=G_k.copy()

        #deal with rebalancers
        #estimate #ri_k, update OD, assign, update costs
        ri_k,G_crt=estimate_ri_k(G_crt,dummy_nodes, ri_smoothing, a_k)
        OD=update_OD(OD,ri_k,a_k,dummy_nodes, G_crt, evolving_bounds)
        G_crt=update_capacities(G_crt,ri_k, dummy_nodes)

        G_crt=assign_rebalancers(G_crt,OD,rebalancer_smoothing,a_k)
        G_crt=update_costs(G_crt,80)
        #disp_costs(G_crt)#debug helper
        
        #perform AON assignment
        y_k=AoN(G_crt,OD,dummy_nodes)
        if step == 'line_search':
            a_k,obj_k=line_search(G_crt,y_k,edge_list)
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

def assign_rebalancers(G,OD, rebalancer_smoothing,a_k):
    #assign rebalancers to shortest path, in a definite manner 
    #i.e. we do not optimize on that, and consider them fixed at every iteration

    if not rebalancer_smoothing:
        beta=1 #we only take the new assignment into account
    else:
        beta=a_k #we take the step size that is given to us, and update the rebalancers accordingly

    for (o,d) in OD.keys():
        N_r=OD[o,d]
        if d=='R' and N_r >10**-6:
            path=nx.shortest_path(G,source=o,target=d,weight='cost')
            for i in range(len(path)-1):
                G[path[i]][path[i+1]]['f_r']=(1-beta)*G[path[i]][path[i+1]]['f_r']+beta*N_r #introduce some kind of smoothing of the rebalancing update
    return G

def estimate_ri_k(G,dummy_nodes,ri_smoothing,a_k):
    #determine whether each node is in excess or deficit of rebalancers
    ri_k=dict()
    ri_k_prev=dict()
    
    if ri_smoothing:
        beta=a_k
    else:
        beta=1 #no smoothing

    for n in G.nodes():
        ri_k[n]=0
        ri_k_prev[n]=G.nodes[n]["ri"]
    for e in G.edges():
        if e[1] not in dummy_nodes.keys():
            ri_k[e[0]]+=G[e[0]][e[1]]['f_m']
            ri_k[e[1]]-=G[e[0]][e[1]]['f_m']
    for n in G.nodes():
        G.nodes[n]["ri"]= (1-beta) * ri_k_prev[n] + beta*ri_k[n]
    return ri_k,G

def update_OD(OD,ri_k, a_k, dummy_nodes, G, evolving_bounds=True):
    
    #update the OD pairs for rebalancers
    eps=10**-6
    for n in ri_k.keys():
        if n!='R' and n not in dummy_nodes.keys():
            if ri_k[n]<-eps: #you are in excess
                OD[(n,'R')]=-ri_k[n]
            else:
                OD[(n,'R')]=0

    #update the bounds for "regular OD pairs"
    #todo: make sure that the costs are a good reflection of the flow at iteration k 
    #tocode: if x_k on inverse demand edge is too close to upper bound or lower bound you 
    #have to decrease or increase those bounds, respectively

    #currently, we keep lower bounds at zero
    if evolving_bounds:
        #parameters, to be given as arguments in the future
        #those parameters should really be given as arguments
        l1=10**-1
        l2=2 #such a high value currently disables it
        alpha=1.2 #I think a moving value on that would be better. And same for the l1/l2. should be based o2
        for (o,d) in OD.keys():
            
            if d in dummy_nodes.keys() and d!='R': #only the nodes that are the dummy nodes
                crt_Ulim=OD[o,d]#currently, we treat only the upper limit 
                crt_p_flow=G[dummy_nodes[d]][d]['f_m']
                rel_U_error=abs(crt_p_flow-crt_Ulim)/crt_p_flow
                # UPPER BOUND
                if rel_U_error <= l1:#x_k too close to U
                    new_Ulim=alpha*crt_p_flow
                elif rel_U_error>=l2:#x_k too far from U
                    new_Ulim=alpha*crt_p_flow
                else:
                    new_Ulim=crt_Ulim
                OD[o,d]=new_Ulim 

    return OD

def update_capacities(G,ri_k, dummy_nodes): 
    eps=10**-6
    for n in G.nodes():
        if n not in dummy_nodes.keys() and n!='R':
            if ri_k[n]>eps: 
                G[n]['R']['k']=ri_k[n]
            else:
                G[n]['R']['k']=eps
    return G

def init_y(G):
    y_k=dict()
    for e in G.edges():
        y_k[e,'f_r']=0
        y_k[e,'f_m']=0
    return y_k



#computes the total cost under the current model
#modified version that takes into account the fact that we compute the total cost
#on passengers only!!
def Total_Cost(G,y_k,a_k,edge_list):
    F_E=0
    for i in range(len(edge_list)):#you know for sure exactly what edge it is for
        e=edge_list[i]
        x_k_e_m=G[e[0]][e[1]]['f_m']
        x_k_e_r = G[e[0]][e[1]]['f_r'] 
        y_k_e_m=y_k[(e[0],e[1]),'f_m'] #we do not retrieve the y_k_r flow, as we do not use it here
        
        flow_tmp=x_k_e_m+a_k*(y_k_e_m-x_k_e_m)#this is only a flow of rebalancers
        
        #retrieve parameters to compute the BPR
        phi=G[e[0]][e[1]]['phi']
        k=G[e[0]][e[1]]['k']

        if k <10**-5:#you eliminate the edges that are considered non-usable
            continue
        if e[1]=='R':#not including the cost of edges 1R and 2R might make sense, as we want to rebalance whatever happens
            continue
        

        ###### IMPORTANT NOTE ####
        #is that really correct... ? 
        #if you consider you have some "exogenous flow" due to rebalancing, 
        #there should be a -BPR_int(...,x_k_e_r,...)
        #Basically because you are integrating over the consumer flow and because
        #the cost function has changed... 
        #
        F_E+=BPR_int(phi,flow_tmp + x_k_e_r,k)#I am assuming there will be syntaxic problems there
        
        #this has to be included because it is directly included in the definition of the cost function
        if G[e[0]][e[1]]['sign']==(-1): #we have a negative edge
            F_E-=(flow_tmp+x_k_e_r)*80#INVERSE_DEMAND_SHIFT
        
        # not entirely sure this needs to be here
        # if 'pot' in G.nodes[e[1]]:
        #     F_E+=G.nodes[e[1]]['pot']*flow_tmp
            
    return F_E


def fixed_step(G,y_k,edge_list,k):
    #here we only compute for 'f_m', because we do not want to update 
    # the rebalancers and the passengers at the same time

    gamma=2/(k+2)

    for i in range(len(edge_list)):
        e=edge_list[i]
        for flag in ['f_m']:
            x_k_e=G[e[0]][e[1]][flag] # retrieve the flow
            y_k_e=y_k[(e[0],e[1]),flag] #retrieve the flow from the manual assignment
            G[e[0]][e[1]][flag]=(1-gamma)*x_k_e+gamma*y_k_e
        
    return gamma, Value_Total_Cost(G)
    


def line_search(G,y_k,edge_list):
    a_k=cp.Variable()
    constraints=[a_k>=0, a_k<=1]
    obj=Total_Cost(G,y_k,a_k,edge_list)
    prob=cp.Problem(cp.Minimize(obj),constraints)
    prob.solve(verbose=False)
    print(prob.status)
    try:
        obj_val=obj.value
    except:
        obj_val=0
    return a_k.value,obj_val

def update_flows(G,y_k,a_k,edge_list):
    
    for i in range(len(edge_list)):
        e=edge_list[i]
        for flag in ['f_m']:#, 'f_r']: we consider the flow of rebalancers fixed during the computation of costs
            x_k_e=G[e[0]][e[1]][flag] # retrieve the flow
            y_k_e=y_k[(e[0],e[1]),flag] #retrieve the flow from the manual assignment
            G[e[0]][e[1]][flag]=(1-a_k)*x_k_e + a_k * y_k_e
    return G   






