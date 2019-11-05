#this document contains a few functions that were tried before in the FW scheme


def init_flows(G,OD):
    for e in G.edges:
        if e[1]=='R':
            G[e[0]][e[1]]['f_m']=0
            G[e[0]][e[1]]['f_r']=0 #G[e[0]][e[1]]['k'] #The initialization cannot work if the dummy edges are full in the beginning? 
        else:
            G[e[0]][e[1]]['f_m']=0
            G[e[0]][e[1]]['f_r']=0 
    return G

def init_flows(G,OD,dummy_nodes,opt=1):
    #Initiliaze the flows, with a feasible solution
    #still unclear whether we need to initialize the rebalancers appropriately too
    if opt==1:
        for (o,d) in OD.keys():
            #assign to dummy edge
            G[o][d]['f_m']+=OD[o,d]
            new_o=dummy_nodes[d]
            path=nx.shortest_path(G,source=new_o,target='R',weight='cost')
            for i in range(len(path)-1):
                G[path[i]][path[i+1]]['f_r']+=OD[o,d] 
    elif opt==2:
        for (o,d) in OD.keys():
            path=nx.shortest_path(G,source=o,target=d,weight='cost')
            for i in range(len(path)-1):
                G[path[i]][path[i+1]]['f_m']+=OD[o,d]
    else:
        print("Not correct option")
        return 
    
    return G


def init_flows2(G,OD):#here no initialization of the R edges, as they 
    #Initiliaze the flows, with a feasible solution
    #still unclear whether we need to initialize the rebalancers appropriately too
    for (o,d) in OD.keys():
        if d!='R':
            path=nx.shortest_path(G,source=o,target=d,weight='cost')
            for i in range(len(path)-1):
                if d=='R':
                    G[path[i]][path[i+1]]['f_r']+=OD[o,d]
                else:
                    G[path[i]][path[i+1]]['f_m']+=OD[o,d]
    return G



#computes the total cost under the current model
# def Total_Cost(G,y_k,a_k,edge_list):
#     F_E=0
#     for i in range(len(edge_list)):#you know for sure exactly what edge it is for
#         e=edge_list[i]
#         x_k_e=G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_r'] # retrieve the flow, total
#         y_k_e=y_k[(e[0],e[1]),'f_m']+y_k[(e[0],e[1]),'f_r'] #retrieve the flow from the manual assignment, total
        
#         flow_tmp=x_k_e+a_k*(y_k_e-x_k_e)
        
#         #retrieve parameters to compute the BPR
#         phi=G[e[0]][e[1]]['phi']
#         k=G[e[0]][e[1]]['k']
#         # sign=G[e[0]][e[1]]['sign']
#         if k <10**-5:#you eliminate the edges that are considered non-usable
#             continue

#         if e[1]=='R':#not including the cost of edges 1R and 2R might make sense, as we want to rebalance whatever happens
#             continue
#         # print(k)
#         F_E+=BPR_int(phi,flow_tmp,k)#I am assuming there will be syntaxic problems there
        
        
#         #this has to be included because it is directly included in the definition of the cost function
#         if G[e[0]][e[1]]['sign']==(-1): #we have a negative edge
#             F_E-=flow_tmp*80#INVERSE_DEMAND_SHIFT
        
#         # not entirely sure this needs to be here
#         # if 'pot' in G.nodes[e[1]]:
#         #     F_E+=G.nodes[e[1]]['pot']*flow_tmp
            
#     return F_E