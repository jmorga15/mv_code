# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 07:52:14 2019

@author: Jonathan
"""
import numpy as np
from functions_mv_v2 import *
import random

def tau_calc(t, t_max, count, mv_lifetime, coordinates, ka, kr, k_01, k_10, k_12, k_23, k_34, k_45, k_54, k_53, k_40, k_20, R, r, ave_full):
    t_signal = -1
    
    ### tau calculation for nth level activation ###
    ### number of MV capable of nth level activation ###
    
    ### k_on: need to figure out how to account for multiple antigen ###
#    k_01_MV = np.sum((coordinates[:,2] == 1))
#    if k_01_MV != 0:
#        tau_k_01 = gen_tau(k_01, k_01_MV)
#    else:
#        tau_k_01 = 1000
    
    ### an array with length of max antigen covered by MV ###
    if np.sum(coordinates[:,2] == 1) != 0:
#    if np.max(coordinates[:,5]) != 0:
        tau_k_01_array = np.zeros((1,int(np.max(coordinates[:,5]))))
        for i in range(int(np.max(coordinates[:,5]))):
            loc_nodes = np.where(coordinates[:,2] == 1)[0]
            k_01_MV = np.sum((coordinates[loc_nodes,5] == i + 1))
            if k_01_MV != 0:
                tau_k_01 = gen_tau(k_01 * (i+1), k_01_MV)
            else:
                tau_k_01 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))
            tau_k_01_array[0,i] = tau_k_01
    else:
        tau_k_01_array = np.zeros((1,5)) + np.random.uniform(low=0.0, high=1.0, size=(5,)) + 100000.0
#        tau_k_01 = np.min(tau_k_01_array[0])
#    else:
#        tau_k_01 = 1000
    #print(tau_k_01)

    k_10_MV = np.sum((coordinates[:,2] == 2))
    if k_10_MV != 0:
        tau_k_10 = gen_tau(k_10, k_10_MV)
    else:
        tau_k_10 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    k_12_MV = np.sum((coordinates[:,2] == 2))
    if k_12_MV != 0:
        tau_k_12 = gen_tau(k_12, k_12_MV)
    else:
        tau_k_12 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

#    k_21_MV = np.sum((coordinates[:,2] == 3))
#    if k_21_MV != 0:
#        tau_k_21 = gen_tau(k_21, k_21_MV)
#    else:
#        tau_k_21 = 1000

    k_23_MV = np.sum((coordinates[:,2] == 3))
    if k_23_MV != 0:
        tau_k_23 = gen_tau(k_23, k_23_MV)
    else:
        tau_k_23 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    k_34_MV = np.sum((coordinates[:,2] == 4))
    if k_34_MV != 0:
        tau_k_34 = gen_tau(k_34, k_34_MV)
        #print("k_34 =" , tau_k_34)
#        if tau_k_34 > t_max - t:
#            #print("hi")
#            t_signal = t_max - t 
    else:
        tau_k_34 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

#    k_43_MV = np.sum((coordinates[:,2] == 5))
#    if k_43_MV != 0:
#        tau_k_43 = gen_tau(k_43, k_43_MV)
#    else:
#        tau_k_43 = 1000
        
    ### an array with length of max antigen covered by MV ###
    if np.sum(coordinates[:,2] == 5) != 0:
#    if np.max(coordinates[:,5]) != 0:
        tau_k_45_array = np.zeros((1,int(np.max(coordinates[:,5]))))
        for i in range(int(np.max(coordinates[:,5]))):
            loc_nodes = np.where(coordinates[:,2] == 5)[0]
            k_45_MV = np.sum((coordinates[loc_nodes,5] == i + 1))
            if k_45_MV != 0:
                tau_k_45 = gen_tau(k_45 * (i+1), k_45_MV)
            else:
                tau_k_45 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))
            tau_k_45_array[0,i] = tau_k_45
    else:
        tau_k_45_array = np.zeros((1,5)) + np.random.uniform(low=0.0, high=1.0, size=(5,)) + 100000.0

#    k_45_MV = np.sum((coordinates[:,2] == 5))
#    if k_45_MV != 0:
#        tau_k_45 = gen_tau(k_45, k_45_MV)
#    else:
#        tau_k_45 = 1000

    k_54_MV = np.sum((coordinates[:,2] == 6))
    if k_54_MV != 0:
        tau_k_54 = gen_tau(k_54, k_54_MV)
    else:
        tau_k_54 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    k_53_MV = np.sum((coordinates[:,2] == 6))
    if k_53_MV != 0:
        tau_k_53 = gen_tau(k_53, k_53_MV)
    else:
        tau_k_53 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    k_40_MV = np.sum((coordinates[:,2] == 5))
    if k_40_MV != 0:
        tau_k_40 = gen_tau(k_40, k_40_MV)
    else:
        tau_k_40 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    k_20_MV = np.sum((coordinates[:,2] == 3))
    if k_20_MV != 0:
        tau_k_20 = gen_tau(k_20, k_20_MV)
        #t_signal = tau_k_20
#        if tau_k_20 > t_max - t:
#            #print("hi")
#            t_signal = t_max - t
    else:
        tau_k_20 = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    ### count number of active MV that are in state 0 and can be removed ###
    on = np.sum(coordinates[:,2] == 0) + np.sum(coordinates[:,2] == 1)



    off = ave_full - len(coordinates)

    if off != 0:
        tau_A = gen_tau(ka, off)
    else:
        tau_A = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    ### tau calculations for removal of 0 state nodes ###
    if on != 0:
        tau_R = gen_tau(kr, on)
    else:
        tau_R = 100000 + np.random.uniform(low=0.0, high=1.0, size=(1,))

    #tau_vec = np.array([tau_A, tau_R, tau_k_01, tau_k_12, tau_k_10, tau_k_23, tau_k_34, tau_k_45, tau_k_54, tau_k_53, tau_k_40, tau_k_20 ])
    tau_vec = np.array([tau_A, tau_R, tau_k_12, tau_k_10, tau_k_23, tau_k_34, tau_k_54, tau_k_53, tau_k_40, tau_k_20 ])
    
#    return tau_vec



#
#    tau_vec = np.array([tau_A, tau_R, tau_k_01, tau_k_12, tau_k_10, tau_k_21, tau_k_23, tau_k_34, tau_k_43, tau_k_45, tau_k_54, tau_k_53, tau_k_40, tau_k_20 ])
#
        ### append the tau_A to the tau_vec ###
#
#
    ### find the location of the minimum tau in tau_vec ###
    min_location1 = np.where(tau_vec == np.min(tau_vec))
    
    ### find the minimum location in the k_on vec ###
    min_location2 = np.where(tau_k_01_array[0] == np.min(tau_k_01_array[0]))
    
    ### find the minimum location in the second instance k_on vec ###
    min_location3 = np.where(tau_k_45_array[0] == np.min(tau_k_45_array[0]))
    
#    
#    print(min_location2[0])
#    print(min_location1[0])
#    print(min_location3[0])
#    print(tau_k_01_array[0, min_location2[0]])
#    print(tau_k_45_array[0, min_location3[0]])
#    print(tau_vec[min_location1[0]])
#    print("coordinates", coordinates)
#    print(tau_vec)
    
    if (tau_vec[min_location1[0]] < tau_k_01_array[0, min_location2[0]]) and (tau_vec[min_location1[0]] < tau_k_45_array[0, min_location3[0]]):
        tau = tau_vec[min_location1[0]]
        min_location = min_location1

        #(x_new,y_new) = gen_coords(coordinates, R, r)
        ### check if the minimum tau is the activation position ###
        if min_location[0] == 0:
            (x_new,y_new) = gen_coords(coordinates, R, r)
    
            ### A node got turned on ###
            loc_node = (x_new, y_new)
            coordinates = np.append(coordinates,[[x_new,y_new, 0,-1, 0, 0]], axis = 0)
            tau = tau_A
            #print(on_off)
        elif min_location[0] == 1:
            ### a node is removed ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(on,1))
            loc_nodes = np.where(coordinates[:,2] < 2)[0]
            loc_node = np.where(np.max(choose_node) == (choose_node))[0]
#            print("loc node =", loc_node)
#            print(choose_node)
            #print(coordinates[loc_nodes[loc_node],4])
            ### get the lifetime of dead mv ###
            mv_lifetime[0, count] = coordinates[loc_nodes[loc_node], 4] + tau_R
            count = count + 1
            #print(np.max(mv_lifetime))
            
            coordinates = np.delete(coordinates, loc_nodes[loc_node],0)
            tau = tau_R
            
#        elif min_location[0] == 2:
#            ### MV is activated to level 2 ###
#            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_01_MV,1))
#            loc_nodes = np.where(coordinates[:,2] == 1)[0]
#            #print(loc_nodes)
#            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
#            #print(loc_node)
#            coordinates[loc_nodes[loc_node],2] = 2
#            #print(coordinates)
#            
            
        elif min_location[0] == 2:
            ### MV is activated to level 3 ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_12_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 2)[0]
#            print(loc_nodes)
#            print(choose_node)
#            print(tau_k_12)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 3
            #print(coordinates)
        elif min_location[0] == 3:
            ### MV is reduced to level 1 from 2 ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_10_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 2)[0]
#            print(loc_nodes)
#            print(choose_node)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 1
            #print(coordinates)
    #    elif min_location[0] == 5:
    #        ### MV is reduced to level 2 from 3 ###
    #        choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_21_MV,1))
    #        loc_nodes = np.where(coordinates[:,2] == 3)[0]
    #        #print(loc_nodes)
    #        loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
    #        #print(loc_node)
    #        coordinates[loc_nodes[loc_node],2] = 2
    #        #print(coordinates)
        elif min_location[0] == 4:
            #print("hi")
            ### T_cell is activated - "DESTROY" ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_23_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 3)[0]
#            print(loc_nodes)
#            print(choose_node)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 4
            #print(coordinates)
        elif min_location[0] == 5:
            ### Entering unbound secure state ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_34_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 4)[0]
#            print(loc_nodes)
#            print(choose_node)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 5
            #print(coordinates)
    #    elif min_location[0] == 7:
    #        ### Entering unbound secure state ###
    #        choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_43_MV,1))
    #        loc_nodes = np.where(coordinates[:,2] == 5)[0]
    #        #print(loc_nodes)
    #        loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
    #        #print(loc_node)
    #        coordinates[loc_nodes[loc_node],2] = 4
    #        #print(coordinates)
#        elif min_location[0] == 6:
#            ### Entering unbound secure state ###
#            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_45_MV,1))
#            loc_nodes = np.where(coordinates[:,2] == 5)[0]
#            #print(loc_nodes)
#            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
#            #print(loc_node)
#            coordinates[loc_nodes[loc_node],2] = 6
#            #print(coordinates)
        elif min_location[0] == 6:
            ### Entering unbound secure state ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_54_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 6)[0]
#            print(loc_nodes)
#            print(choose_node)
#            print(tau_k_54)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 5
            #print(coordinates)
        elif min_location[0] == 7:
            ### Entering unbound secure state ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_54_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 6)[0]
#            print(loc_nodes)
#            print(choose_node)
#            print(tau_k_54)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 4
            #print(coordinates)
        elif min_location[0] == 8:
            ### Entering unbound secure state ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_40_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 5)[0]
#            print(loc_nodes)
#            print(choose_node)
#            print(tau_k_40)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 1
            #print(coordinates)
        elif min_location[0] == 9:
            ### Entering unbound secure state ###
            choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_20_MV,1))
            loc_nodes = np.where(coordinates[:,2] == 3)[0]
#            print(loc_nodes)
#            print(tau_k_20)
            loc_node = np.where(np.max(coordinates[loc_nodes,2] * choose_node) == (coordinates[loc_nodes,2] * choose_node))[0]
            #print(loc_node)
            coordinates[loc_nodes[loc_node],2] = 1
            #print(coordinates)
            #t_signal = tau_k_20
    elif tau_k_01_array[0, min_location2[0]] < tau_vec[min_location1[0]] and  tau_k_01_array[0, min_location2[0]] < tau_k_45_array[0, min_location3[0]]:
  #    if tau_vec[min_location1[0]] < tau_k_01_array[0, min_location2[0]] and tau_vec[min_location1[0]] < tau_k_45_array[0, min_location3[0]]:   ### kon occurred in first instance ###
        min_location = min_location2

        ### find the tau that corresponds to the kon event ###
        tau = tau_k_01_array[0, min_location2[0]]

        
        ### location of nodes that can level up ###
        loc_nodes1 = np.where(coordinates[:,2] == 1)[0]
        
        ### number of nodes capable of being upgraded that match the min_location2 ###
        ### min_location2 should also be the number of antigen under MV ###
        k_on1_MV = np.sum((coordinates[loc_nodes1,5] == min_location2[0] + 1 ))
        
        ### random vector to decide which mv is promoted ###
        choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_on1_MV,))
        
        ### location of nodes in level 1 and have the right number of antigen covered ###
        ### need to add on to the location because indexing starts at 0, but the indexing plus 1 is ###
        ### the number of antigen present covered by the microvilli ###
        
        if len(coordinates[loc_nodes1,5] == min_location2[0] + 1) == 0:
            
            print("tau_k01_array =",tau_k_01_array[0] )
            print("min1 =",min_location1[0] )
            print("min2 =",min_location2[0] )
            print("min3 =",min_location3[0] )
            print("loc_nodes1  =",loc_nodes1 )
            print("coordinates  =",coordinates )

            
        
        loc_nodes = np.where((coordinates[loc_nodes1,5] == min_location2[0] + 1))[0]
        
        ### now we determine which MV is promoted by taking a max of random vector broadcasted ###
        ### onto the coordinate array constants###
        loc_node = np.where(np.max(coordinates[loc_nodes1[loc_nodes],2] * choose_node) == (coordinates[loc_nodes1[loc_nodes],2] * choose_node))[0]
        
        ### absurd indexing becuase i suck ###
        coordinates[loc_nodes1[loc_nodes[loc_node]],2] = 2
        
    else :
        ### kon occurred in second instance ###
        min_location = min_location3

        ### find the tau that corresponds to the kon event ###
        tau = tau_k_45_array[0, min_location3[0]]

        
        ### location of nodes that can level up ###
        loc_nodes1 = np.where(coordinates[:,2] == 5)[0]
        
        ### number of nodes capable of being upgraded that match the min_location2 ###
        ### min_location2 should also be the number of antigen under MV ###
        k_on2_MV = np.sum((coordinates[loc_nodes1,5] == min_location3[0] + 1 ))
        
        ### random vector to decide which mv is promoted ###
        choose_node = np.random.uniform(low=0.0, high=1.0, size=(k_on2_MV,))
        
        ### location of nodes in level 1 and have the right number of antigen covered ###
        ### need to add on to the location because indexing starts at 0, but the indexing plus 1 is ###
        ### the number of antigen present covered by the microvilli ###
        loc_nodes = np.where((coordinates[loc_nodes1,5] == min_location3[0] + 1))[0]
        
        ### now we determine which MV is promoted by taking a max of random vector broadcasted ###
        ### onto the coordinate array constants###
        loc_node = np.where(np.max(coordinates[loc_nodes1[loc_nodes],2] * choose_node) == (coordinates[loc_nodes1[loc_nodes],2] * choose_node))[0]
        
        ### absurd indexing becuase i suck ###
        coordinates[loc_nodes1[loc_nodes[loc_node]],2] = 6
        
        
        
        
    return (coordinates, tau, count, mv_lifetime, t_signal)
