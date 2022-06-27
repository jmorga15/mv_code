# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 15:19:39 2019

@author: Jonathan
"""
import random
import numpy as np

def gen_coords(coordinates, R, r):
    if len(coordinates) == 0:
        #print("hello")
        ### Try a new position for the next MV ###
        x_new = random.uniform(-(R - r), (R - r))
        y_new = random.uniform(-(R - r), (R - r))

        ### determine the new coordinates distance from the origin ###
        d_org =  np.sqrt( (x_new - 0)**2 + (y_new - 0)**2)
        d_org = d_org <= (R - r)

        counter = 0
        ### keep generating new coordinates until a new MV fits or made a certain number of attempts ###
        skip = 0
        while d_org == False and counter < 500000 and skip == 0:
            x_new = random.uniform(-(1 - r), (1 - r))
            y_new = random.uniform(-(1 - r), (1 - r))

            d_org =  np.sqrt( (x_new - 0)**2 + (y_new - 0)**2)
            d_org = d_org <= (1 - r)

            counter = counter + 1
            if counter >= 500000:
                #print("hi")
                skip = 1
    else:
        ### Try a new position for the next MV ###
        #x_new = random.uniform(-(R - r), (R - r))
        #y_new = random.uniform(-(R - r), (R - r))
        
        theta_new = random.uniform(0, 2*np.pi)
        r_new = np.sqrt(random.uniform(0,(R - r)))
        
        x_new = r_new * np.cos(theta_new)
        y_new = r_new * np.sin(theta_new)
        ### determine the new coordinates distance from the origin ###
        d_org =  np.sqrt( (x_new - 0)**2 + (y_new - 0)**2)
        d_org = d_org <= (R - r)

        ### determine the new coordinates distance from all other MV ###
        d = np.sqrt( (x_new - coordinates[:,0])**2 + (y_new - coordinates[:,1])**2)
        d = d >= 2*r

        counter = 0
        ### keep generating new coordinates until a new MV fits or made a certain number of attempts ###
        skip = 0
        while ((False in d) == True or d_org == False) and counter < 200000 and skip == 0:
            x_new = random.uniform(-(1 - r), (1 - r))
            y_new = random.uniform(-(1 - r), (1 - r))
            d = np.sqrt( (x_new - coordinates[:,0])**2 + (y_new - coordinates[:,1])**2)
            d = d >= 2*r

            d_org =  np.sqrt( (x_new - 0)**2 + (y_new - 0)**2)
            d_org = d_org <= (1 - r)

            counter = counter + 1
            if counter >= 200000:
                skip = 1
                print("hi")

    return(x_new,y_new)

def gen_tau(kx, l):
    ### generate a random vector of size l ###
    u = np.random.uniform(low=0.0, high=1.0, size=(l,))

    ###calculate the tau vec ###
    tau = 1/kx*np.log(1.0/u)

    return(np.min(tau))
    
def gen_tau_kon(kx, l, coordinates):
    
    loc_nodes = np.where(coordinates[:,2] == 1)
    ### generate a random vector of size l ###
    u = np.random.uniform(low=0.0, high=1.0, size=(l,))

    ###calculate the tau vec ###
    tau = 1/kx*np.log(1.0/u)

    return(np.min(tau))

#def MV_antigen_dist(theta, r_antigen, coordinates, r):
#    x_antigen = r_antigen*np.cos(theta)
#    y_antigen = r_antigen*np.sin(theta)
#
#    ### determine MV distance from antigen ###
#    d = np.sqrt( (x_antigen - coordinates[:,0])**2 + (y_antigen - coordinates[:,1])**2)
#    MV_activation = (d <= r)
#    hit = np.where(MV_activation == True)[0]
#    if coordinates[hit,2] < 2:
#        coordinates[hit,2] = 1
#
#    return(coordinates)
##    if True in MV_activation:
##        print(coordinates)

def antigen_antigen_dist(antigen_coordinates, r):
    r_antigen = np.sqrt(np.random.uniform(low=0.0, high=1, size=(1,)))
    theta = np.random.uniform(low=0.0, high=1.0, size=(1,))*np.pi*2

    x_antigen = r_antigen*np.cos(theta)
    y_antigen = r_antigen*np.sin(theta)

    ## determine distance from other antigens, needed if preventing overlap ###
    d = np.sqrt( (x_antigen - antigen_coordinates[:,0]*np.cos(antigen_coordinates[:,1]))**2 + (y_antigen - antigen_coordinates[:,0]*np.sin(antigen_coordinates[:,1]))**2)
    d = d >= 0.01*r
#
    counter = 0
    ### keep generating new coordinates until a new antigen fits or made a certain number of attempts, needed if preventing overlap ###
    skip = 0
    while ((False in d) == True) and counter < 50000 and skip == 0:
        r_antigen = np.random.uniform(low=0.0, high=1.0, size=(1,))
        theta = np.random.uniform(low=0.0, high=1.0, size=(1,))*np.pi*2
        x_antigen = r_antigen*np.cos(theta)
        y_antigen = r_antigen*np.sin(theta)

        ### determine distance from other antigens ###
        d = np.sqrt( (x_antigen - antigen_coordinates[:,0]*np.cos(antigen_coordinates[:,1]))**2 + (y_antigen - antigen_coordinates[:,0]*np.sin(antigen_coordinates[:,1]))**2)
        d = d >= 0.01*r

    return(r_antigen, theta)




def MV_antigen_dist2(antigen_coordinates, coordinates, r):
    ## find where MV cover antigen ##
    
    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])
    
    if len(coordinates) != 0:
        
        for i in range(len(coordinates)):
            
            d = np.sqrt( (x_antigen - coordinates[i,0])**2 + (y_antigen - coordinates[i,1])**2)
            #print(d)
            MV_activation = (d <= r)
            
            #print (MV_activation)
            ## Find antigen covered by MV ##
            if True in MV_activation:
                antigen_coordinates[np.where(MV_activation==True),2] = 2 
            
#                if coordinates[i,2] < 2 and coordinates[i,5] == 0:
#                    coordinates[i,2] = 1
#                    coordinates[i,5] = 1
    #print (np.max(coordinates[:,5]))
        
#    for i in antigen_coordinates:
#
#        r_antigen = i[0]
#        theta = i[1]
#
#        x_antigen = r_antigen*np.cos(theta)
#        y_antigen = r_antigen*np.sin(theta)
#
#        ### determine MV distance from antigen ###
#        d = np.sqrt( (x_antigen - coordinates[:,0])**2 + (y_antigen - coordinates[:,1])**2)
#        MV_activation = (d <= r)
#        hit = np.where(MV_activation == True)[0]
#        if len(coordinates) != 0:
#            #print(coordinates)
#            #print(hit)
#            
#            if len(hit) > 1:
#                ### something went wrong and there is mv overlap ###
#                print(coordinates[hit])
#                
##            elif coordinates[hit,2] < 2 and coordinates[hit,5] > 0:
##                ### there are more than one antigen under mv ###
##                coordinates[hit,5] = coordinates[hit,5] + 1
#                
#            elif coordinates[hit,2] < 2:
#                coordinates[hit,2] = 1

    return(antigen_coordinates)
    #    if True in MV_activation:
    #        print(coordinates)
def MV_antigen_dist3(antigen_coordinates, coordinates, r, loc, P, N):
    ## working on removing MV and dissociating corresponding antigen ##
    
    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])
    #print("loc",loc)
    d = np.sqrt( (x_antigen - coordinates[loc,0])**2 + (y_antigen - coordinates[loc,1])**2)
    #print(d)
    MV_activation = (d <= r)
            
    #print (MV_activation)
    if True in MV_activation:
        #print("P",P)
        #print("zeroed antigen full",antigen_coordinates[np.where(MV_activation==True),2])
        #print("zeroed antigen",antigen_coordinates[np.where(MV_activation==True),2][0])
        ant_zero = antigen_coordinates[np.where(MV_activation==True),2][0]
        #print(int(ant_zero))
        for i in ant_zero:
            P[int(i),0] = P[int(i),0] - 1
        antigen_coordinates[np.where(MV_activation==True),2] = 0 
        #print("P",P)
    #coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])
    
    return (antigen_coordinates,P)


def MV_antigen_dist1(antigen_coordinates, coordinates, r,P):
    
    ### this function finds the number of antigen that are within an MV contact Boundary ###
    
    ### takes the state matrix for antigen and MV (antigen_coordinates and coordinates), as well as MV radius and population vector P ###
    
    ### modifies the antigen states in antigen_coordinates to reflect coverage by MV (state 2) ###
    
    
    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])
    
    d = np.sqrt( (x_antigen - coordinates[-1,0])**2 + (y_antigen - coordinates[-1,1])**2)
    #print(d)
    MV_activation = (d <= r)
            
    #print (MV_activation)
    if True in MV_activation:
        #print("P",P)
        #print("zeroed antigen full",antigen_coordinates[np.where(MV_activation==True),2])
        #print("zeroed antigen",antigen_coordinates[np.where(MV_activation==True),2][0])
        ant_zero = antigen_coordinates[np.where(MV_activation==True),2][0]
        #print(int(ant_zero))
        for i in ant_zero:
            P[2,0] = P[2,0] + 1
        antigen_coordinates[np.where(MV_activation==True),2] = 2 
        #print("P",P)
    #coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])
    
    return antigen_coordinates,P


def MV_antigen_dist4(antigen_coordinates, coordinates, r, loc, P, N):
    ##antigen transition that allows mv to become stabilized ##
    # find which MV is affected
    x_antigen = antigen_coordinates[loc,0] * np.cos(antigen_coordinates[loc,1])
    y_antigen = antigen_coordinates[loc,0] * np.sin(antigen_coordinates[loc,1])
    #print("loc",loc)
    d = np.sqrt( (x_antigen - coordinates[:,0])**2 + (y_antigen - coordinates[:,1])**2)
    #print(d)
    MV_activation = (d <= r)
            
    #print (MV_activation)
    if True in MV_activation:
        
        #print("P",P)
        #print("zeroed antigen full",antigen_coordinates[np.where(MV_activation==True),2])
        #print("zeroed antigen",antigen_coordinates[np.where(MV_activation==True),2][0])
        # if this is the first activated antigen and TCR under the MV
        # need to update the MV state such that it can now be stabilized
        if coordinates[np.where(MV_activation==True),2] < N+4:
            coordinates[np.where(MV_activation==True),2]=N+4
        #print(int(ant_zero))
            #update MV states
            P[N+4,0] = P[N+4,0] + 1
            P[1,0] = P[1,0] - 1
            
        # if the mv was stabilized already and no tcr was signalling, this transition could prevent destabilization ku #
        elif coordinates[np.where(MV_activation==True),2] == N+6:
            coordinates[np.where(MV_activation==True),2]= N+5
            #print("hi")
            P[N+6,0] = P[N+6,0] - 1
            P[N+5,0] = P[N+5,0] + 1
            #coordinates[loc2,2]=15
            # update antigen states
#            P[transit[0][0],0] = P[transit[0][0],0] - 1
#            P[transit[1][0],0] = P[transit[1][0],0] + 1

        
        #print("P",P)
    #coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])
    
    return (coordinates,P)

def MV_antigen_dist5(antigen_coordinates, coordinates, r, loc, P, N):
    ## antigen transition that possibly causes changes in mv state ##
    
    x_antigen = antigen_coordinates[loc,0] * np.cos(antigen_coordinates[loc,1])
    y_antigen = antigen_coordinates[loc,0] * np.sin(antigen_coordinates[loc,1])
    #print("loc",loc)
    d = np.sqrt( (x_antigen - coordinates[:,0])**2 + (y_antigen - coordinates[:,1])**2)
    #print(d)
    MV_activation = (d <= r)
    
    # get the mv that covers the antigen transition
    loc2 = np.where(MV_activation==True)[0]
    #print ("loc2",loc2)
    # find the states of antigen covered in the mv that transitioned
    if True in MV_activation:
        x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
        y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])
        d = np.sqrt( (x_antigen - coordinates[loc2,0])**2 + (y_antigen - coordinates[loc2,1])**2)
        #print(d)
        MV_activation = (d <= r)
#        print(MV_activation)
#        print(len(antigen_coordinates[np.where(MV_activation==True),2]==13))
#        print(len(antigen_coordinates))
#        print(antigen_coordinates[np.where(MV_activation==True),2]==13)
        
        #if not stabilized, the mv is removed from the possibility of becoming stabilized#
        # if the sum is equal to one, then that means the only antigen/tcr signalling is dissociating #
        if np.sum(antigen_coordinates[np.where(MV_activation==True),2]==N+3) == 1 and coordinates[loc2,2]==N+4:
            #print("hi")
            P[N+4,0] = P[N+4,0] - 1
            P[1,0] = P[1,0] + 1
            coordinates[loc2,2]=1
        # if stabilized and signalling, the mv now has the possibility of being destabilized
        elif np.sum(antigen_coordinates[np.where(MV_activation==True),2]==N+3) == 1 and coordinates[loc2,2]==N+5:
            #print("hi")
            P[N+5,0] = P[N+5,0] - 1
            P[N+6,0] = P[N+6,0] + 1
            coordinates[loc2,2]=N+6
        #if antigen_coordinates[np.where(MV_activation==True),2]==13:
            
        #print(int(ant_zero))
        
        #P[14,0] = P[14,0] + 1
        
        #print("P",P)
    #coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])
    
    return (coordinates,P)

#def MV_antigen_dist6(antigen_coordinates, coordinates, r, loc, P):
#    ## Only a MV transition, no antigen transition ##
#    
#    # Get locations of MV capable of stabilizing #
#    
#    loc2 = np.where(coordinates[:,2]==14)
#    
#    
#    print ("loc2",loc2)
#    # find the states of antigen covered in the mv that transitioned
#    if True in MV_activation:
#        x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
#        y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])
#        d = np.sqrt( (x_antigen - coordinates[loc2,0])**2 + (y_antigen - coordinates[loc2,1])**2)
#        #print(d)
#        MV_activation = (d <= r)
#        #print("P",P)
#        #print(np.sum(antigen_coordinates[np.where(MV_activation==True),2]==13))
#        if np.sum(antigen_coordinates[np.where(MV_activation==True),2]==13) == 1:
#            #print("hi")
#            P[14,0] = P[14,0] - 1
#            P[1,0] = P[1,0] + 1
#            coordinates[loc2,2]=1
#        #if antigen_coordinates[np.where(MV_activation==True),2]==13:
#            
#        #print(int(ant_zero))
#        
#        #P[14,0] = P[14,0] + 1
#        
#        #print("P",P)
#    #coordinates[-1,2] = coordinates[-1,5] + 1
#    #print("state = ", coordinates[-1,2])
#    
#    return (coordinates,P)

#
#def signal_accumulator(P, tau, accumulator, k_act, t_signal):
#    ### not accounting for the change back ??? ###
#    ### not accumulating signal while in stage 3 ###
#    ### make something cool happen when t_cell is activated ###
#
#    #loc4 = np.where((coordinates[:,2] == 4) == True)[0]
#    #loc3 = np.where((coordinates[:,2] == 3) == True)[0]
#
#    ### first, add up the nodes that are already accumulating signal ###
#    ### 0 indicates the node has been accumulating signal ###
##    signaling_nodes4 = np.sum((coordinates[loc4,3] == 0))
##    signaling_nodes3 = np.sum((coordinates[loc3,3] == 0))
#    signaling_nodes = P[13,0]
#    #print(accumulator)
#    ### for single mv ###
#    ### there can be very large time jumps so must increment in small steps ###
##    if t_signal > 0:
##        print("hi")
##        l = 0
##        while l < t_signal:
##        #print("t_signal =",t_signal)
##            accumulator = accumulator + 0.01*(signaling_nodes)*(k_act * (1-accumulator))
##            l = l + 0.01
##    else: 
##        accumulator = accumulator + tau*(signaling_nodes)*(k_act * (1-accumulator))
#    #accumulator = accumulator + tau*(signaling_nodes3 + signaling_nodes4)*(k_act * (1-accumulator))
#    
#    ### might be a better approach, take into account nodes that were just shut off ###
#    accumulator = accumulator + tau*(signaling_nodes)*(k_act * (1-accumulator))
#    #print(tau)
#    #print(accumulator)
#    #print("3 and 4", signaling_nodes3,signaling_nodes4)
#    #print("tau = ", tau)
#    #print("accumulator =",accumulator)
#    ### second, turn on new signaling nodes ###
#    ### -1 indicates a node was previously not accumulating signal ###
#    if len(loc4 != 0):
#        for i in loc4:
#            if coordinates[i,3] == -1:
#                coordinates[i,3] = 0
#    if len(loc3 != 0):
#        for i in loc3:
#            if coordinates[i,3] == -1:
#                coordinates[i,3] = 0
#
##    non_signaling_nodes4 = np.where((coordinates[loc4,3] == -1) == True)[0]
##    print("non_signal_nodes4", non_signaling_nodes4)
##    coordinates[non_signaling_nodes4,3] = 0
#
#    #non_signaling_nodes3 = np.where((coordinates[loc3,3] == -1))
#    #coordinates[non_signaling_nodes3,3] = 0
#    #print(coordinates)
#
#    for i in range(len(coordinates)):
#        if coordinates[i,2] != 3 and coordinates[i,2] != 4 and coordinates[i,3] == 0:
#            coordinates[i,3] = -1
#
#    return(coordinates, accumulator)
    
#def signal_accumulator_dir(P, tau, accumulator, k_act):
#    ### not accounting for the change back ??? ###
#    ### not accumulating signal while in stage 3 ###
#    ### make something cool happen when t_cell is activated ###
#
#    signaling_mv = P[4,0]+P[5,0]
#    #signaling_mv = P[4,0]
#
#    ### first, add up the nodes that are already accumulating signal ###
#    ### 0 indicates the node has been accumulating signal ###
##    signaling_nodes4 = np.sum((coordinates[loc4,3] == 0))
##    signaling_nodes3 = np.sum((coordinates[loc3,3] == 0))
#    #signaling_nodes = np.sum((coordinates[:,3] == 0))
#    #print(accumulator)
#    ### for single mv ###
#    ### there can be very large time jumps so must increment in small steps ###
#    if tau > 1.0:
#        #print("hi")
#        #print(tau)
#        l = 0
#        while l < tau:
#        #print("t_signal =",t_signal)
#            accumulator = accumulator + 0.01*(signaling_mv)*(k_act * (1-accumulator))
#            l = l + 0.01
#    else: 
#        accumulator = accumulator + tau*(signaling_mv)*(k_act * (1-accumulator))
#    #accumulator = accumulator + tau*(signaling_nodes3 + signaling_nodes4)*(k_act * (1-accumulator))
#    
#    ### might be a better approach, take into account nodes that were just shut off ###
#    #accumulator = accumulator + tau*(signaling_mv)*(k_act * (1-accumulator))
#    #print(tau)
#    #print(accumulator)
#    #print("3 and 4", signaling_nodes3,signaling_nodes4)
#    #print("tau = ", tau)
#    #print("accumulator =",accumulator)
#    ### second, turn on new signaling nodes ###
#    ### -1 indicates a node was previously not accumulating signal ###
##    if len(loc4 != 0):
##        for i in loc4:
##            if coordinates[i,3] == -1:
##                coordinates[i,3] = 0
##    if len(loc3 != 0):
##        for i in loc3:
##            if coordinates[i,3] == -1:
##                coordinates[i,3] = 0
#
##    non_signaling_nodes4 = np.where((coordinates[loc4,3] == -1) == True)[0]
##    print("non_signal_nodes4", non_signaling_nodes4)
##    coordinates[non_signaling_nodes4,3] = 0
#
#    #non_signaling_nodes3 = np.where((coordinates[loc3,3] == -1))
#    #coordinates[non_signaling_nodes3,3] = 0
#    #print(coordinates)
#
##    for i in range(len(coordinates)):
##        if coordinates[i,2] != 3 and coordinates[i,2] != 4 and coordinates[i,3] == 0:
##            coordinates[i,3] = -1
#
#    return(accumulator)
def signal_accumulator_dir(P, tau, accumulator, k_act, gam, N, Xact, t):
    #number of activated TCRs
    signaling_mv = P[N+3,0]


    ### be careful of blowups. the update rule for signal only works in small time increments.
    ### break up larger time increments to smaller time increments
    if tau > 0.1:
        #print("hi")
        #print(tau)
        l = 0
        while l < tau and accumulator < Xact:
        #added a decay 0.05
            #accumulator = accumulator - 0.0001*0.05*accumulator
            accumulator = accumulator + 0.01*(signaling_mv)*(k_act * (1.0-accumulator)) - 0.01*gam*accumulator
            #accumulator = accumulator + 0.0001*(signaling_mv)
            l = l + 0.01
            t = t + 0.01
    else:
        #accumulator = accumulator*(1.0 - tau*0.05)
        accumulator = accumulator + tau*(signaling_mv)*(k_act * (1.0-accumulator)) - tau*gam*accumulator
    #print("accum after update", accumulator)
    return(accumulator,t)
#def signal_accumulator_dir(P, tau, accumulator, k_act):
#    ### not accounting for the change back ??? ###
#    ### not accumulating signal while in stage 3 ###
#    ### make something cool happen when t_cell is activated ###
#
#    signaling_mv = np.sum(P[42:82,0])
#
#    ### first, add up the nodes that are already accumulating signal ###
#    ### 0 indicates the node has been accumulating signal ###
##    signaling_nodes4 = np.sum((coordinates[loc4,3] == 0))
##    signaling_nodes3 = np.sum((coordinates[loc3,3] == 0))
#    #signaling_nodes = np.sum((coordinates[:,3] == 0))
#    #print(accumulator)
#    ### for single mv ###
#    ### there can be very large time jumps so must increment in small steps ###
#    if tau > 1.0:
#        #print("hi")
#        #print(tau)
#        l = 0
#        while l < tau:
#        #print("t_signal =",t_signal)
#            accumulator = accumulator + 0.01*(signaling_mv)*(k_act * (1-accumulator))
#            l = l + 0.01
#    else: 
#        accumulator = accumulator + tau*(signaling_mv)*(k_act * (1-accumulator))
#    #accumulator = accumulator + tau*(signaling_nodes3 + signaling_nodes4)*(k_act * (1-accumulator))
#    
#    ### might be a better approach, take into account nodes that were just shut off ###
#    #accumulator = accumulator + tau*(signaling_mv)*(k_act * (1-accumulator))
#    #print(tau)
#    #print(accumulator)
#    #print("3 and 4", signaling_nodes3,signaling_nodes4)
#    #print("tau = ", tau)
#    #print("accumulator =",accumulator)
#    ### second, turn on new signaling nodes ###
#    ### -1 indicates a node was previously not accumulating signal ###
##    if len(loc4 != 0):
##        for i in loc4:
##            if coordinates[i,3] == -1:
##                coordinates[i,3] = 0
##    if len(loc3 != 0):
##        for i in loc3:
##            if coordinates[i,3] == -1:
##                coordinates[i,3] = 0
#
##    non_signaling_nodes4 = np.where((coordinates[loc4,3] == -1) == True)[0]
##    print("non_signal_nodes4", non_signaling_nodes4)
##    coordinates[non_signaling_nodes4,3] = 0
#
#    #non_signaling_nodes3 = np.where((coordinates[loc3,3] == -1))
#    #coordinates[non_signaling_nodes3,3] = 0
#    #print(coordinates)
#
##    for i in range(len(coordinates)):
##        if coordinates[i,2] != 3 and coordinates[i,2] != 4 and coordinates[i,3] == 0:
##            coordinates[i,3] = -1
#
#    return(accumulator)
    
def lt_bin_calc(mv_lifetime,mv_lifetime_index):
    lt_minus = mv_lifetime_index - mv_lifetime
    index = np.where(lt_minus == np.min(lt_minus[lt_minus > 0]))
    return (index)
    
