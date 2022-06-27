# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 16:21:36 2019

@author: Jonathan
"""
import numpy as np
from initial_mv_conditions_v2 import *
from functions_mv_v2 import *
from numpy import linalg as LA
import time
import h5py

def mv_sim_num_direct( antigen_count, k_off, k_s, Nkp, k_d ):
    
    # running sum of the concentration of mv in stabilized states #
    #uconc = 0
    #uconc_list = np.zeros((1,1))
    #tau_list = np.zeros((1,1))
    mv_lt = np.array([])
    Xmax = 0.0
    ### array for recording average area covered by MV ###
    #area_event = np.zeros((1,1))
    #area_event[0] = 0
    
    # a variable for storing accumulated signal #
    accumulator = 0.0
    
    # an array for collecting lifetime intervals #
    # used for histogram of MV lifetimes #
    #mv_lifetime_index = np.cumsum(np.ones((1,340)) / 4.0)

    lt_bin = np.zeros((1,340))
    
    
    r_is = 1.0
    r = 0.05345224838248488
    #r = 0.0313727
    
    # This is the average number of MV found with a set radius that fill the immunological synapse #
    # This is determined by multiple simulations of randomly placing MV in a circle #
    ave_max = 170
    #ave_max = 512
    
    coordinates = mv0(r, ave_max)
    
    ### Add antigen to simulation ###
    
    num_ant = int(antigen_count)
    
    antigen_coordinates = np.zeros((num_ant,3))
   
    
    i = 0
    # distribute the antigen across the immunological synapse #
    while i < num_ant:
    
            r_antigen,theta = antigen_antigen_dist(antigen_coordinates, r)
    
            antigen_coordinates[i,0] = r_antigen
            
            antigen_coordinates[i,1] = theta
            
            i = i + 1
    #print("len", len(antigen_coordinates))
    ### Find antigen in contact proximity of MV and change the state vector ###       
    antigen_coordinates = MV_antigen_dist2(antigen_coordinates, coordinates, r)
    #print(len(antigen_coordinates))
    #print("len", len(antigen_coordinates))
    ### Setting the relevant rates ###
    kd = k_d
    kon = 0.1

    ka = kd * 4.666667
    kt = 0.1

    koff = k_off
    ku = 0.03
    ks = k_s
    k_act = 0.2
    
    N=Nkp
    ### Initializing the population vector ###
    P = np.zeros((N+7,1))
    
    # Set a separate state for every group of MVs covering a specific number of antigen #
    #for i in range(len(P[0:23]) - 1):
        # the fifth column of the coordinate matrix denotes the number of antigen within binding range of the MV #
        #P[i+1][0] = np.sum(coordinates[:,5] == 1)
        
    # The number of MV still capable of being added to the simulation #
    
    P[1,0] = len(coordinates)
    P[0][0] = ave_max - P[1,0]
    P[2,0] = len(np.where(antigen_coordinates[:,2]>=2)[0])
    #print(len(np.where(antigen_coordinates[:,2]>=2)[0]))
    #print(np.where(antigen_coordinates[:,2]>=2))
    #print("im here",antigen_coordinates)

    
    # set the coordinate states according to their number of antigen covered #
    #coordinates[:,2] = coordinates[:,5] + 1

    
    
    # T is the transition matrix, sort of # 
    # MV that cover a certain number of antigen are treated as having their own pathway. This is because their kon rates will differ #
    # They are given their own pathway in order to track placement and transitions. For example, a MV that covers two antigen cannot #
    # transition to a state where it only covers one antigen. This is because the antigen in this model are fixed #
    # The simulation allows for the possibility of 20 antigen being covered. The placement of antigen takes this into consideration #
#    T = np.zeros((482,482))
#    T[0][1] = ka
#    T[1][0] = kd
#    for i in range(20):
#        T[i+2][0] = kd
#        T[i+2][i+20+2] = kon*(i+1)
#    for i in range(20):
#        T[i+2+20][0] = kd
#        T[i+2+20][i+2] = koff
#        T[i+2+20][i+2+40] = kt
#    for i in range(20):
#        T[i+2+40][0] = kd
#        T[i+2+40][i+2] = koff
#        T[i+2+40][i+2+60] = kt
#    for i in range(20):
#        T[i+2+60][0] = kd
#        T[i+2+60][i+2] = koff
#        T[i+2+60][i+2+80] = kt
#    for i in range(20):
#        T[i+2+80][0] = kd
#        T[i+2+80][i+2] = koff
#        T[i+2+80][i+2+100] = kt
#    for i in range(20):
#        T[i+2+100][0] = kd
#        T[i+2+100][i+2] = koff
#        T[i+2+100][i+2+120] = kt
#    for i in range(20):
#        T[i+2+120][0] = kd
#        T[i+2+120][i+2] = koff
#        T[i+2+120][i+2+140] = kt
#    for i in range(20):
#        T[i+2+140][0] = kd
#        T[i+2+140][i+2] = koff
#        T[i+2+140][i+2+160] = kt
#    for i in range(20):
#        T[i+2+160][0] = kd
#        T[i+2+160][i+2] = koff
#        T[i+2+160][i+2+180] = kt
#    for i in range(20):
#        T[i+2+180][0] = kd
#        T[i+2+180][i+2] = koff
#        T[i+2+180][i+2+200] = kt
#    for i in range(20):
#        T[i+2+200][0] = kd
#        T[i+2+200][i+2] = koff
#        T[i+2+200][i+2+220] = kt
#    for i in range(20):
#        T[i+2+220][0] = kd
#        T[i+2+220][i+2] = koff
#        T[i+2+220][i+2+240] = ks
#    for i in range(20):
#        T[i+2+240][i+2+260] = koff
#    for i in range(20):
#        T[i+2+260][i+2] = ku
#        T[i+2+260][i+2+280] = kon2*(i+1)
#    for i in range(20):
#        T[i+2+280][i+2+20] = ku
#        T[i+2+280][i+2+300] = kt2
#        T[i+2+280][i+2+260] = koff
#    for i in range(20):
#        T[i+2+300][i+2+40] = ku
#        T[i+2+300][i+2+320] = kt2
#        T[i+2+300][i+2+260] = koff
#    for i in range(20):
#        T[i+2+320][i+2+60] = ku
#        T[i+2+320][i+2+340] = kt2
#        T[i+2+320][i+2+260] = koff
#    for i in range(20):
#        T[i+2+340][i+2+80] = ku
#        T[i+2+340][i+2+360] = kt2
#        T[i+2+340][i+2+260] = koff
#    for i in range(20):
#        T[i+2+360][i+2+100] = ku
#        T[i+2+360][i+2+380] = kt2
#        T[i+2+360][i+2+260] = koff
#    for i in range(20):
#        T[i+2+380][i+2+120] = ku
#        T[i+2+380][i+2+400] = kt2
#        T[i+2+380][i+2+260] = koff
#    for i in range(20):
#        T[i+2+400][i+2+140] = ku
#        T[i+2+400][i+2+420] = kt2
#        T[i+2+400][i+2+260] = koff
#    for i in range(20):
#        T[i+2+420][i+2+160] = ku
#        T[i+2+420][i+2+440] = kt2
#        T[i+2+420][i+2+260] = koff
#    for i in range(20):
#        T[i+2+440][i+2+180] = ku
#        T[i+2+440][i+2+460] = kt2
#        T[i+2+440][i+2+260] = koff
#    for i in range(20):
#        T[i+2+460][i+2+200] = ku
#        T[i+2+460][i+2+240] = kt2
#        T[i+2+460][i+2+260] = koff
#    
    #P = np.zeros((16,1))

    # The number of TCRs in a simulation. In the control, this is equal to the number of antigen #
    #P[0][0] = av
    #P[15][0] = int(antigen_count)
    coordinates[:,2]=1
    
    T = np.zeros((N+7,N+7))
    # For MV dynamics or some other circumstance where TCRs can be removed or not be engaged with antigen #
    T[0,1] = ka
    T[1,0] = kd
    #T[1,14] = ks
    #T[2,0] = kd
    T[2,3] = kon
    i=3
    while i < 3+N:
        T[i,2] = koff
        T[i,i+1] = kt
        i=i+1

        
    T[N+3,2] = koff
    T[N+4,N+5] = ks
    T[N+4,0] = kd
    T[N+6,1] = ku

    
    t = 0
    coordinates[:,5] = 0.0
    gam = 0.0
    maxT = k_act/(k_act + gam)
    xact = 0.90*maxT
    while t < 600 and accumulator < xact:
#    while t < 60:
        # The Direct Gillespie Algorithm #
        # The first random number #
        r1 = np.random.uniform(0,1)

        # P is the population vector of MV in all states #
        # T is the transition matrix #
        # P is broadcasted, element wise onto each row of T #
        # Every nonzero element in R is where transitions may occur #
        # Rij = Pi * Tij where Pi is the number of MV capable of the transition from state i to state j #                               
        R = P*T
        #print("P",P)
        # Cumulative sum of the elements in R #
        # (R11, R11 + R12, ..., R11 + R12 + ... + R21, ...) #
        # Rcs is a vector in R[122*122] #
        Rcs = np.cumsum(R)
        
        # Reshape Rcs into matrix #
        Rcs = np.reshape(Rcs,(N+7,N+7))
        
        # since Rcs has many zeros, there are a lot of repeat elements from the cumulative sum #
        # set those elements back to zero #
        Rcs[R == 0] = 0

        # Get the 1-norm of all transition elements #
        w1 = np.sum(np.absolute(R))
        
        # generate the next time #
        tau = 1.0 / w1 * np.log(1.0/r1)
        
        if t + tau > 600.0:
            #print("hi2")
            tau = 600.0 - t
            #print(tau)
            #break
        
        # Second random number #
        r2 =  np.random.uniform(0,1)
        
        # Find which event occured #
        # This is done by subtracting the event value (r2*w1) from all elements in the cumulative sum matrix #
        # The new matrix Rminus will contain one minimum positive value #
        # The matrix coordinates of this value is the transition that occurred #
        # For example, if the location of the minimum positive value is located in the ith row and jth column, #
        # then a MV has transitioned from the ith state to the jth state #
        Rminus = Rcs - r2*w1
        transit = np.where(Rminus == np.min(Rminus[Rminus > 0]))
    
        # update the amount of signal accumulated #
        #print("accum before function =",accumulator)
        (accumulator, tstar) = signal_accumulator_dir(P, tau, accumulator, k_act, gam, N, xact, t)
        if accumulator > Xmax:
            Xmax = accumulator
            
        t = t + tau
        if accumulator > xact:
            t = tstar
        #print("accum after function =",accumulator)
        #accumulator = XX + 0.0
        
        # the 4th column is the individual MV lifetime #
        # add the next reaction time to the MV lifetimes #
        coordinates[:,4] = coordinates[:,4] + tau
        
        ### Update the State Matrix ###
        # Handling the addition of a new MV #
        # transit[0][0] is the state a MV is transitioning from #
        # the first element in P is the population of MV capable of being placed in the IS (immunological synapse) #
        # if transit[0][0] is zero, then a new MV is appearing on the surface #
        if transit[0][0] == 0:
            
            # position the MV using the gen_coords function #
            # Takes in the arguments of: IS radius r_is, MV radius r, and the state matrix coordinates #
            (x_new,y_new) = gen_coords(coordinates, r_is, r)

            # update the state matrix with new MV coordinates #
            coordinates = np.append(coordinates,[[x_new,y_new, 1, -1, 0, 0]], axis = 0)
            
            # Determine the number of antigen within binding distance of new MV#
            antigen_coordinates,P = MV_antigen_dist1(antigen_coordinates, coordinates, r,P)

            # update the population of states by subtracting 1 from where the MV came from #
            # and adding one to the state where the MV transitioned to #
            P[0,0] = P[0,0] - 1
    
            # the second column in the state matrix corresponds to the individual states of each MV #
            # in the case of adding a new MV to the IS, the state of the new MV is dependent upon #
            # how many antigen are covered, if any #
            P[1,0] = P[1,0] + 1

            
        # updating the state with new transition #
        # not coming from zero state and not going to zero state #
        # also, not worrying about active and inactive MV#
        elif transit[1][0] > 0 and transit[0][0] > 1:

            # a random number for every antigen that is capable of the transition
            #print("transit",transit)
#            ru = np.random.uniform(0,1, size = (int(np.sum(antigen_coordinates[:,2] == transit[0][0]))))
##            print("transit",transit)
##            print("ru",ru)
##            print("antigen coord",antigen_coordinates)
##            print("pop",P)
#            # find the max value #
#            rmax = np.where(ru == np.max(ru))[0]

            # get the location of antigen capable of making the transition #
            #loc_mv = np.where(antigen_coordinates[:,2] == transit[0][0])[0]
            
            # an antigen transition that allows MV to become stabilized, so possible 2 transitions, if mv is not in pre-stab
            if transit[1][0] == N+3:
                ru = np.random.uniform(0,1, size = (int(np.sum(antigen_coordinates[:,2] == transit[0][0]))))

                rmax = np.where(ru == np.max(ru))[0]
                loc_mv = np.where(antigen_coordinates[:,2] == transit[0][0])[0]
                coordinates,P = MV_antigen_dist4(antigen_coordinates, coordinates, r, loc_mv[rmax],P, N)
                antigen_coordinates[loc_mv[rmax],2] = transit[1][0]
                
                # inside MVant4, only mv states were adjusted, now adjust tcr states #
                P[transit[0][0],0] = P[transit[0][0],0] - 1
                P[transit[1][0],0] = P[transit[1][0],0] + 1
                
            # possible antigen transition that removes mv from stab intermediate #
            elif transit[0][0] == N+3:
                #print("transit",transit)
                #print('why')
                ru = np.random.uniform(0,1, size = (int(np.sum(antigen_coordinates[:,2] == transit[0][0]))))

                rmax = np.where(ru == np.max(ru))[0]
                loc_mv = np.where(antigen_coordinates[:,2] == transit[0][0])[0]
                #print("hi")
                
                # update MV information
                coordinates,P = MV_antigen_dist5(antigen_coordinates, coordinates, r, loc_mv[rmax],P, N)
                
                # now update new antigen information
                antigen_coordinates[loc_mv[rmax],2] = transit[1][0]
                P[transit[0][0],0] = P[transit[0][0],0] - 1
                P[transit[1][0],0] = P[transit[1][0],0] + 1
                #print("P2",P)
                
            # MV contact is stabilized #
            ###  ###
            elif transit[0][0] == N+4 and transit[1][0] == N+5:
                #print("transit",transit)
                #print(P)
                ru = np.random.uniform(0,1, size = (int(np.sum(coordinates[:,2] == transit[0][0]))))

                rmax = np.where(ru == np.max(ru))[0]
                loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]
                
                #print("hi")
                #coordinates,P = MV_antigen_dist5(antigen_coordinates, coordinates, r, loc_mv[rmax],P)
                coordinates[loc_mv[rmax],2] = transit[1][0]
                P[transit[0][0],0] = P[transit[0][0],0] - 1
                P[transit[1][0],0] = P[transit[1][0],0] + 1
                #print("P2",P)
                
            # MV is removed from intermediate stab state #
#            elif transit[0][0] == 14 and transit[1][0] == 1:
#                #print("transit",transit)
#                #print(P)
#                ru = np.random.uniform(0,1, size = (int(np.sum(coordinates[:,2] == transit[0][0]))))
#
#                rmax = np.where(ru == np.max(ru))[0]
#                loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]
#                
#                #print("hi")
#                #coordinates,P = MV_antigen_dist5(antigen_coordinates, coordinates, r, loc_mv[rmax],P)
#                coordinates[loc_mv[rmax],2] = transit[1][0]
#                P[transit[0][0],0] = P[transit[0][0],0] - 1
#                P[transit[1][0],0] = P[transit[1][0],0] + 1
#                #print("P2",P)
#                
#                
            # mv is destabilized, need to see if there are any deactivated TCRs underneath
            elif transit[0][0] == N+6:
                #print("transit",transit)
                #print(P)
                ru = np.random.uniform(0,1, size = (int(np.sum(coordinates[:,2] == transit[0][0]))))

                rmax = np.where(ru == np.max(ru))[0]
                loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]
                #print("hi")
                #coordinates,P = MV_antigen_dist5(antigen_coordinates, coordinates, r, loc_mv[rmax],P)
                coordinates[loc_mv[rmax],2] = transit[1][0]
                P[transit[0][0],0] = P[transit[0][0],0] - 1
                P[transit[1][0],0] = P[transit[1][0],0] + 1
                #print("P2",P)
            # just an antigen transition    
            else:
                ru = np.random.uniform(0,1, size = (int(np.sum(antigen_coordinates[:,2] == transit[0][0]))))

                rmax = np.where(ru == np.max(ru))[0]
                loc_mv = np.where(antigen_coordinates[:,2] == transit[0][0])[0]
                antigen_coordinates[loc_mv[rmax],2] = transit[1][0]
            
                P[transit[0][0],0] = P[transit[0][0],0] - 1
                P[transit[1][0],0] = P[transit[1][0],0] + 1
            
#        elif transit[1][0] == 0 and transit[0][0]==14: 
#            # a random number for every antigen that is capable of the transition
#            ru = np.random.uniform(0,1, size = (int(np.sum(antigen_coordinates[:,2] == transit[0][0]))))
##            print("transit",transit)
##            print("ru",ru)
##            print("antigen coord",antigen_coordinates)
##            print("pop",P)
#            # find the max value #
#            rmax = np.where(ru == np.max(ru))[0]
#
#            # get the location of antigen capable of making the transition #
#            loc_mv = np.where(antigen_coordinates[:,2] == transit[0][0])[0]
#            
#
#            # just a MV transition
#            elif transit[0][0] == 14:
#                print("hi")
#                coordinates,P = MV_antigen_dist5(antigen_coordinates, coordinates, r, loc_mv[rmax],P)
        # handling the removal of an mv #       
        elif transit[1][0] == 0: 
            
            ### Active MV can either be in state 1 (active) or 2 (stabilized) ###
            ru = np.random.uniform(low=0.0, high=1.0, size=(int(np.sum(coordinates[:,2] == transit[0][0]))))
            #print("transit",transit)
#            print("ru",ru)
#            print("MV coord",coordinates[:,2]==1)
#            print("pop",P)
            loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]
            rmax = np.where(ru == np.max(ru))[0]
            
            # record the lifetime of the removed MV #
            #mv_lifetime = coordinates[loc_mv[rmax],4]
            antigen_coordinates,P = MV_antigen_dist3(antigen_coordinates, coordinates, r, loc_mv[rmax],P, N)
            mv_lt = np.append(mv_lt,coordinates[loc_mv[rmax],5])
            coordinates = np.delete(coordinates, loc_mv[rmax],0)
            #antigen_coordinates = MV_antigen_dist1(antigen_coordinates, coordinates, r)
            P[transit[0][0],0] = P[transit[0][0],0] - 1
            P[transit[1][0],0] = P[transit[1][0],0] + 1
        
            # Record MV lifetimes into bins with a width of 0.25s #
            #lt_index = lt_bin_calc(mv_lifetime,mv_lifetime_index)
            #lt_bin[0,lt_index] = lt_bin[0,lt_index] + 1
          
        # a mv is being removed from the state of pre-stab    
#        elif transit[1][0] == 0 and transit[0][0] == 5: 
#            
#            ### Active MV can either be in state 1 (active) or 2 (stabilized) ###
#            ru = np.random.uniform(low=0.0, high=1.0, size=(int(np.sum(coordinates[:,2] == 5))))
#            #print("transit",transit)
##            print("ru",ru)
##            print("MV coord",coordinates[:,2]==1)
##            print("pop",P)
#            loc_mv = np.where(coordinates[:,2] == 5)[0]
#            rmax = np.where(ru == np.max(ru))[0]
#            
#            # record the lifetime of the removed MV #
#            #mv_lifetime = coordinates[loc_mv[rmax],4]
#            antigen_coordinates,P = MV_antigen_dist3(antigen_coordinates, coordinates, r, loc_mv[rmax],P)
#            coordinates = np.delete(coordinates, loc_mv[rmax],0)
#            #antigen_coordinates = MV_antigen_dist1(antigen_coordinates, coordinates, r)
#            P[transit[0][0],0] = P[transit[0][0],0] - 1
#            P[transit[1][0],0] = P[transit[1][0],0] + 1
        
            # Record MV lifetimes into bins with a width of 0.25s #
            #lt_index = lt_bin_calc(mv_lifetime,mv_lifetime_index)
            #lt_bin[0,lt_index] = lt_bin[0,lt_index] + 1
        
        
        
        coordinates[:,5] = coordinates[:,5] + tau
        #print(np.sum(P[2:14,0]))
        #tau_list = np.append(tau_list,tau)
#        if P[13,0]>0 and P[14,0]+P[15,0]==0:
#            print(P)
#        if P[16,0]==1:
#            print(P)
        #print(accumulator)
        #print(tau_list[1])
        # Record the instantaneous area covered by MV #
        #area_event = np.append(area_event,r**2*np.sum(P[1:8,0]))
         

        #uconc_list = np.append(uconc_list,(np.sum(P[102:122,0]) + np.sum(P[122:142,0])))
    
    #uconc = 0
    #i = 0
    #while i < len(uconc_list) - 1:
        #uconc = uconc + uconc_list[i] * tau_list[i+1]
        #i = i + 1
    #uconc = uconc / 60.0
    #area_event = np.sum(area_event) / len(area_event)
    #mv_lt = np.append(mv_lt,coordinates[:,5])
    #mv_lt_ave = np.sum(mv_lt)/len(mv_lt)
    return (accumulator, t, Xmax, mv_lt)

hm_matrix = np.zeros((2,80))
nkp = 3
kd = np.logspace(-3.0, -1.0, num=60)
act_sum = 0
koffi = np.logspace(-2.0, -1.5, num=10)
koff = 0.01
ind = 1
for k in koffi:
  for j in range(len(kd)):
      t_act = 0
      trial = 0
      act_sum = 0
      while trial < 50:
          (x, t, Xmax, mv_lt) = mv_sim_num_direct( 1,  k, 100.0, nkp, kd[j] )
          #(x,c,d,af2,cf2,t2)=mv_sim_num_direct( 0, 50.5, 1,0.0 )
          t_act = t_act + t
          trial = trial + 1
          if x >= 0.9:
              act_sum = act_sum + 1
      #    fig = plt.figure()
      #    ax = plt.subplot(111)
      #    ax.plot(t, af, label='Instantaneous Area Fraction')
      #    ax.plot(t, cf, label='Cumulative Area Fraction IS')
      #    ax.plot(t2, cf2, label='Cumulative Area Fraction')
      #    ax.axes.set_xlim([0,60])
      #    ax.axes.set_ylim([0,1])
      #    plt.title('Cumulative and Instantaneous Area Fractions')
      #    plt.ylabel('Fraction of Area Covered')
      #    plt.xlabel('s')
      #    ax.legend()
      #    plt.show()
      act_sum = act_sum / 100.0
      t_act = t_act / 100.0
      print( kd[j],t_act,act_sum)
      hm_matrix[0,j] = t_act
      hm_matrix[1,j] = act_sum 
  
  np.save('n1_hm_mat'+str(ind)+'.npy', hm_matrix)
  koff = koff + 0.1
  ind = ind + 1
#plt.hist(mv_lt, bins = 20) 
#plt.show()
#print(np.sum(mv_lt)/len(mv_lt))
#print(np.shape(pop_sig))
#antigen = [1, 2, 4, 6, 8, 12, 18, 26, 37, 54, 78, 112, 162, 233, 335, 483, 695, 1000]
#antigen = [0.01, 0.1, 1.0, 10] 
#f = h5py.File("act_tcr_MV_stab_t60koff001.hdf5", "w")
#f.create_dataset("act_tcr_ant10000", data = pop_sig)
##f.create_dataset("X", data = X)
##for j in range(len(antigen)):
##    i=0
##    while i < 1:
##        (X, t, popdis, act_tcr) = mv_sim_num_direct( 1000, antigen[j], 100.0 )
##        #act_tcr_array[j] = act_tcr[0]
##        i = i + 1
##        print(act_tcr.shape)
##        print(j)
##        plt.plot(act_tcr[0]) 
##        plt.show()
##        
##    f.create_dataset("act_tcr_ant"+str(antigen[j]), data = act_tcr)
##                     
##params = ['kd = 1/7','kon = 100', 'ka = kd*4.67','kt = 3.5', 'koff = 0.75','r = 0.0534','not averaged']
##asciiList = [n.encode("ascii", "ignore") for n in params]
##f.create_dataset('params', (len(asciiList),1),'S20', asciiList)
##
#f.close()



