# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 08:55:01 2019

@author: Jonathan
"""

import numpy as np
from initial_mv_conditions_v2 import *
from cumulative_functions import *
from numpy import linalg as LA
#import time
#import h5py

import matplotlib.pyplot as plt
def mv_sim_num_direct( antigen_count, k_off, k_s, k_on, rr, avemax, t_sim):
    o=0
    # running sum of the concentration of mv in stabilized states #
    #uconc = 0
    #uconc_list = np.zeros((1,1))
    X_list = np.array([])
    tau_list = np.zeros((1,500000))
    tau_list2 = np.zeros((1,500000))
    pop = np.zeros((16,500000))
    pop_loc = np.zeros((1,500000))
    t_list=np.array([])
    af_list=np.array([])
    cf_list=np.array([])
    ### array for recording average area covered by MV ###
    #area_event = np.zeros((1,1))
    #area_event[0] = 0
    rebind_time = np.array([])
    # a variable for storing accumulated signal #
    accumulator = 0.0

    # an array for collecting lifetime intervals #
    # used for histogram of MV lifetimes #
    #mv_lifetime_index = np.cumsum(np.ones((1,340)) / 4.0)

    lt_bin = np.zeros((1,340))


    r_is = 1.0
    r = rr

    # This is the average number of MV found with a set radius that fill the immunological synapse #
    # This is determined by multiple simulations of randomly placing MV in a circle #
    ave_max = int(avemax)

    coordinates = np.array([])
    (x_new,y_new) = gen_coords(coordinates, r_is, r)

    coordinates = mv0(r, ave_max, 4.666667)

    ### Add antigen to simulation ###

    num_ant = int(antigen_count)

    antigen_coordinates = np.zeros((num_ant,2))
    antigen_coordinates2 = np.zeros((2000,2))
    #coordinates = np.array([[x_new,y_new, 0,-1, 0, 0]])
    i = 0
    golden_angle = np.pi * (3 - np.sqrt(5))
    # distribute the antigen across the immunological synapse #
    while i < 2000:

            theta = i * golden_angle
            r_antigen = np.sqrt(i) / np.sqrt(2000)


            #r_antigen,theta = antigen_antigen_dist(antigen_coordinates, r)

            #antigen_coordinates[i,0] = r_antigen

            #antigen_coordinates[i,1] = theta

            antigen_coordinates2[i,0] = r_antigen

            antigen_coordinates2[i,1] = theta

            i = i + 1
    i=0
    while i < num_ant:

            r_antigen,theta = antigen_antigen_dist(antigen_coordinates, r)

            antigen_coordinates[i,0] = r_antigen

            antigen_coordinates[i,1] = theta

            i = i + 1
    ### Find antigen in contact proximity of MV and change the state vector ###
    (coordinates,antigen_coordinates,antigen_coordinates2) = MV_antigen_dist2(antigen_coordinates,antigen_coordinates2, coordinates, r)

    ### Setting the relevant rates ###
    kd = 1/7.0
    kon = k_on
    kon2 = k_on
    ka = kd * 4.666667
    kt = 3.5
    kt2 = 3.5
    koff = k_off
    ku = 0.03
    ks = k_s
    k_act = 2.0

    ### Initializing the population vector ###
    P = np.zeros((482,1))

    # Set a separate state for every group of MVs covering a specific number of antigen #
    for i in range(len(P[0:23]) - 1):
        # the fifth column of the coordinate matrix denotes the number of antigen within binding range of the MV #
        P[i+1][0] = np.sum(coordinates[:,5] == i)

    # The number of MV still capable of being added to the simulation #
    P[0][0] = ave_max - np.sum(P[1:483,0])


    # set the coordinate states according to their number of antigen covered #
    coordinates[:,2] = coordinates[:,5] + 1



    # T is the transition matrix, sort of #
    # MV that cover a certain number of antigen are treated as having their own pathway. This is because their kon rates will differ #
    # They are given their own pathway in order to track placement and transitions. For example, a MV that covers two antigen cannot #
    # transition to a state where it only covers one antigen. This is because the antigen in this model are fixed #
    # The simulation allows for the possibility of 20 antigen being covered. The placement of antigen takes this into consideration #
    T = np.zeros((482,482))
    T[0][1] = ka
    T[1][0] = kd
    for i in range(20):
        T[i+2][0] = kd
        T[i+2][i+20+2] = kon*(i+1)
    for i in range(20):
        T[i+2+20][0] = kd
        T[i+2+20][i+2] = koff
        T[i+2+20][i+2+40] = kt
    for i in range(20):
        T[i+2+40][0] = kd
        T[i+2+40][i+2] = koff
        T[i+2+40][i+2+60] = kt
    for i in range(20):
        T[i+2+60][0] = kd
        T[i+2+60][i+2] = koff
        T[i+2+60][i+2+80] = kt
    for i in range(20):
        T[i+2+80][0] = kd
        T[i+2+80][i+2] = koff
        T[i+2+80][i+2+100] = kt
    for i in range(20):
        T[i+2+100][0] = kd
        T[i+2+100][i+2] = koff
        T[i+2+100][i+2+120] = kt
    for i in range(20):
        T[i+2+120][0] = kd
        T[i+2+120][i+2] = koff
        T[i+2+120][i+2+140] = kt
    for i in range(20):
        T[i+2+140][0] = kd
        T[i+2+140][i+2] = koff
        T[i+2+140][i+2+160] = kt
    for i in range(20):
        T[i+2+160][0] = kd
        T[i+2+160][i+2] = koff
        T[i+2+160][i+2+180] = kt
    for i in range(20):
        T[i+2+180][0] = kd
        T[i+2+180][i+2] = koff
        T[i+2+180][i+2+200] = kt
    for i in range(20):
        T[i+2+200][0] = kd
        T[i+2+200][i+2] = koff
        T[i+2+200][i+2+220] = kt
    for i in range(20):
        T[i+2+220][0] = kd
        T[i+2+220][i+2] = koff
        T[i+2+220][i+2+240] = ks
    for i in range(20):
        T[i+2+240][i+2+260] = koff
    for i in range(20):
        T[i+2+260][i+2] = ku
        T[i+2+260][i+2+280] = kon2*(i+1)
    for i in range(20):
        T[i+2+280][i+2+20] = ku
        T[i+2+280][i+2+300] = kt2
        T[i+2+280][i+2+260] = koff
    for i in range(20):
        T[i+2+300][i+2+40] = ku
        T[i+2+300][i+2+320] = kt2
        T[i+2+300][i+2+260] = koff
    for i in range(20):
        T[i+2+320][i+2+60] = ku
        T[i+2+320][i+2+340] = kt2
        T[i+2+320][i+2+260] = koff
    for i in range(20):
        T[i+2+340][i+2+80] = ku
        T[i+2+340][i+2+360] = kt2
        T[i+2+340][i+2+260] = koff
    for i in range(20):
        T[i+2+360][i+2+100] = ku
        T[i+2+360][i+2+380] = kt2
        T[i+2+360][i+2+260] = koff
    for i in range(20):
        T[i+2+380][i+2+120] = ku
        T[i+2+380][i+2+400] = kt2
        T[i+2+380][i+2+260] = koff
    for i in range(20):
        T[i+2+400][i+2+140] = ku
        T[i+2+400][i+2+420] = kt2
        T[i+2+400][i+2+260] = koff
    for i in range(20):
        T[i+2+420][i+2+160] = ku
        T[i+2+420][i+2+440] = kt2
        T[i+2+420][i+2+260] = koff
    for i in range(20):
        T[i+2+440][i+2+180] = ku
        T[i+2+440][i+2+460] = kt2
        T[i+2+440][i+2+260] = koff
    for i in range(20):
        T[i+2+460][i+2+200] = ku
        T[i+2+460][i+2+240] = kt2
        T[i+2+460][i+2+260] = koff


    t = 0
    n = 0
    gam = 0.005
    maxT = k_act/(k_act + gam)
    while t <= t_sim and accumulator < 1.0:

        # The Direct Gillespie Algorithm #
        # The first random number #
        r1 = np.random.uniform(0,1)

        # P is the population vector of MV in all states #
        # T is the transition matrix #
        # P is broadcasted, element wise onto each row of T #
        # Every nonzero element in R is where transitions may occur #
        # Rij = Pi * Tij where Pi is the number of MV capable of the transition from state i to state j #
        R = P*T

        # Cumulative sum of the elements in R #
        # (R11, R11 + R12, ..., R11 + R12 + ... + R21, ...) #
        # Rcs is a vector in R[122*122] #
        Rcs = np.cumsum(R)

        # Reshape Rcs into matrix #
        Rcs = np.reshape(Rcs,(482,482))

        # since Rcs has many zeros, there are a lot of repeat elements from the cumulative sum #
        # set those elements back to zero #
        Rcs[R == 0] = 0

        # Get the 1-norm of all transition elements #
        w1 = np.sum(np.absolute(R))

        # generate the next time #
        tau = 1.0 / w1 * np.log(1.0/r1)

        if t + tau > t_sim:
            #print("hi2")
            tau = t_sim - t
            #print(tau)
            break

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
        (accumulator) = signal_accumulator_dir(P, tau, accumulator, k_act,gam)
        X_list = np.append(X_list, accumulator)
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
            coordinates = np.append(coordinates,[[x_new,y_new, 0, -1, 0, 0]], axis = 0)

            # Determine the number of antigen within binding distance of new MV#
            (coordinates,antigen_coordinates,antigen_coordinates2) = MV_antigen_dist1(antigen_coordinates,antigen_coordinates2, coordinates, r)

            # update the population of states by subtracting 1 from where the MV came from #
            # and adding one to the state where the MV transitioned to #
            P[transit[0][0],0] = P[transit[0][0],0] - 1

            # the second column in the state matrix corresponds to the individual states of each MV #
            # in the case of adding a new MV to the IS, the state of the new MV is dependent upon #
            # how many antigen are covered, if any #
            P[int(coordinates[-1,2]),0] = P[int(coordinates[-1,2]),0] + 1


        # updating the state with new transition #
        # not coming from zero state and not going to zero state #
        elif transit[1][0] != 0 and transit[0][0] != 0:

            # a random number for every MV that is capable of the transition
            ru = np.random.uniform(0,1, size = (int(np.sum(coordinates[:,2] == transit[0][0]))))

            # find the max value #
            rmax = np.where(ru == np.max(ru))[0]

            # get the location of MV capable of making the transition #
            loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]

            # randomly choosing a MV to make the transition #
            coordinates[loc_mv[rmax],2] = transit[1][0]

            P[transit[0][0],0] = P[transit[0][0],0] - 1
            P[transit[1][0],0] = P[transit[1][0],0] + 1


        # handling the removal of an mv #
        elif transit[1][0] == 0:

            ru = np.random.uniform(low=0.0, high=1.0, size=(int(np.sum(coordinates[:,2] == transit[0][0]))))

            loc_mv = np.where(coordinates[:,2] == transit[0][0])[0]
            rmax = np.where(ru == np.max(ru))[0]

            # record the lifetime of the removed MV #
            mv_lifetime = coordinates[loc_mv[rmax],4]

            coordinates = np.delete(coordinates, loc_mv[rmax],0)

            P[transit[0][0],0] = P[transit[0][0],0] - 1
            P[transit[1][0],0] = P[transit[1][0],0] + 1

            # Record MV lifetimes into bins with a width of 0.25s #
            #lt_index = lt_bin_calc(mv_lifetime,mv_lifetime_index)
            #lt_bin[0,lt_index] = lt_bin[0,lt_index] + 1
        t = t + tau
        t_list=np.append(t_list,t)
        af_list=np.append(af_list,len(coordinates)*r**2)
        cf_list=np.append(cf_list,(2000-len(antigen_coordinates2))*1/2000)
        
        
        
        pop[0,n] = P[0,0]
        pop[1,n] = P[1,0]
        pop[2,n] = np.sum(P[2:22,0])
        pop[3,n] = np.sum(P[22:42,0])
        pop[4,n] = np.sum(P[42:62,0])
        pop[5,n] = np.sum(P[62:82,0])
        pop[6,n] = np.sum(P[82:102,0])
        pop[7,n] = np.sum(P[102:122,0])
        pop[8,n] = np.sum(P[122:142,0])
        pop[9,n] = np.sum(P[142:162,0])
        pop[10,n] = np.sum(P[162:182,0])
        pop[11,n] = np.sum(P[182:202,0])
        pop[12,n] = np.sum(P[202:222,0])
        pop[13,n] = np.sum(P[222:242,0])


        tau_list[0,n] = t
        if np.sum(pop[2:14,n]!=0)>0:
            ploc = np.where(pop[2:14,n]!=0)
            #print(ploc)
            pop_loc[0,n] = ploc[0][0]+1
        else:
            ploc = 0
            pop_loc[0,n] = 0
        #pop_loc[0,n] = ploc[1][0]
        n=n+1

        #fig, (ax1, ax2) = plt.subplots(1, 2)

#        if o % 1000 == 0:
#            ax1 = plt.subplot(121, polar=True)
#            ax1.grid(False)
#            ax1.set_xticklabels([])
#            ax1.set_yticklabels([])
#            ax1.set_rmax(1)
#            ax1.set_rmin(0)
#            if len(antigen_coordinates)!=0:
#                for i in range(len(antigen_coordinates)):
#                    ant  = plt.Circle((antigen_coordinates[i,0]*np.cos(antigen_coordinates[i,1]) , antigen_coordinates[i,0]*np.sin(antigen_coordinates[i,1])), r/10, transform=ax1.transData._b, color="red", alpha=1.0)
#                    ax1.add_artist(ant)
#            for i in range(len(coordinates)):
#                circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax1.transData._b, color="blue", alpha=0.6)
#
#                ax1.add_artist(circle)
#
#            #plt.title("AF ="+str(len(coordinates)*r**2) + ", AC = "+str((1000-len(antigen_coordinates))*1/1000)+", t = "+str(t))
#            #print("WTF",t_list)
#            if t>1:
#                #print(len(t_list))
#                ax2 = plt.subplot(122)
#                plt.plot(t_list,af_list)
#                plt.plot(t_list,cf_list)
#                ax2.axes.set_xlim([0,60])
#                ax2.axes.set_ylim([0,1])
#                plt.title('Cumulative and Instantaneous Area Fractions')
#                plt.ylabel('Fraction of Area Covered')
#                plt.xlabel('s')
#
#            if o<1:
#                plt.pause(5.05)
#            else:
#                plt.pause(0.05)
#                plt.clf()
#        o=o+1
        #ax2 = plt.subplot(tau_list[0,:])

    #
    #rebinding time#
#    i=0
#    while i < len(pop_loc[0]):
#        if pop_loc[0,i] == 1:
#            ti = tau_list[0,i]
#
#            while pop_loc[0,i] > 0 and pop_loc[0,i] < 12:
#                i = i + 1
#                #print(i)
#            if pop_loc[0,i] == 12 or pop_loc[0,i] == 0 and tau_list[0,i]!=0:
#                tf = tau_list[0,i]
#                #rebind_time = np.append(rebind_time,(tf - ti))
#                rebind_time = np.append(rebind_time,((tf - ti)-4.44659319327913)**2)
#
#        else:
#            i = i + 1

    #print(tau_list[0,0:10])
    #print(pop_loc[0,0:10])


        #tau_list = np.append(tau_list,tau)
        #print(tau)
        #print(P)
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
    count = MV_antigen_dist3(antigen_coordinates2, coordinates, r)
       
    return (accumulator, coordinates, antigen_coordinates,af_list,cf_list,t_list, count)


r = np.array([0.01, 0.0313727, 0.0555684])
ave_max = [5037, 512, 162]
k_s = 0.0
koff = 1.0
count_sum = 0
capture_mlist = np.zeros((3,30))
i=1
while i <= 30:
    t_sim = i
    for k in range(len(r)):
        j=0
        count_sum = 0
        while j < 20:
            (x,c,d,af,cf,t,count)=mv_sim_num_direct( 0, koff, 1, 10.0, r[k], ave_max[k], t_sim )
            #(x,c,d,af2,cf2,t2)=mv_sim_num_direct( 0, 50.5, 1,0.0 )
            count_sum = count_sum + count
            #print(count)
            j=j+1
        #print("r,t,count_ave =",(r[k],t_sim,count_sum / j))
        capture_mlist[k,t_sim] = count_sum / j
        #k=k+1
    #print(capture_mlist)
    i=i+1
print(capture_mlist)
np.save('capture_mlist.npy', capture_mlist)
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

#
#    np.save('AF_r' + str("%.3f" % round(r[i],3)) +'.npy', af)
#    np.save('CF_r' + str("%.3f" % round(r[i],3)) +'.npy', cf)
#    np.save('t_r' + str("%.3f" % round(r[i],3)) +'.npy', t)
# plt.plot(t,af,label='$y = numbers')
# plt.ylabel('Fraction of Area Covered')
# plt.xlabel('s')
# plt.show()
#print(x,d)
#start = time.time()
#conc_list = np.zeros((1,30))
#p = 0
#accumulator_sum = 0
#while p < 1:
#    (accumulator, lt_bin, pop, tau_list) = mv_sim_num_direct( 1, 1.0, 0.0000001 )
#    conc_list[0,p] = uconc
#    #print(uconc)
#    print(lt_bin)
#    #accumulator_sum = accumulator_sum + accumulator
#    #print(area_event)
#    p = p + 1
##print(accumulator_sum / 30.0)
#end = time.time()
#print(end - start)
