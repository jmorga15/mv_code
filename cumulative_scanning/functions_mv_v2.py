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
                print(len(coordinates))
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
    coordinates[loc_nodes,]
    ### generate a random vector of size l ###
    u = np.random.uniform(low=0.0, high=1.0, size=(l,))

    ###calculate the tau vec ###
    tau = 1/kx*np.log(1.0/u)

    return(np.min(tau))

def MV_antigen_dist(theta, r_antigen, coordinates, r):
    x_antigen = r_antigen*np.cos(theta)
    y_antigen = r_antigen*np.sin(theta)

    ### determine MV distance from antigen ###
    d = np.sqrt( (x_antigen - coordinates[:,0])**2 + (y_antigen - coordinates[:,1])**2)
    MV_activation = (d <= r)
    hit = np.where(MV_activation == True)[0]
    if coordinates[hit,2] < 2:
        coordinates[hit,2] = 1

    return(coordinates)
#    if True in MV_activation:
#        print(coordinates)

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

    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])

    if len(coordinates) != 0:

        for i in range(len(coordinates)):

            d = np.sqrt( (x_antigen - coordinates[i,0])**2 + (y_antigen - coordinates[i,1])**2)
            #print(d)
            MV_activation = (d <= r)

            #print (MV_activation)
            if True in MV_activation:

                if coordinates[i,2] < 2 and coordinates[i,5] == 0:
                    coordinates[i,2] = 1
                    coordinates[i,5] = np.sum(MV_activation)




    return(coordinates)
    #    if True in MV_activation:
    #        print(coordinates)

def MV_antigen_dist1(antigen_coordinates, coordinates, r):

    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])

    d = np.sqrt( (x_antigen - coordinates[-1,0])**2 + (y_antigen - coordinates[-1,1])**2)
    #print(d)
    MV_activation = (d <= r)

    #print (MV_activation)
    if True in MV_activation:
        coordinates[-1,5] = np.sum(MV_activation)
    coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])

    return coordinates


def signal_accumulator_dir(P, tau, accumulator, k_act, gam, Xact, t):
    #accumulator = np.float128(accumulator)
    #tau = np.float128(tau)
    ### not accounting for the change back ??? ###
    ### not accumulating signal while in stage 3 ###
    ### make something cool happen when t_cell is activated ###
    #print("accum before update", accumulator)
    signaling_mv = np.sum(P[7,0])
    #signaling_mv = np.sum(P[222:262,0])
    #print("sig_mv =",signaling_mv)

    ### first, add up the nodes that are already accumulating signal ###
    ### 0 indicates the node has been accumulating signal ###
#    signaling_nodes4 = np.sum((coordinates[loc4,3] == 0))
#    signaling_nodes3 = np.sum((coordinates[loc3,3] == 0))
    #signaling_nodes = np.sum((coordinates[:,3] == 0))
    #print(accumulator)
    ### for single mv ###
    ### there can be very large time jumps so must increment in small steps ###
    #print(tau)
    #print(accumulator + tau*(signaling_mv)*(k_act * (1.0-accumulator)))
    lj=0
    if tau > 0.1:
        #print("hi")
        #print(tau)
        l = 0
        while l < tau and accumulator < Xact:
        #added a decay 0.05
            #accumulator = accumulator - 0.0001*0.05*accumulator
            accumulator = accumulator + 0.001*(signaling_mv)*(k_act * (1.0-accumulator)) - 0.001*gam*accumulator
            #accumulator = accumulator + 0.0001*(signaling_mv)*(k_act) - 0.0001*gam*accumulator*(signaling_mv+1)
            l = l + 0.001
            if accumulator > Xact:
                long_shot = 1
                #print("t,l = ",(t,l))
                t = t + l
                lj=1
                #print("long jump")
                #break
    else:
        #accumulator = accumulator*(1.0 - tau*0.05)
        accumulator = accumulator + tau*(signaling_mv)*(k_act * (1.0-accumulator)) - tau*gam*accumulator
        #accumulator = accumulator + tau*(signaling_mv)*(k_act) - tau*gam*accumulator*(signaling_mv+1)
    #print("accum after update", accumulator)
    #if accumulator < 0.012365:
    #    accumulato = 0.012365
    return(accumulator, t, lj)

def lt_bin_calc(mv_lifetime,mv_lifetime_index):
    lt_minus = mv_lifetime_index - mv_lifetime
    index = np.where(lt_minus == np.min(lt_minus[lt_minus > 0]))
    return (index)
