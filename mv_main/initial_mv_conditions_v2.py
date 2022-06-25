#!/usr/bin/python

import numpy as np
import random
from functions_mv2_v2 import *
import sys

def mv0(r, ave_full):
    ### initialize coordinates array ###
    coordinates = np.array([])
    
    ### an array to track lifetimes and a count for indexing ###
    mv_lifetime = np.zeros((1,100000))
    count = 1
    
    ave_mv_lifetime = np.zeros((1,100000))
    
    ### Radius of Cell Interface ###
    R = 1
    
    ### Radius of MV ###
    count = 0
    ### Rate of MV Arrival ###
    kd = 1/7.0
    ka = kd * 4.666667
    
    
    
     ### Get coordinates for first MV ###
    (x_new,y_new) = gen_coords(coordinates, R, r)
    
    ### Update the time for the first MV appearance ###
    tau = gen_tau(ka, ave_full)
    
    ### Place the first MV ###
    coordinates = np.array([[x_new,y_new, 0,-1, 0, 0]])
    
    
    t_max = 100
    
    t = tau
    #print(t)
    
    while t < 20:
    
        ### There are some MV present, multiple events are possible, use tau calc function ###
        if len(coordinates) != 0:
    
            #(coordinates, tau, count, mv_lifetime, t_signal) = tau_calc(t, t_max, count, mv_lifetime, coordinates, ka, kd, k_on, k_10, k_t, k_s, k_34, k_45, k_54, k_53, k_u, k_20, R, r, ave_full)
            
            ### Update the time for the first MV appearance ###
            tau_A = gen_tau(ka, ave_full - len(coordinates))
            tau_D = gen_tau(kd, len(coordinates))
            if tau_A < tau_D:
                # MV addition #
                (x_new,y_new) = gen_coords(coordinates, R, r)
                coordinates = np.append(coordinates,[[x_new,y_new, 0, -1, 0, 0]], axis = 0)
            else:
                # MV Removal #
                ru = random.randint(0, len(coordinates) - 1)
                coordinates = np.delete(coordinates, ru,0)
        ###else: there are no nodes so just add one ###
        else:
            #print(len(coordinates))
            (x_new,y_new) = gen_coords(coordinates, R, r)
    
            ### Update the time for the first MV appearance ###
            tau = gen_tau(ka, ave_full)
    
            coordinates = np.append(coordinates,[[x_new,y_new, 0, -1, 0, 0]], axis = 0)
    
        t = t + tau
        
        coordinates[:,4] = coordinates[:,4] + tau

        
    return (coordinates)
# r = 0.05345224838248488

# mv0(r, 60)

