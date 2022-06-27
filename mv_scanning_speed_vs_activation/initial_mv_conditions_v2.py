#!/usr/bin/python

import numpy as np
import random
from functions_mv_v2 import *
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
    r = r
    count = 0
    ### Rate of MV Arrival ###
#    kd = 0.0001
#    ka = 1.0
    kd = 1/7.0
    ka = kd * 4.666667
    #k_off = 0.75
    k_on = 0
    k_10 = 0
    k_t = 0
    #k_21 = 0.9
    k_s = 0
    k_20 = 0
    #k_32 = 0.5
    k_34 = 0
    #k_43 = 0.5
    k_45 = 0
    k_u = 0
    k_54 = 0
    k_53 = 0
    
    #area_fraction = np.zeros((1,1000))
    ### The average number of circles fitted at given radius ###
    #ave_full = 170
    
    
    
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
    
            (coordinates, tau, count, mv_lifetime, t_signal) = tau_calc(t, t_max, count, mv_lifetime, coordinates, ka, kd, k_on, k_10, k_t, k_s, k_34, k_45, k_54, k_53, k_u, k_20, R, r, ave_full)
    
    
        ###else: there are no nodes so just add one ###
        else:
            #print(len(coordinates))
            (x_new,y_new) = gen_coords(coordinates, R, r)
    
            ### Update the time for the first MV appearance ###
            tau = gen_tau(ka, ave_full)
    
            coordinates = np.append(coordinates,[[x_new,y_new, 0, -1, 0, 0]], axis = 0)
    
        t = t + tau
        
        
        coordinates[:,4] = coordinates[:,4] + tau
        
        #ave_mv_lifetime[0,count - 1] = np.sum(mv_lifetime) / count
        
        #print(ave_mv_lifetime[0,count - 1])
        #print(np.max(mv_lifetime))
    
#        ax = plt.subplot(111, polar=True)
#        ax.grid(False)
#        ax.set_xticklabels([])
#        ax.set_yticklabels([])
#        ax.set_rmax(1)
#        ax.set_rmin(0)
#        for i in range(len(coordinates)):
#    
#         if coordinates[i,2] == 0:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="green", alpha=0.4)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 2:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="red", alpha=0.4)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 1:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="red", alpha=0.2)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 3:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="purple", alpha=0.6)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 4:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="blue", alpha=0.6)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 5:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="orange", alpha=0.4)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 6:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="grey", alpha=0.2)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 7:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="purple", alpha=0.6)
#             ax.add_artist(circle)
#    
#         elif coordinates[i,2] == 4:
#             circle  = plt.Circle((coordinates[i,0], coordinates[i,1]), r, transform=ax.transData._b, color="blue", alpha=0.6)
#             ax.add_artist(circle)
#    
#    
#        if count == 0:
#            plt.pause(0.1)
#        else:
#            plt.pause(0.1)
#        plt.clf()
#        count = count + 1
        
        #print(r**2*len(coordinates))
        
    return (coordinates)
#coordinates = mv0(0.1, 45)
#print(coordinates)