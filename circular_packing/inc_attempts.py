# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 14:15:27 2019

@author: Jonathan
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:07:53 2019

@author: Jonathan
"""

import numpy as np
#import matplotlib.pyplot as plt
#import pylab as pl
import random
#import time

#from IPython import display
from functions_mv import *

def fill_circle(t_sim):
    radii = np.array([0.01, 0.0313727,0.0555684])
    attempts=np.logspace(2.0, 6.0, num=20)
    ### number of sims per parameter set ###
    num_sim = 100
    ### Radius of Cell Interface ###
    R = 1
    
    ### Radius of MV ###
    r = radii[0]
    
    ### Rate of MV Arrival ###
    ka = 1.0
    kr = 0.000000001
    #array_of_counters = np.array([1.0])
    ### List of center points for MV ###
    #coordinates = np.zeros
    #coordinates = np.append(coordinates,[[1,1]], axis = 0)
    
    ### Position the first circle ###

    #r = 0.2 + np.random.randn(N)*0.1
    
    ### average number of circles per radius r ###
    ave_circle_array = np.zeros((3,20))
    master_array_of_counters = np.zeros((3, 20, num_sim, 6000))
    ### Loop through the three radii ###
    for i in range(len(radii)):
        r = radii[i]
        
        ### Loop through the different attempt numbers
        for j in range(len(attempts)):
            
            ### store the 100 simulations. Each row is a single sim. 
            #master_array_of_counters = np.zeros((10,6000))
            master_circle_sim = 0
            
           
            # number of MV #
            circle_count = 0
            
            ### How many simulations with the same parameters ###
            while master_circle_sim < num_sim:
                array_of_counters = np.zeros((1,6000))
    
                coordinates = np.array([])
                (x_new, y_new, prob_count, skip) = gen_coords(coordinates, R, r,attempts[j])
            
                coordinates = np.array([[x_new,y_new]])
            
                t = 0
                while t < 6000:
                    
                    ### Draw a random number ###
                    u = random.uniform(0, 1)
                    

                    
                    (x_new, y_new, prob_count, skip) = gen_coords(coordinates, R, r, attempts[j])
                    
                    array_of_counters[0,t] = prob_count
                    
                    if skip == 1:
                        t = 6200
                    else:
                        coordinates = np.append(coordinates,[[x_new,y_new]], axis = 0)
                        t = t + 1
    
                    #coordinates = np.append(coordinates,[[x_new,y_new]], axis = 0)
                    #print(t)
        #            
                
        #    
                circle_count = circle_count + len(coordinates)

                print(circle_count)
                master_array_of_counters[i,j,master_circle_sim] = array_of_counters[0]
                master_circle_sim = master_circle_sim + 1
                
            ### average number of MV in area per simulation ###
            ave_circle_array[i,j] = circle_count / num_sim
        
        
  
        print(ave_circle_array)
    np.save('master_array_of_counters_r' + str("%.3f" % round(r,3)) + '.npy', master_array_of_counters)
    np.save('avemax_circle_array.npy', ave_circle_array)
        #
fill_circle(100)