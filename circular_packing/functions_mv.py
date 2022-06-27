# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 15:19:39 2019

@author: Jonathan
"""
import random
import numpy as np

def gen_coords(coordinates, R, r, attempts):
    if len(coordinates) == 0:
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
        
        skip = 0

        counter = 0
        
                
        prob_count = 0
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
        while ((False in d) == True or d_org == False) and counter < attempts and skip == 0:
            x_new = random.uniform(-(1 - r), (1 - r))
            y_new = random.uniform(-(1 - r), (1 - r))
            d = np.sqrt( (x_new - coordinates[:,0])**2 + (y_new - coordinates[:,1])**2)
            d = d >= 2*r

            d_org =  np.sqrt( (x_new - 0)**2 + (y_new - 0)**2)
            d_org = d_org <= (1 - r)

            counter = counter + 1
            if counter >= attempts:
                skip = 1
                #print("hi")
                
        prob_count = counter + 1

    return(x_new,y_new, prob_count, skip)

    
    
