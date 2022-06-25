#!/usr/bin/python
import numpy as np
import random

from functions_mv_v2 import *

from mv_numerical_kp2 import*

def fn_of_koff(arg):
    koff_index = arg[0]
    ks = arg[1]
    tc_list = np.zeros((20,100))
    tc_sim = np.zeros((1,20))
    #area_covered_sim = np.zeros((1,20))
    accumulator_sim = np.zeros((1,20))
    num_antigens = np.zeros((1,20))
    #stab_conc_mean_sim = np.zeros((1,20))
    indexing_array = np.logspace(-2.0, 1.0, num=40)
    antigen_array = np.array([1,2,4,6,11,18,29,48,78,127,206,335,545,885,1438,2335,3792,6158,10000])
    lt_bin_sim = np.zeros((20,340))
    k_off = indexing_array[koff_index]
    X = np.zeros((20,100))
    Xmax = np.zeros((20,100))
    var = np.zeros((1,20))
    
#    i = 10
    m = 0
    for i in antigen_array:
    
        area_iteration = 0
        accumulator_iteration = 0
        tc_iteration = 0
        j = 1
        nkp=15
        ### number of iterations per antigen number
        while j < 101:
            #N=55#
            (accumulator, tc, X_max) = mv_sim_num_direct( i, k_off, ks, 15 )
            
            #stab_conc_mean_sim[0,m] = stab_conc_mean_sim[0,m] +  uconc
            
            #lt_bin_sim[m] = lt_bin_sim[m] + lt_bin[0]
            
            #area_iteration = area_iteration + area_event
            Xmax[m,j-1] = X_max
            X[m,j-1] = accumulator
            tc_list[m,j-1] = tc
            tc_iteration =  tc_iteration + tc
            accumulator_iteration = accumulator_iteration + accumulator
            
            j = j+1
            #print(j)
    
        #area_iteration = area_iteration / 100.0
        accumulator_iteration = accumulator_iteration / 100.0
        tc_iteration =  tc_iteration / 100.0
        #stab_conc_mean_sim[0,m] = stab_conc_mean_sim[0,m] / 100.0
    
        #area_covered_sim[0][m] = area_iteration
        accumulator_sim[0][m] = accumulator_iteration
        tc_sim[0][m] = tc_iteration
        num_antigens[0][m] = i
        var[0,m] = np.sum((X[m] - accumulator_iteration)**2)/ 99.0
        m = m + 1
        
        #print("i, koff =", (i,k_off))

    #print("Done is", k_off)
    
    
    
    #print(np.max(mv_lifetime_sim[m-1]))
    #print(accumulator_sim)
    
   ### Write data to HDF5


    #np.save('stab_conc_mean_sim' + str(koff_index) + 'ks' + str("%.2f" % round(ks,2)) + '.npy', stab_conc_mean_sim)
    np.save('signal_data'+ str(koff_index) + '.npy', accumulator_sim)
    np.save('time_data_ave_'+ str(koff_index) + 'ks' + str("%.2f" % round(ks,2)) + '.npy', tc_sim)
    #np.save('area_covered_sim'+ str(koff_index) + 'kt' + str(kt) + '.npy', area_covered_sim)
    #np.save('lt_bin_sim'+ str(koff_index) + 'kt' + str(kt) + '.npy', lt_bin_sim)
    #np.save('var'+ str(koff_index) + 'kt' + str(kt) +  '.npy', var)
    np.save('X'+ str(koff_index) + 'ks' + str("%.2f" % round(ks,2)) +  '.npy', X)
    np.save('Xmax'+ str(koff_index) + 'ks' + str("%.2f" % round(ks,2)) +  '.npy', Xmax)
    np.save('tc_i'+ str(koff_index) + 'ks' + str("%.2f" % round(ks,2)) +  '.npy', tc_list)
    


#    plt.plot(antigen_array,X[:,0])
#    plt.show()
#   
    
    
    
    return (accumulator_sim)


#fn_of_koff((18,1.0))