#!/usr/bin/python
import numpy as np
import random
import time
from functions_mv_v2 import *
import sys
from mv_numerical_kp2 import*
import multiprocessing
import os
from mv_multi_sim_directg import *
from multiprocessing import Process



#print(os.cpu_count())

#import Pool
from multiprocessing import Pool

ks = 1.11
procs = []
while ks < 1.12:
  k_ind = int(30)
  while k_ind <= 39:
    #Assuming you want to use 11 processors
    if __name__ ==  '__main__':
        ### num_processors = 2
        numbers = [(0 + k_ind,ks), (1 + k_ind,ks), (2 + k_ind,ks), (3 + k_ind,ks), (4 + k_ind,ks), (5 + k_ind,ks), (6 + k_ind,ks), (7 + k_ind,ks), (8 + k_ind,ks), (9 + k_ind,ks)]
        for index, number in enumerate(numbers):
            proc = Process(target=fn_of_koff, args=(number,))
            procs.append(proc)
            proc.start()

        for job in procs:
            job.join()
    k_ind = k_ind + int(10)
  ks = ks + 0.01            
    #    for proc in procs:
    #        proc.join()
            
    #    #Create a pool of processors
    #    p=Pool(processes = num_processors)
    #    #get them to work in parallel
    #    output = p.map(fn_of_koff,[i for i in range(0,2)])
    #    print(output)
    #    #p.close()
