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




def MV_antigen_dist2(antigen_coordinates,antigen_coordinates2, coordinates, r):

    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])

    x_antigen2 = antigen_coordinates2[:,0] * np.cos(antigen_coordinates2[:,1])
    y_antigen2 = antigen_coordinates2[:,0] * np.sin(antigen_coordinates2[:,1])

    if len(coordinates) != 0:

        for i in range(len(coordinates)):

            x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
            y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])

            x_antigen2 = antigen_coordinates2[:,0] * np.cos(antigen_coordinates2[:,1])
            y_antigen2 = antigen_coordinates2[:,0] * np.sin(antigen_coordinates2[:,1])

            d = np.sqrt( (x_antigen - coordinates[i,0])**2 + (y_antigen - coordinates[i,1])**2)

            d2 = np.sqrt( (x_antigen2 - coordinates[i,0])**2 + (y_antigen2 - coordinates[i,1])**2)
            #print(d)
            MV_activation = (d <= r)

            MV_activation2 = (d2 <= r)

            #print (MV_activation)
            if True in MV_activation:
                if coordinates[i,2] < 2 and coordinates[i,5] == 0:
                    coordinates[i,2] = 1
                    coordinates[i,5] = np.sum(MV_activation)
            if True in MV_activation2:

                loc = np.where(MV_activation2==True)
                #print(loc[0])
                #print(len(loc[0]))
                # print(len(antigen_coordinates))
                # print(antigen_coordinates)
                #numpy.delete(antigen_coordinates, (loc[0][0]), axis=0)
                j=len(loc[0])-1
                while j >= 0:
                    antigen_coordinates2=np.delete(antigen_coordinates2, (loc[0][j]), axis=0)
                    j = j - 1


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

    return(coordinates,antigen_coordinates,antigen_coordinates2)
    #    if True in MV_activation:
    #        print(coordinates)
    
def MV_antigen_dist3(antigen_coordinates2, coordinates, r):
    ### determining how many spacing points have a full MV radius open around them ###
    ### This function uses vogels method around every undiscovered point to determine an approximation #
    ### of the space available for MV placement in which the mv can capture the point ###
    
    # the count will be the total number of Vogel's points where a MV can be placed to capture the undiscovered location #
    count = 0


    if len(coordinates) != 0:            
        
        # antigen_coordinates2 is a list of undiscovered locations in the IS, spacing by Vogel's Method #
        for i in range(len(antigen_coordinates2)):
            
            ### distribute the possible capture points within 1 MV radius of the point ###
            
            # n is the number of capture points to use. More of these points yield higher accuracy, more computationally expensive. #
            n = 200
            antigen_coordinates3 = np.zeros((n,2))
            j = 0
            golden_angle = np.pi * (3 - np.sqrt(5))
            
            # The while loop spaces the points using Vogel's Method #
            while j < n:
                    
                    theta = j * golden_angle
                    r_antigen = antigen_coordinates2[i,0] + r*(np.sqrt(j) / np.sqrt(n))
        
        
                    #r_antigen,theta = antigen_antigen_dist(antigen_coordinates, r)
        
                    #antigen_coordinates[i,0] = r_antigen
        
                    #antigen_coordinates[i,1] = theta
        
                    antigen_coordinates3[j,0] = r_antigen
        
                    antigen_coordinates3[j,1] = theta
        
                    j = j + 1

            
            for k in range(len(antigen_coordinates3)):
                ### for each possible capture point created above, this section checks whether or not a MV can be placed there. ###
                ### the process for this check is to determine if any other MV are within 2*r of this point, i.e. a new MV would be centered here ###
                x_antigen2 = antigen_coordinates3[k,0] * np.cos(antigen_coordinates3[k,1])
                y_antigen2 = antigen_coordinates3[k,0] * np.sin(antigen_coordinates3[k,1])
    
    
                # distance check with all MV coordinates. #
                d2 = np.sqrt( (coordinates[:,0] - x_antigen2)**2 + (coordinates[:,1] - y_antigen2)**2)
                
                 
                # its 2*r because each MV is tracked by its center. Assuming that a MV is placed with its new center at this new capture point #
                # then 1 MV radius reaches this MV's boundary, and any MV center within 1 MV radius of this boundary would be overlapping #
                MV_activation2 = (d2 <= 2*r)
    
                if True in MV_activation2:
                    ### means there is at least one mv within r distance of point ###
                    count = count
                # Here, no MV overlap, however, still need to check if the capture point overlaps the IS boundary. #
                elif r_antigen < 1.0 and (1 - r_antigen) >= r:
                    count = count + 1

        #print("r's = ",(antigen_coordinates2[i,0],r_antigen))


    return(count/len(antigen_coordinates2))
    
    
    
def MV_antigen_dist1(antigen_coordinates,antigen_coordinates2, coordinates, r):

    x_antigen = antigen_coordinates[:,0] * np.cos(antigen_coordinates[:,1])
    y_antigen = antigen_coordinates[:,0] * np.sin(antigen_coordinates[:,1])

    x_antigen2 = antigen_coordinates2[:,0] * np.cos(antigen_coordinates2[:,1])
    y_antigen2 = antigen_coordinates2[:,0] * np.sin(antigen_coordinates2[:,1])

    d = np.sqrt( (x_antigen - coordinates[-1,0])**2 + (y_antigen - coordinates[-1,1])**2)

    d2 = np.sqrt( (x_antigen2 - coordinates[-1,0])**2 + (y_antigen2 - coordinates[-1,1])**2)
    #print(d)
    MV_activation = (d <= r)

    MV_activation2 = (d2 <= r)


    #print (MV_activation)
    if True in MV_activation2:

        loc = np.where(MV_activation2==True)
        #print(loc[0])
        # print(len(antigen_coordinates))
        # print(antigen_coordinates)
        #numpy.delete(antigen_coordinates, (loc[0][0]), axis=0)
        i=len(loc[0])-1
        while i >= 0:
            antigen_coordinates2=np.delete(antigen_coordinates2, (loc[0][i]), axis=0)
            i = i - 1
        #antigen_coordinates=np.delete(antigen_coordinates, loc[0][0])
        # print(antigen_coordinates)
        # print(len(antigen_coordinates))
        # print(loc)
    if True in MV_activation:
        coordinates[-1,5] = np.sum(MV_activation)
    coordinates[-1,2] = coordinates[-1,5] + 1
    #print("state = ", coordinates[-1,2])

    return (coordinates,antigen_coordinates,antigen_coordinates2)

def signal_accumulator(coordinates, tau, accumulator, k_act, t_signal):
    ### not accounting for the change back ??? ###
    ### not accumulating signal while in stage 3 ###
    ### make something cool happen when t_cell is activated ###

    loc4 = np.where((coordinates[:,2] == 4) == True)[0]
    loc3 = np.where((coordinates[:,2] == 3) == True)[0]

    ### first, add up the nodes that are already accumulating signal ###
    ### 0 indicates the node has been accumulating signal ###
#    signaling_nodes4 = np.sum((coordinates[loc4,3] == 0))
#    signaling_nodes3 = np.sum((coordinates[loc3,3] == 0))
    signaling_nodes = np.sum((coordinates[:,3] == 0))
    #print(accumulator)
    ### for single mv ###
    ### there can be very large time jumps so must increment in small steps ###
#    if t_signal > 0:
#        print("hi")
#        l = 0
#        while l < t_signal:
#        #print("t_signal =",t_signal)
#            accumulator = accumulator + 0.01*(signaling_nodes)*(k_act * (1-accumulator))
#            l = l + 0.01
#    else:
#        accumulator = accumulator + tau*(signaling_nodes)*(k_act * (1-accumulator))
    #accumulator = accumulator + tau*(signaling_nodes3 + signaling_nodes4)*(k_act * (1-accumulator))

    ### might be a better approach, take into account nodes that were just shut off ###
    accumulator = accumulator + tau*(signaling_nodes)*(k_act * (1-accumulator))
    #print(tau)
    #print(accumulator)
    #print("3 and 4", signaling_nodes3,signaling_nodes4)
    #print("tau = ", tau)
    #print("accumulator =",accumulator)
    ### second, turn on new signaling nodes ###
    ### -1 indicates a node was previously not accumulating signal ###
    if len(loc4 != 0):
        for i in loc4:
            if coordinates[i,3] == -1:
                coordinates[i,3] = 0
    if len(loc3 != 0):
        for i in loc3:
            if coordinates[i,3] == -1:
                coordinates[i,3] = 0

#    non_signaling_nodes4 = np.where((coordinates[loc4,3] == -1) == True)[0]
#    print("non_signal_nodes4", non_signaling_nodes4)
#    coordinates[non_signaling_nodes4,3] = 0

    #non_signaling_nodes3 = np.where((coordinates[loc3,3] == -1))
    #coordinates[non_signaling_nodes3,3] = 0
    #print(coordinates)

    for i in range(len(coordinates)):
        if coordinates[i,2] != 3 and coordinates[i,2] != 4 and coordinates[i,3] == 0:
            coordinates[i,3] = -1

    return(coordinates, accumulator)

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
def signal_accumulator_dir(P, tau, accumulator, k_act, gam):
    #accumulator = np.float128(accumulator)
    #tau = np.float128(tau)
    ### not accounting for the change back ??? ###
    ### not accumulating signal while in stage 3 ###
    ### make something cool happen when t_cell is activated ###
    #print("accum before update", accumulator)
    signaling_mv = np.sum(P[222:242,0])
    #signaling_mv = np.sum(P[7,0])
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
    if tau > 0.001:
        #print("hi")
        #print(tau)
        l = 0
        while l < tau:
        #added a decay 0.05
            #accumulator = accumulator - 0.0001*0.05*accumulator
            accumulator = accumulator + 0.0001*(signaling_mv)*(k_act * (1.0-accumulator)) - 0.0001*gam*accumulator
            #accumulator = accumulator + 0.0001*(signaling_mv)*(k_act) - 0.0001*gam*accumulator*(signaling_mv+1)
            l = l + 0.0001
    else:
        #accumulator = accumulator*(1.0 - tau*0.05)
        accumulator = accumulator + tau*(signaling_mv)*(k_act * (1.0-accumulator)) - tau*gam*accumulator
        #accumulator = accumulator + tau*(signaling_mv)*(k_act) - tau*gam*accumulator*(signaling_mv+1)
    #print("accum after update", accumulator)
    return(accumulator)
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
##    if len(loc3 != 0):13
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
