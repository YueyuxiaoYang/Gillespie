# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 21:17:32 2018

@author: kevin1024
"""
import numpy as np

def compute_reaction_rate(x_list,k_list,V):
    '''
        a mass function: lambda(r_i) = x
        V: column->substrate, row->reaction, ("state change matrix")^T

    '''  
    reaction_rates = []
    for r_i in range(len(k_list)):
        sub_ind = [i for i,v in enumerate(V[r_i]) if v<0]
        r_rate = k_list[r_i]
        for n in x_list[sub_ind]:
            r_rate = r_rate*n
        reaction_rates.append(r_rate)
    return np.array(reaction_rates)

def Gillespie(x_list, k_list, T, V):
    '''
        x_list: substrate vector
        k_list: reaction rates
        T: stop time
        V: state change matrix
        
        func: compute_reaction_rate
    '''
    x = x_list
    t = 0
    x_trace = np.zeros([1000,len(x_list)])
    count = 0
    while t < T:
        if count>=500:
            break
        # compute reaction rate
        reaction_rates = compute_reaction_rate(x,k_list,V)
        if sum(reaction_rates)==0: # no reaction is ongoing
            break     
        if len([i for i in reaction_rates if i<0])>0:
            print 'negative reaction rates'            
            continue
        
        print reaction_rates, count
        R_tot = sum(reaction_rates)                
        # next reaction time
        delta_t = np.random.exponential(1/R_tot)
        # choose a reaction             
        try:        
            r_i = np.random.choice(range(k_list.shape[0]),p=reaction_rates/reaction_rates.sum())
        except:
                print '!!!wroonnngg!!',r_i,reaction_rates/reaction_rates.sum()
                continue
        x = x+V[r_i]
        if len([i for i in x if i<0])>0:
            x += V[r_i]
            print 'negative x'
            continue
        
        x_trace[count,:] += x
        count += 1
        t += delta_t        
        
    return x_trace[:count,:]

# A + B →_k1 C ->_k2 D + E   
V = [[-1,-1,1],[1,1,-1]]
x_list = np.array([11,10,10])
k_list = np.array([1.0,1.0])

x_tot = Gillespie(x_list,k_list,30,V)
plot(x_tot)

