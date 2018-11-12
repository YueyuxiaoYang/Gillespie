# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 21:17:32 2018

@author: kevin1024
"""
import numpy as np
import matplotlib.pyplot as plt

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
    x_trace = np.zeros([10000,len(x_list)])
    count = 0
    while t < T:
        if count>=10000:
            break
        # compute reaction rate
        reaction_rates = compute_reaction_rate(x,k_list,V)
        if sum(reaction_rates)==0: # no reaction is ongoing
            break     
        if len([i for i in reaction_rates if i<0])>0:
            print ('negative reaction rates')            
            continue
        
        #print (reaction_rates, count)
        R_tot = sum(reaction_rates)                
        # next reaction time
        delta_t = np.random.exponential(1/R_tot)
        # choose a reaction             
        try:        
            r_i = np.random.choice(range(k_list.shape[0]),p=reaction_rates/reaction_rates.sum())
        except:
                print ('!!!wroonnngg!!',r_i,reaction_rates/reaction_rates.sum())
                continue
        x = x+V[r_i]
        if len([i for i in x if i<0])>0:
            x += V[r_i]
            print ('negative x')
            continue
        
        x_trace[count,:] += x
        count += 1
        t += delta_t        
        
    return x_trace[:count,:]

'''
# A + B →_k1 C ->_k2 D + E   
V = [[-1,-1,1],[1,1,-1]]
x_list = np.array([11,10,10])
k_list = np.array([1.0,1.0])

x_tot = Gillespie(x_list,k_list,30,V)
plot(x_tot)
'''

'''
# ------- 3 Gene expression---------
# dna_close ->(k1) dna
# dna ->(k) rna + dna
# rna ->(k) prot + rna
# rna ->(k) 
# prot -> (k)

subs = ['dna_close', 'dna', 'rna', 'prot']
V = [[-1,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,-1,0],[0,0,0,-1]]
x_list = np.array([10,0,0,0])

k1 = 5 # dna open rate
k_list = np.array([k1,1,1,1,1])
x_tot = Gillespie(x_list,k_list,50,V)

for i in range(4):
    
    line =plt.plot(x_tot[:,i])
plt.legend(['dna_close', 'dna', 'rna', 'prot'])
plt.title('k1 = '+str(k1))
'''

# ----------prey and predator----------
# -> prey
# prey+ pd -> 2 pd
# pd ->

# in this case, we can not use V directly, since r2 is py+pd -> 2pd
subs = ['prey','pd']
V = [[1,0],[-1,1],[0,-1]]
x_list = np.array([0,0])
k1 = 1
k2 = 0.1
k3 = 1.5
k_list = np.array([k1,k2,k3])
x_tot  = Gillespie(x_list,k_list,100000,V)
for i in range(2):
    line =plt.plot(x_tot[:,i])
plt.legend(['prey','pd'])
plt.title('k = '+str(k))

# compare to ODE, eular algo
# dx/dt = k1-k2x*y, dy/dt = k2*x*y - k3y
# x_(t+delta(t)) = delta(t)*(k1-k2*x*y)+x_(t)

n = 100000 # number of iters
x = np.zeros(n)
y = np.zeros(n)
k1 = 1
k2 = 0.1
k3 = 1.5

x[0] = 10
y[0] = 10
dt = 0.01

for i in range(n-1):
    x[i+1] = dt*(k1-k2*x[i]*y[i]) + x[i] 
    y[i+1] = dt*(k2*x[i]*y[i]-k3*y[i]) + y[i]

plt.plot(x,label='prey')
plt.plot(y,label='predator')
#plt.title('k1=',k2,k3 = '+str(k))





















