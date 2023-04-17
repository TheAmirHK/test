# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %TheAmirHK
"""

import numpy as np
import dimod
from dimod import Binary, Integer, ConstrainedQuadraticModel, SimulatedAnnealingSampler

import numpy as np
import urllib.request

url = "https://github.com/TheAmirHK/test/raw/main/matrix_form.npy"
filename = "matrix_form.npy"

# Download the file
urllib.request.urlretrieve(url, filename)

# Load the file using np.load()
matrix_form = np.load(filename)
# In[] 
import numpy as np
import dimod
from dimod import BinaryQuadraticModel

def optimization_code(matrix_form):
    admissble_KTE = 23
    rang = 1000
    x = {(i,j): 'x({},{})'.format(i,j) for i in range(rang) for j in range(rang)}
    y = {(i,j): 'y({},{})'.format(i,j) for i in range(rang) for j in range(rang)}
    bqm = BinaryQuadraticModel.empty(dimod.BINARY)
    
    for i in range(rang):
        for j in range(i+1,rang):
            bqm.add_variable(x[(i,j)], matrix_form[i,j])
            bqm.add_variable(y[(i,j)], 1000)
            bqm.add_interaction(x[(i,j)], y[(i,j)], 1000)

    # Constraints
    for i in range(rang):
        x_sum = dimod.quicksum(x[(i,j)] for j in range(rang))
        y_sum = dimod.quicksum(x[(i,j)] for j in range(rang))
        bqm.add_linear_equality_constraint([(x_sum, 1), (y_sum, 1)], admissble_KTE, 10)
        
        x_sum = dimod.quicksum(x[(j,i)] for j in range(rang))
        y_sum = dimod.quicksum(y[(j,i)] for j in range(rang))
        bqm.add_linear_equality_constraint([(x_sum, 1), (y_sum, 1)], admissble_KTE, 10)
        
        bqm.add_linear_equality_constraint([(x[(i,j)], 1), (y[(i,j)], 1)], 1, 10)

    # Objective function
    response = dimod.SimulatedAnnealingSampler().sample(bqm)
    energies = response.data_vectors['energy']
    indices = np.argwhere(energies == energies.min())[0]
    x_min = {k: response.record.sample[0][i] for i, k in enumerate(bqm.variables)}
    y_min = {k: response.record.sample[0][i] for i, k in enumerate(bqm.variables)}
    m = 0
    n = 0
    kte_list = []
    spur_ = []
    crown_ = []
    for i in range (rang):
        for j in range(i+1, rang):
            if x_min[x[(i,j)]] == 1:
                print (i, j, matrix_form[i,j])
                m += 1
                kte_list.append(matrix_form[i,j])
            if y_min[y[(i,j)]] == 1:
                n += 1
                spur_.append(i)
                crown_.append(j)

    print ("Residual spur gears", np.asarray(spur_)+1001)
    print("Residual crown wheels",np.asarray(crown_)+1001)
    print ("Number of residuals %0.f"%n)
    print("number of pairs %0.f"%m)
    print ("Mean KTE value = %0.2f"%((energies.min() - n*1000)/m))
    return kte_list


run = optimization_code(matrix_form)
