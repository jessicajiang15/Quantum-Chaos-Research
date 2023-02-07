from pyemd import emd
import numpy as np
import pandas as pd
import time

mat0=np.genfromtxt("mat0.csv", delimiter=',')
supply=np.abs(np.genfromtxt("supplyamount.csv", delimiter=','))
demand=np.abs(np.genfromtxt("demandamount.csv", delimiter=','))
demandzeros=np.zeros(np.shape(demand));
supplyzeros=np.zeros(np.shape(supply));
thesupply=np.append(supply,demandzeros)
thedemand=np.append(supplyzeros,demand)
first_histogram = thesupply.copy(order='C')
second_histogram = thedemand.copy(order='C')
distance_matrix = mat0.copy(order='C')
start = time.time()
themed=emd(first_histogram, second_histogram, distance_matrix)
end =time.time()
print(end-start)
print(themed)
