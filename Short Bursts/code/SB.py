import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  
from scipy import stats
import csv
from math import pi
from astropy.coordinates import Angle
from astropy import units as u
from scipy.stats import poisson
from scipy.misc import factorial
from scipy.optimize import curve_fit

data = pd.read_csv('values.csv')             #read csv file
data= data[data['t90']<=2]                    #constraint
data= data[data['t90']>0.1]

angles1 = np.array(data['lii'])
angles2 = np.array(data['bii'])  

longitude =  np.radians(angles1)                   
latitude =   np.radians(angles2)                 


l = len(longitude)

for i in range(l):
	if(longitude[i]>=pi):
		longitude[i]=longitude[i]-2*pi


length = len(angles1)

for i in range(length):
	if(angles1[i]>=180):
		angles1[i]=angles1[i]-360

region = np.array([0]*8)

# Function to check quadrant 
def quadrant(x,y): 
    if (x <= -90 and -30<=y<=30): 
        region[0] = region[0]+1 
  
    elif (-90<=x<=0 and -30<=y<=30): 
        region[1] = region[1]+1 
          
    elif ( 0<=x<=90 and y <= 0): 
        region[2] = region[2]+1 
      
    elif (x >= 90 and -30<=y<=30): 
        region[3] = region[3]+1 
          
    elif (-180<=x<=0 and y>=30): 
        region[4] = region[4]+1 
      
    elif (0<=x<=180 and y>=30): 
        region[5] = region[5]+1 
      
    elif (-180<=x<=0 and y<=-30): 
        region[6] = region[6]+1 
      
    elif (0<=x<=180 and y<=-30): 
        region[7] = region[7]+1 
      
for f, b in zip(angles1,angles2):
    quadrant(f,b)


def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)

mu = length/8
x = np.linspace(0,80, 1000)

plt.figure(1)
plt.hist(region,bins=25)
plt.plot(x, poisson(x,mu)*10, 'r--', label='Poisson fit')
plt.xlim(0,80)
plt.xlabel('Number of GRBs bursts inside a region')
plt.ylabel('Number of regions')
plt.title('448 Gamma Bursts Poison distribution')


plt.figure(2,figsize=(8,4.2))
plt.subplot(111,projection='aitoff')
plt.grid(True)
plt.title('BATSE SB (0.1s<T90<=2s) GRB EVENTS')
plt.xlabel('Galactic Coordinates - 8 regions')
plt.scatter(longitude,latitude,color='black')
plt.subplots_adjust(top=0.95,bottom=0.0)
plt.show()
