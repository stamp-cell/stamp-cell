# Open the postprcoessing.txt. Next read the J_A from each saved ones. After that, add all the J_A s together.
# Do the same for other species. 
# If needed, start thinking about adding key to convert A to the specific name of the speicies.

import numpy as np
from pandas import Series, DataFrame
import pandas as pd
import json
from datetime import datetime
from dateutil.parser import parse # dateutil is capable of parsing almost any human-intelligible date representation
import pytz
import seaborn
import matplotlib.pyplot as plt
import scipy.stats 

'''
def file_read(fname):
        content_array = []
        with open(fname) as f:
                #Content_list is the list that contains the read lines.     
                for line in f:
                        content_array.append(line[:-1])
                print(content_array)

file_read('postprocessing.txt')
#print (content_array)
postreading = open("postprocessing.txt","r+") 
print postreading.readlines()
print "Output of Read function is "
#print postreading.read() 

# define an empty list
places = []

# open file and read the content in a list
with open('postprocessing.txt', 'r') as filehandle:
    for line in filehandle:
        # remove linebreak which is the last character of the string
        currentPlace = line[:-1]

        # add item to the list
        places.append(currentPlace)
print currentPlace

# define empty list
places = []
'''
# open file and read the content in a list
with open('postprocessing.txt', 'r') as filehandle:
    places = [current_place.rstrip() for current_place in filehandle.readlines()]

print (places[0])

#f1 = open(places[0],'r') 

#____________________________________________________________
# Read the file with absolute humidity values: 
dff1 = pd.read_csv(places[0])
MColumns = ['J_AMN', 'J_BMN','J_CMN']
print ('the minimum of A flux is : ')
print (min(dff1['J_AMN']))
dff2= pd.read_csv(places[1])
MColumns2 = ['J_AMN', 'J_BMN']
print (dff1['J_AMN']+dff2['J_AMN'])





