#! /usr/bin/python

# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
#%matplotlib inline


# In[2]:

x1=pd.read_csv('fixA1.dat',header=None)
x2=pd.read_csv('fixA2.dat',header=None)


# In[6]:

plt.figure()
n1,bins1,patches1=plt.hist(x1,100,normed=True)
#plt.axis([0,400,0,500])
plt.savefig('fixA1.png')
plt.close()


# In[7]:

plt.figure()
n2,bins2,patches2=plt.hist(x2,100,normed=True)
#plt.axis([0,400,0,500])
plt.savefig('fixA2.png')
plt.close()


# In[5]:

print("fixA1 {}".format(x1.describe()))
print("fixA2 {}".format(x2.describe()))


# In[ ]:




# In[ ]:



