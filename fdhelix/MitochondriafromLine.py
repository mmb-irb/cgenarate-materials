#!/usr/bin/env python
# coding: utf-8

# In[1]:


inputfile = 'line_mitochondria.pdb'
outputfile = 'mitochondriaCG.pdb'


# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


ylist=[]
with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        ylist.append(float(line[46-8:54-8]))


# In[4]:


zlist=[]
with open(inputfile, 'r')  as fin:
    for index, line in enumerate(fin):
        zlist.append(float(line[46:55]))


# In[5]:


dist=max(zlist)-min(zlist)
dist


# In[6]:


y=np.array(ylist)
z=np.array(zlist)


# In[7]:


r=y+dist/(2*np.pi)
theta= z/dist*(2*np.pi)
r,theta


# In[8]:


y2, z2 = r*np.cos(theta), r*np.sin(theta)
y2, z2


# In[9]:


plt.plot(z2)
plt.plot(y2)


# In[10]:


a=0
b=a+10
d=16568
c=d-10
off1=-7000
off2=-39500
plt.scatter(z2[a:b],y2[a:b])
plt.scatter(z2[c:d],y2[c:d])
# plt.plot(y[a:b]+off1,z[a:b]+off2)
plt.axis('equal')


# In[11]:


with open(inputfile, 'r')  as fin:
    with open(outputfile, 'w')  as fout:
        for index, line in enumerate(fin):
#             print(line.strip())
            
            y=float(line[38:46])
            z=float(line[46:54])
            
            r=y+dist/(2*np.pi)
            theta= z/dist*2*np.pi
            y2, z2 = r*np.cos(theta), r*np.sin(theta)
            
            line2 = line[:4] + "  {:5d}".format(index+1) + line[11:38]
            line2 += "{:8.2f}".format(y2) + "{:8.2f}".format(z2) + line[54:]
            
            print(line2.strip())
            fout.writelines(line2)


# In[ ]:




