# coding: utf-8

# In[ ]:

#get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np

import ROOT
from ROOT import larlite

import sys


# In[ ]:

mgr = larlite.storage_manager()
mgr.set_io_mode(mgr.kREAD)

infile = sys.argv[1]
#infile = "output_simpleflash.root"

flashproducer = str(sys.argv[2])
#flashproducer = "SimpleFlashFinder"
mgr.add_in_filename(infile)

mgr.open()

b = []
print mgr    
while mgr.next_event() :

    a = mgr.get_data(larlite.data.kOpFlash,flashproducer)

    #z = mgr.get_data(larlite.data.kOpHit,"OpHitFinder")

    if a.size() <= 0: # no flash ok fine
        continue
        
    #print a.size()
    #print z.size()
    
    for aa in xrange(a.size()):
        b.append(a[aa].TotalPE())

b = np.array(b)


# In[ ]:

fix, ax = plt.subplots(figsize=(10,6))
plt.hist(b,bins=np.linspace(0,2000,500),histtype='step')
plt.show()


# In[ ]:



