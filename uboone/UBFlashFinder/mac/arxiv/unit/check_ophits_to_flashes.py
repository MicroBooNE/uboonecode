import numpy as np

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib.dates as dts

import matplotlib.patches as patches

import numpy as np
import ROOT  as rt
from larlite import larlite as fmwk


event_ophit = rt.larlite.event_ophit()

channel = int(np.floor(32*np.random.random()))
readout_window = 23.4

t_starts = []
ophits   = []    

#add some random flashes...
for i in xrange(0,10): #10 flashes
    t_peak = 23.4*np.random.random_sample()
    width = 10.0
    area = 100*np.random.random_sample()
    peakheight = np.random.uniform(1,15)
    pe = peakheight/2.0
    fastototal = 0.0;
        
    hit = rt.larlite.ophit(channel,
                           t_peak,
                           0,
                           0,
                           width,
                           area,
                           peakheight,
                           pe,
                           fastototal)
    ophits.append(hit)

#add a train of pulses
for i in xrange(0,20):
    t_peak = 10 + i*np.random.uniform(0,0.1)
    width = 3
    area = 100;
    peakheight = 30 - float(i)/2.0 
    pe = 30 - float(i)/2.0 
    fastototal = 0.0;
        
    hit = rt.larlite.ophit(channel,
                           t_peak,
                           0,
                           0,
                           width,
                           area,
                           peakheight,
                           pe,
                           fastototal)
    ophits.append(hit)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[52]:

mgr = fmwk.storage_manager()
mgr.set_io_mode(fmwk.storage_manager.kWRITE)
mgr.set_out_filename("temp.root")
mgr.open()
mgr.set_id(1,0,1)
mgr.event_id()

event_ophit = mgr.get_data(fmwk.data.kOpHit,'OpHitFinder')

for h in ophits:
    event_ophit.push_back(h)

mgr.next_event()
mgr.close()


# In[53]:

my_proc = fmwk.ana_processor()
my_proc.add_input_file ("temp.root")
my_proc.set_output_file("temp2.root")
my_proc.set_io_mode(fmwk.storage_manager.kBOTH) 
#my_proc.set_data_to_write(fmwk.data.kOpFlash,'funny')
my_module = fmwk.SimpleFlashFinder()
my_module.Configure("../simpleflashfindermodule.fcl")
my_proc.add_process(my_module)
my_proc.run()

mgr = fmwk.storage_manager()
mgr.set_io_mode(mgr.kREAD)
mgr.add_in_filename("temp2.root")
mgr.open()
mgr.next_event()
event_opflash = mgr.get_data(fmwk.data.kOpFlash,'SimpleFlashFinder')



plt.rcParams.update({'font.size': 16}) 
fig,ax = plt.subplots(figsize=(15,6))

for i in np.linspace(0.1,23.4,234):
    ax.vlines(i,0,100,lw=1,color='red',linestyles='dashed',alpha=0.7)

for h in ophits:
    ax.vlines(h.PeakTime(),0,h.PE(),lw=4,color='black',alpha=1)
    
    
for i in xrange(event_opflash.size()):
    flash = event_opflash[i]
    ax.vlines(flash.Time(),0,flash.TotalPE(),lw=2,color='orange',alpha=0.8)

plt.xlim(0,24)
plt.xlabel('Time [us]')
plt.ylabel('PE')
plt.show()

import os
os.system('rm -r temp.root')
os.system('rm -r temp2.root')


