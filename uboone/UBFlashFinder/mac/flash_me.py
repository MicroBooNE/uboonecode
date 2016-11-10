import ROOT
from ROOT import larlite
import sys

mgr = larlite.storage_manager()
mgr.set_io_mode(mgr.kREAD)

mgr.add_in_filename(sys.argv[1])

mgr.open()

b = []
    
while mgr.next_event() :

    a = mgr.get_data(larlite.data.kOpFlash,"SimpleFlashFinder")

    if a.size() <= 0:
        continue
    
    for aa in xrange(a.size()):
        b.append(a[aa].TotalPE())
