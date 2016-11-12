import sys
from ROOT import gSystem
gSystem.Load("libOpticalRecoTool_SimpleFlashFinder")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing SimpleFlashFinder..."

sys.exit(0)

